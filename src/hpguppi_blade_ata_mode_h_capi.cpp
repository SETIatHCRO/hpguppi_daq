#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/pipelines/ata/mode_b.hh"
#include "blade/pipelines/generic/mode_h.hh"
#include "blade/plan.hh"
#include "blade/runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_h_capi.h"
}

using namespace Blade;
using namespace Blade::Pipelines::ATA;
using namespace Blade::Pipelines::Generic;

using BladePipelineB = ModeB<CF32>;
using BladePipelineH = ModeH<CF32, F32>;

static struct {
    struct {
        BladePipelineB::Config B;
        BladePipelineH::Config H;
    } RunnersConfig;

    struct {
        std::unique_ptr<Runner<BladePipelineB>> B; 
        std::unique_ptr<Runner<BladePipelineH>> H; 
    } RunnersInstances;
} State;


bool blade_ata_h_initialize(
    struct blade_ata_mode_h_config ata_h_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    if (State.RunnersInstances.B || State.RunnersInstances.H) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    std::vector<XYZ> antennaPositions(ata_h_config.inputDims.NANTS);
    std::vector<RA_DEC> beamCoordinates(ata_h_config.beamformerBeams);
    std::vector<std::complex<double>> antennaCalibrationsCpp(
            ata_h_config.inputDims.NANTS*\
            ata_h_config.inputDims.NCHANS*\
            ata_h_config.inputDims.NPOLS);
    int i;
    for(i = 0; i < ata_h_config.inputDims.NANTS; i++){
        antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
        antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
        antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
    }
    for(i = 0; i < ata_h_config.beamformerBeams; i++){
        beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
        beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
    }
    memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
            antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

    State.RunnersConfig.B = {  
        .preBeamformerChannelizerRate = ata_h_config.channelizerRate,

        .phasorObservationFrequencyHz = observationMeta->rfFrequencyHz,
        .phasorChannelBandwidthHz = observationMeta->channelBandwidthHz,
        .phasorTotalBandwidthHz = observationMeta->totalBandwidthHz,
        .phasorFrequencyStartIndex = observationMeta->frequencyStartIndex,
        .phasorReferenceAntennaIndex = observationMeta->referenceAntennaIndex,
        .phasorArrayReferencePosition = {
            .LON = arrayReferencePosition->LON,
            .LAT = arrayReferencePosition->LAT,
            .ALT = arrayReferencePosition->ALT
        },
        .phasorBoresightCoordinate = {
            .RA = obs_phase_center_radecrad[0],
            .DEC = obs_phase_center_radecrad[1]
        },
        .phasorAntennaPositions = antennaPositions,
        .phasorAntennaCalibrations = {0},
        .phasorBeamCoordinates = beamCoordinates,
        
        .beamformerNumberOfAntennas  = ata_h_config.inputDims.NANTS,
        .beamformerNumberOfFrequencyChannels = ata_h_config.inputDims.NCHANS,
        .beamformerNumberOfTimeSamples  = ata_h_config.inputDims.NTIME,
        .beamformerNumberOfPolarizations  = ata_h_config.inputDims.NPOLS,
        .beamformerNumberOfBeams = ata_h_config.beamformerBeams,
        .beamformerIncoherentBeam = false,

        .detectorEnable = false,
        .detectorIntegrationSize = 1,
        .detectorNumberOfOutputPolarizations = 1,

        .outputMemWidth = ata_h_config.outputMemWidth,
        .outputMemPad = ata_h_config.outputMemPad,

        .castBlockSize = ata_h_config.castBlockSize,
        .channelizerBlockSize = ata_h_config.channelizerBlockSize,
        .beamformerBlockSize = ata_h_config.beamformerBlockSize,
    };
    
    
    State.RunnersConfig.B.phasorAntennaCalibrations.resize(
        State.RunnersConfig.B.beamformerNumberOfAntennas *
        State.RunnersConfig.B.beamformerNumberOfFrequencyChannels *
        State.RunnersConfig.B.preBeamformerChannelizerRate *
        State.RunnersConfig.B.beamformerNumberOfPolarizations
    );

    const size_t calAntStride = 1;
    const size_t calPolStride = State.RunnersConfig.B.beamformerNumberOfAntennas * calAntStride;
    const size_t calChnStride = State.RunnersConfig.B.beamformerNumberOfPolarizations * calPolStride;

    const size_t weightsPolStride = 1;
    const size_t weightsChnStride = State.RunnersConfig.B.beamformerNumberOfPolarizations * weightsPolStride;
    const size_t weightsAntStride = State.RunnersConfig.B.beamformerNumberOfFrequencyChannels * State.RunnersConfig.B.preBeamformerChannelizerRate * weightsChnStride;
    BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", State.RunnersConfig.B.beamformerNumberOfFrequencyChannels, State.RunnersConfig.B.preBeamformerChannelizerRate);

    for (U64 antIdx = 0; antIdx < State.RunnersConfig.B.beamformerNumberOfAntennas; antIdx++) {
        for (U64 chnIdx = 0; chnIdx < State.RunnersConfig.B.beamformerNumberOfFrequencyChannels; chnIdx++) {
            for (U64 polIdx = 0; polIdx < State.RunnersConfig.B.beamformerNumberOfPolarizations; polIdx++) {
                for (U64 fchIdx = 0; fchIdx < State.RunnersConfig.B.preBeamformerChannelizerRate; fchIdx++) {
                    const auto inputIdx = chnIdx * calChnStride +
                                          polIdx * calPolStride + 
                                          antIdx * calAntStride;

                    const auto frqIdx = chnIdx * State.RunnersConfig.B.preBeamformerChannelizerRate + fchIdx;
                    const auto outputIdx = antIdx * weightsAntStride +
                                           polIdx * weightsPolStride +
                                           frqIdx * weightsChnStride;

                    State.RunnersConfig.B.phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                }
            }
        }
    }

    State.RunnersConfig.H = {
        .accumulateRate = BLADE_ATA_MODE_H_ACCUMULATE_RATE, 
        
        .channelizerNumberOfBeams = State.RunnersConfig.B.beamformerNumberOfBeams,
        .channelizerNumberOfFrequencyChannels = State.RunnersConfig.B.beamformerNumberOfFrequencyChannels * 
                                     State.RunnersConfig.B.preBeamformerChannelizerRate,
        .channelizerNumberOfTimeSamples = State.RunnersConfig.B.beamformerNumberOfTimeSamples / 
                               State.RunnersConfig.B.preBeamformerChannelizerRate,
        .channelizerNumberOfPolarizations = State.RunnersConfig.B.beamformerNumberOfPolarizations,

        .detectorNumberOfOutputPolarizations = ata_h_config.numberOfOutputPolarizations,

        // .channelizerBlockSize = 512,
        // .detectorBlockSize = 512,
    };

    State.RunnersInstances.B = Runner<BladePipelineB>::New(
        numberOfWorkers, 
        State.RunnersConfig.B
    );

    State.RunnersInstances.H = Runner<BladePipelineH>::New(
        numberOfWorkers, 
        State.RunnersConfig.H
    );

    return true;
}

void blade_ata_h_terminate() {
    if (!State.RunnersInstances.B || !State.RunnersInstances.H) {
        BL_FATAL("Can't terminate because Blade Runner isn't initialized.");
        throw Result::ASSERTION_ERROR;
    }

    State.RunnersInstances.B.reset();
    // State.RunnersInstances.H->applyToAllWorkers([&](auto& worker){
    //     worker.resetAccumulatorSteps();
    //     return Result::SUCCESS;
    // });
    State.RunnersInstances.H.reset();
}

U64 blade_ata_h_get_input_size() {
    assert(State.RunnersInstances.B);
    return State.RunnersInstances.B->getWorker().getInputSize();
}

U64 blade_ata_h_get_output_size() {
    assert(State.RunnersInstances.H);
    return State.RunnersInstances.H->getWorker().getOutputSize();
}

size_t blade_ata_h_accumulator_counter() {
    return State.RunnersInstances.H->getNextWorker().getCurrentAccumulatorStep();
}

bool blade_ata_h_enqueue_b(void* input_ptr, const U64 b_id, double time_mjd, double dut1) {
    assert(State.RunnersInstances.B);
    assert(State.RunnersInstances.H);

    if (!State.RunnersInstances.H->slotAvailable()) {
        return false;
    }

    return State.RunnersInstances.B->enqueue([&](auto& worker){
        auto input = Vector<Device::CPU, CI8>(input_ptr, worker.getInputSize());
        auto& next = State.RunnersInstances.H->getNextWorker();
        auto time_mjd_vec = Vector<Device::CPU, F64>(&time_mjd, 1);
        auto dut1_vec = Vector<Device::CPU, F64>(&dut1, 1);

        BL_CHECK_THROW(worker.run(time_mjd_vec, dut1_vec, input, next));

        return b_id;
    });
}

bool blade_ata_h_dequeue_b(U64* b_id) {
    assert(State.RunnersInstances.B);
    return State.RunnersInstances.B->dequeue(b_id);
}

bool blade_ata_h_enqueue_h(void* output_ptr, const U64 h_id) {
    assert(State.RunnersInstances.B);
    assert(State.RunnersInstances.H);

    return State.RunnersInstances.H->enqueue([&](auto& worker) {
        auto output = Vector<Device::CPU, F32>(output_ptr, worker.getOutputSize());

        BL_CHECK_THROW(worker.run(output));

        return h_id;
    });
}

bool blade_ata_h_dequeue_h(U64* id) {
    assert(State.RunnersInstances.H);
    return State.RunnersInstances.H->dequeue(id);
}
