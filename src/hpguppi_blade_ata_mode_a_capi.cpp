#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/pipelines/ata/mode_b.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_a_capi.h"
}

using namespace Blade;
using namespace Blade::Pipelines::ATA;

using BladePipeline = ModeB<BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;
static std::unique_ptr<Runner<BladePipeline>> runner;

bool blade_ata_a_initialize(
    struct blade_ata_mode_a_config ata_a_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    if (runner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    std::vector<XYZ> antennaPositions(ata_a_config.inputDims.NANTS);
    std::vector<RA_DEC> beamCoordinates(ata_a_config.beamformerBeams);
    std::vector<std::complex<double>> antennaCalibrationsCpp(
            ata_a_config.inputDims.NANTS*\
            ata_a_config.inputDims.NCHANS*\
            ata_a_config.inputDims.NPOLS);
    int i;
    for(i = 0; i < ata_a_config.inputDims.NANTS; i++){
        antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
        antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
        antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
    }
    for(i = 0; i < ata_a_config.beamformerBeams; i++){
        beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
        beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
    }
    memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
            antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

    BladePipeline::Config config = {
        .preBeamformerChannelizerRate = ata_a_config.channelizerRate,

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
        .phasorAntennaCalibrations = antennaCalibrationsCpp,
        .phasorBeamCoordinates = beamCoordinates,
        
        .beamformerNumberOfAntennas  = ata_a_config.inputDims.NANTS,
        .beamformerNumberOfFrequencyChannels = ata_a_config.inputDims.NCHANS,
        .beamformerNumberOfTimeSamples  = ata_a_config.inputDims.NTIME,
        .beamformerNumberOfPolarizations  = ata_a_config.inputDims.NPOLS,
        .beamformerNumberOfBeams = ata_a_config.beamformerBeams,
        .beamformerIncoherentBeam = BLADE_ATA_MODE_A_OUTPUT_INCOHERENT_BEAM,

        .detectorEnable = true,
        .detectorIntegrationSize = ata_a_config.integrationSize,
        .detectorNumberOfOutputPolarizations = ata_a_config.numberOfOutputPolarizations,

        .outputMemWidth = ata_a_config.outputMemWidth,
        .outputMemPad = ata_a_config.outputMemPad,

        .castBlockSize = ata_a_config.castBlockSize,
        .channelizerBlockSize = ata_a_config.channelizerBlockSize,
        .beamformerBlockSize = ata_a_config.beamformerBlockSize,
        .detectorBlockSize = ata_a_config.detectorBlockSize
    };

    
    runner = Runner<BladePipeline>::New(numberOfWorkers, config);

    return true;
}

void blade_ata_a_terminate() {
    if (!runner) {
        BL_FATAL("Can't terminate because Blade Runner isn't initialized.");
        throw Result::ASSERTION_ERROR;
    }
    runner.reset();
}

size_t blade_ata_a_get_input_size() {
    assert(runner);
    return runner->getWorker().getInputSize();
}

size_t blade_ata_a_get_output_size() {
    assert(runner);
    return runner->getWorker().getOutputSize();
}

bool blade_ata_a_enqueue(void* input_ptr, void* output_ptr, size_t id, double time_mjd, double dut1) {
    assert(runner);
    return runner->enqueue([&](auto& worker){
        auto input = Vector<Device::CPU, CI8>(input_ptr, worker.getInputSize());
        auto output = Vector<Device::CPU, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>(output_ptr, worker.getOutputSize());
        auto time_mjd_vec = Vector<Device::CPU, F64>(&time_mjd, 1);
        auto du1_vec = Vector<Device::CPU, F64>(&dut1, 1);

        worker.run(time_mjd_vec, du1_vec, input, output);

        return id;
    });
}

bool blade_ata_a_dequeue(size_t* id) {
    assert(runner);
    return runner->dequeue(id);
}
