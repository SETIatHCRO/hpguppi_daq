#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/logger.hh"
#include "blade/runner.hh"
#include "blade/bundles/ata/mode_b.hh"
#include "blade/bundles/generic/mode_h.hh"
#include "blade/modules/stacker.hh"

#include "hpguppi_blade_runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_h_capi.h"
}

using namespace Blade;
using namespace Blade::Bundles::ATA;
using namespace Blade::Bundles::Generic;

using BladePipelineB = ModeB<CI8, CF32>;
using BladePipelineH = ModeH<CF32, F32>;

template<typename IT, typename OT>
class ModeHRunner : public ModeRunner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeHRunner(
        const Config& config,
        struct blade_ata_mode_h_config* ata_h_config,
        struct blade_ata_observation_meta* observationMeta,
        struct LonLatAlt* arrayReferencePosition,
        double* obs_phase_center_radecrad,
        double* beamCoordinates_radecrad,
        double* antennaPositions_xyz,
        double _Complex* antennaCalibrations
    )
        : inputBuffer(config.inputShape),
          outputBuffer(config.outputShape),
          blockJulianDate({1}),
          blockDut1({1})
    {
        std::vector<XYZ> antennaPositions(config.inputShape.numberOfAspects());
        std::vector<RA_DEC> beamCoordinates(ata_h_config->beamformerBeams);
        std::vector<std::complex<double>> antennaCalibrationsCpp(
                config.inputShape.numberOfAspects()*\
                config.inputShape.numberOfFrequencyChannels()*\
                config.inputShape.numberOfPolarizations());
        int i;
        for(i = 0; i < config.inputShape.numberOfAspects(); i++){
            antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
            antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
            antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
        }
        for(i = 0; i < ata_h_config->beamformerBeams; i++){
            beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
            beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
        }
        memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
                antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

        auto beamformerOutputShape = ArrayShape({
            ata_h_config->beamformerBeams,
            config.inputShape.numberOfFrequencyChannels(),
            config.inputShape.numberOfTimeSamples(),
            config.inputShape.numberOfPolarizations(),
        });

        auto phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
            config.inputShape.numberOfAspects(),
            config.inputShape.numberOfFrequencyChannels() * ata_h_config->channelizerRate,
            1,
            config.inputShape.numberOfPolarizations(),
        });

        const size_t calAntStride = 1;
        const size_t calPolStride = config.inputShape.numberOfAspects() * calAntStride;
        const size_t calChnStride = config.inputShape.numberOfPolarizations() * calPolStride;

        const size_t weightsPolStride = 1;
        const size_t weightsChnStride = config.inputShape.numberOfPolarizations() * weightsPolStride;
        const size_t weightsAntStride = config.inputShape.numberOfFrequencyChannels() * ata_h_config->channelizerRate * weightsChnStride;
        BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", config.inputShape.numberOfFrequencyChannels(), ata_h_config->channelizerRate);

        U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
        for (antIdx = 0; antIdx < config.inputShape.numberOfAspects(); antIdx++) {
            for (chnIdx = 0; chnIdx < config.inputShape.numberOfFrequencyChannels(); chnIdx++) {
                for (polIdx = 0; polIdx < config.inputShape.numberOfPolarizations(); polIdx++) {
                    inputIdx = chnIdx * calChnStride +
                        polIdx * calPolStride + 
                        antIdx * calAntStride;
                    for (fchIdx = 0; fchIdx < ata_h_config->channelizerRate; fchIdx++) {
                        frqIdx = chnIdx * ata_h_config->channelizerRate + fchIdx;
                        outputIdx = antIdx * weightsAntStride +
                            polIdx * weightsPolStride +
                            frqIdx * weightsChnStride;

                        phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                    }
                }
            }
        }

        this->inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(config.inputShape);

        BladePipelineB::Config configB = {  
            .inputShape = config.inputShape,
            .outputShape = beamformerOutputShape,

            .preBeamformerChannelizerRate = ata_h_config->channelizerRate,
            // .preBeamformerPolarizerConvertToCircular = BLADE_ATA_MODE_H_CIRCULAR_POLARIZATION,

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
            .phasorAntennaCalibrations = phasorAntennaCalibrations,
            .phasorBeamCoordinates = beamCoordinates,

            .beamformerIncoherentBeam = false,

            .detectorEnable = false,
            .detectorIntegrationRate = 1,
            .detectorNumberOfOutputPolarizations = 1,

            .casterBlockSize = ata_h_config->castBlockSize,
            .channelizerBlockSize = ata_h_config->channelizerBlockSize,
            .phasorBlockSize = 512,
            .beamformerBlockSize = ata_h_config->beamformerBlockSize,
            .detectorBlockSize = 512,
        };
        this->connect(
            this->pipelineB,
            configB,
            {
                .dut = blockDut1,
                .julianDate = blockJulianDate,
                .buffer = this->inputBuffer,
            }
        );
        
        this->connect(
            this->stackerFrequency,
            {
                .axis = 2,
                .multiplier = ata_h_config->accumulateRate,
            },
            {
                .buf = this->pipelineB->getOutputBuffer(),
            }
        );

        BL_INFO("Stacker output shape {}.", this->stackerFrequency->getOutputBuffer().shape());

        BladePipelineH::Config configH = {
            .inputShape = this->stackerFrequency->getOutputBuffer().shape(),
            .outputShape = config.outputShape,

            .detectorIntegrationRate = ata_h_config->integrationSize,
            .detectorNumberOfOutputPolarizations = BLADE_ATA_MODE_H_OUTPUT_NPOL,
        };
        this->connect(
            this->pipelineH,
            configH,
            {
                .buffer = this->stackerFrequency->getOutputBuffer(),
            }
        );
    }


    Result transferIn(const ArrayTensor<Device::CPU, IT>& cpuInputBuffer) override {
        BL_CHECK(this->copy(inputBuffer, cpuInputBuffer));
        return Result::SUCCESS;
    }

    Result transferResult() override {
        BL_CHECK(this->copy(outputBuffer, pipelineH->getOutputBuffer()));
        return Result::SUCCESS;
    }
    Result transferOut(void* cpuOutputBuffer) override {
        auto output = ArrayTensor<Device::CPU, OT>(cpuOutputBuffer, State.outputShape);
        BL_CHECK(this->copy(output, outputBuffer));
        return Result::SUCCESS;
    }

    size_t outputByteSize() override {
        return this->outputBuffer.at(0).shape().size()*sizeof(OT);
    }

    void setJulianDate(F64 value) {
        this->blockJulianDate.data()[0] = value;
    }

    void setDut1(F64 value) {
        this->blockDut1.data()[0] = value;
    }

 private:
    std::shared_ptr<Runner> pipelineRunner;
    std::shared_ptr<BladePipelineB> pipelineB;
    std::shared_ptr<Modules::Stacker<CF32, CF32>> stackerFrequency;
    std::shared_ptr<BladePipelineH> pipelineH;
    
    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
    Tensor<Device::CPU, F64> blockJulianDate;
    Tensor<Device::CPU, F64> blockDut1;
};

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
    using ModeHRunner = ModeHRunner<CI8, BLADE_ATA_MODE_H_OUTPUT_ELEMENT_T>;

    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    State.outputShape = ArrayShape({
        ata_h_config.beamformerBeams,
        State.inputShape.numberOfFrequencyChannels() * ata_h_config.channelizerRate * State.inputShape.numberOfTimeSamples() * ata_h_config.accumulateRate,
        1,
        BLADE_ATA_MODE_H_OUTPUT_NPOL,
    });

    ModeHRunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape
    };
    State.pipelineRunner = std::make_shared<ModeHRunner>(
        config,
        &ata_h_config,
        observationMeta,
        arrayReferencePosition,
        obs_phase_center_radecrad,
        beamCoordinates_radecrad,
        antennaPositions_xyz,
        antennaCalibrations
    );
    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    return true;
}

void blade_ata_h_set_block_time_mjd(double mjd) {
    using ModeHRunner = ModeHRunner<CI8, BLADE_ATA_MODE_H_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeHRunner>(State.pipelineRunner)->setJulianDate(mjd);
}

void blade_ata_h_set_block_dut1(double dut1) {
    using ModeHRunner = ModeHRunner<CI8, BLADE_ATA_MODE_H_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeHRunner>(State.pipelineRunner)->setDut1(dut1);
}

size_t blade_ata_h_accumulator_counter() {
    using ModeHRunner = ModeHRunner<CI8, BLADE_ATA_MODE_H_OUTPUT_ELEMENT_T>;
    return std::static_pointer_cast<ModeHRunner>(State.pipelineRunner)->computeCurrentStepCount();
}
