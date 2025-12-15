#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/bundles/ata/mode_b.hh"

#include "hpguppi_blade_runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_b_capi.h"
}

using namespace Blade;
using namespace Blade::Bundles::ATA;

using BladePipeline = ModeB<CI8, BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>;

template<typename IT, typename OT>
class ModeBRunner : public ModeRunner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeBRunner(
        const Config& config,
        struct blade_ata_mode_b_config* ata_b_config,
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
        std::vector<RA_DEC> beamCoordinates(ata_b_config->beamformerBeams);
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
        for(i = 0; i < ata_b_config->beamformerBeams; i++){
            beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
            beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
        }
        memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
                antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

        auto phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
            config.inputShape.numberOfAspects(),
            config.inputShape.numberOfFrequencyChannels() * ata_b_config->channelizerRate,
            1,
            config.inputShape.numberOfPolarizations(),
        });

        const size_t calAntStride = 1;
        const size_t calPolStride = config.inputShape.numberOfAspects() * calAntStride;
        const size_t calChnStride = config.inputShape.numberOfPolarizations() * calPolStride;

        const size_t weightsPolStride = 1;
        const size_t weightsChnStride = config.inputShape.numberOfPolarizations() * weightsPolStride;
        const size_t weightsAntStride = config.inputShape.numberOfFrequencyChannels() * ata_b_config->channelizerRate * weightsChnStride;
        BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", config.inputShape.numberOfFrequencyChannels(), ata_b_config->channelizerRate);

        U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
        for (antIdx = 0; antIdx < config.inputShape.numberOfAspects(); antIdx++) {
            for (chnIdx = 0; chnIdx < config.inputShape.numberOfFrequencyChannels(); chnIdx++) {
                for (polIdx = 0; polIdx < config.inputShape.numberOfPolarizations(); polIdx++) {
                    inputIdx = chnIdx * calChnStride +
                        polIdx * calPolStride + 
                        antIdx * calAntStride;
                    for (fchIdx = 0; fchIdx < ata_b_config->channelizerRate; fchIdx++) {
                        frqIdx = chnIdx * ata_b_config->channelizerRate + fchIdx;
                        outputIdx = antIdx * weightsAntStride +
                            polIdx * weightsPolStride +
                            frqIdx * weightsChnStride;

                        phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                    }
                }
            }
        }

        this->inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(config.inputShape);

        BladePipeline::Config pipeline_config = {      
            .inputShape = config.inputShape,
            .outputShape = config.outputShape,

            .preBeamformerChannelizerRate = ata_b_config->channelizerRate,
            // .preBeamformerPolarizerConvertToCircular = BLADE_ATA_MODE_B_CIRCULAR_POLARIZATION,

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

            .casterBlockSize = ata_b_config->castBlockSize,
            .channelizerBlockSize = ata_b_config->channelizerBlockSize,
            .beamformerBlockSize = ata_b_config->beamformerBlockSize,
            .detectorBlockSize = ata_b_config->beamformerBlockSize
        };
        this->connect(
            this->pipeline,
            pipeline_config,
            {
                .dut = blockDut1,
                .julianDate = blockJulianDate,
                .buffer = this->inputBuffer,
            }
        );
    }


    Result transferIn(const ArrayTensor<Device::CPU, IT>& cpuInputBuffer) override {
        BL_CHECK(this->copy(inputBuffer, cpuInputBuffer));
        return Result::SUCCESS;
    }

    Result transferResult() override {
        BL_CHECK(this->copy(outputBuffer, pipeline->getOutputBuffer()));
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
    std::shared_ptr<BladePipeline> pipeline;
    
    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
    Tensor<Device::CPU, F64> blockJulianDate;
    Tensor<Device::CPU, F64> blockDut1;
};

bool blade_ata_b_initialize(
    struct blade_ata_mode_b_config ata_b_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    using ModeBRunner = ModeBRunner<CI8, BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>;

    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }
    BL_INFO("Initializing...");

    State.inputShape = ArrayShape({
        N_INPUT_ASPECTS,
        N_INPUT_CHANNELS,
        N_INPUT_BLOCK_TIME,
        N_INPUT_POL,
    });

    State.outputShape = ArrayShape({
        ata_b_config.beamformerBeams,
        N_INPUT_CHANNELS * ata_b_config.channelizerRate,
        N_INPUT_BLOCK_TIME / ata_b_config.channelizerRate,
        2 // detector disabled
    });

    ModeBRunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape
    };
    State.pipelineRunner = std::make_shared<ModeBRunner>(
        config,
        &ata_b_config,
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

void blade_ata_b_set_block_time_mjd(double mjd) {
    using ModeBRunner = ModeBRunner<CI8, BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeBRunner>(State.pipelineRunner)->setJulianDate(mjd);
}

void blade_ata_b_set_block_dut1(double dut1) {
    using ModeBRunner = ModeBRunner<CI8, BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeBRunner>(State.pipelineRunner)->setDut1(dut1);
}
