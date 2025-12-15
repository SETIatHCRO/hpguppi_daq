#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/bundles/ata/mode_b.hh"
#include "blade/modules/caster.hh"
#include "blade/modules/kurtosis.hh"

#include "hpguppi_blade_runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_a_capi.h"
}

using namespace Blade;
using namespace Blade::Bundles::ATA;

using BladePipeline = ModeB<CF32, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;

template<typename IT, typename OT>
class ModeARunner : public ModeRunner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeARunner(
        const Config& config,
        struct blade_ata_mode_a_config* ata_a_config,
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
        std::vector<RA_DEC> beamCoordinates(ata_a_config->beamformerBeams);
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
        for(i = 0; i < ata_a_config->beamformerBeams; i++){
            beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
            beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
        }
        memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
                antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

        auto phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
            config.inputShape.numberOfAspects(),
            config.inputShape.numberOfFrequencyChannels() * ata_a_config->channelizerRate,
            1,
            config.inputShape.numberOfPolarizations(),
        });

        const size_t calAntStride = 1;
        const size_t calPolStride = config.inputShape.numberOfAspects() * calAntStride;
        const size_t calChnStride = config.inputShape.numberOfPolarizations() * calPolStride;

        const size_t weightsPolStride = 1;
        const size_t weightsChnStride = config.inputShape.numberOfPolarizations() * weightsPolStride;
        const size_t weightsAntStride = config.inputShape.numberOfFrequencyChannels() * ata_a_config->channelizerRate * weightsChnStride;
        BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", config.inputShape.numberOfFrequencyChannels(), ata_a_config->channelizerRate);

        U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
        for (antIdx = 0; antIdx < config.inputShape.numberOfAspects(); antIdx++) {
            for (chnIdx = 0; chnIdx < config.inputShape.numberOfFrequencyChannels(); chnIdx++) {
                for (polIdx = 0; polIdx < config.inputShape.numberOfPolarizations(); polIdx++) {
                    inputIdx = chnIdx * calChnStride +
                        polIdx * calPolStride + 
                        antIdx * calAntStride;
                    for (fchIdx = 0; fchIdx < ata_a_config->channelizerRate; fchIdx++) {
                        frqIdx = chnIdx * ata_a_config->channelizerRate + fchIdx;
                        outputIdx = antIdx * weightsAntStride +
                            polIdx * weightsPolStride +
                            frqIdx * weightsChnStride;

                        phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                    }
                }
            }
        }

        this->inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(config.inputShape);
        this->connect(
            this->inputCaster,
            {},
            {
                .buf = this->inputBuffer,
            }
        );
        if (ata_a_config->kurtosisEnabled) {
            std::string kurtosisMaskOutputFilepath(
                ata_a_config->kurtosisMaskOutputFilepath
            );
            Kurtosis::Config kurtosis_cfg = {
                .debugMode = false,
                .kurtosisChannelLength = (int) ata_a_config->kurtosisChannelLength,
                .numberOfKurtosisStddev = (int) ata_a_config->kurtosisSigma,
                .numberOfMaskRuns = (int) ata_a_config->kurtosisNumberOfMaskRuns,
                .maskFilePath = kurtosisMaskOutputFilepath,
            };
            if (ata_a_config->kurtosisSigma > 5 || ata_a_config->kurtosisSigma < 3) {
                BL_FATAL("Kurtosis sigma ({}) must be in range: [3, 5].", ata_a_config->kurtosisSigma);
                throw Result::ASSERTION_ERROR;
            }
            BL_DEBUG("Instantiating Kurtosis module.");
            this->connect(
                this->kurtosis,
                kurtosis_cfg,
                {
                    .buf = this->inputCaster->getOutputBuffer(),
                }
            );
        }
        else {
            BL_DEBUG("Bypassing Kurtosis module.");
        }
        
        BladePipeline::Config pipeline_config = {
            .inputShape = config.inputShape,
            .outputShape = config.outputShape,

            .preBeamformerChannelizerRate = ata_a_config->channelizerRate,
            // .preBeamformerPolarizerConvertToCircular = BLADE_ATA_MODE_A_CIRCULAR_POLARIZATION,

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

            .beamformerIncoherentBeam = BLADE_ATA_MODE_A_OUTPUT_INCOHERENT_BEAM,

            .detectorEnable = true,
            .detectorIntegrationRate = ata_a_config->integrationSize,
            .detectorNumberOfOutputPolarizations = ata_a_config->numberOfOutputPolarizations,

            .casterBlockSize = ata_a_config->castBlockSize,
            .channelizerBlockSize = ata_a_config->channelizerBlockSize,
            .beamformerBlockSize = ata_a_config->beamformerBlockSize,
            .detectorBlockSize = ata_a_config->detectorBlockSize
        };
        if (ata_a_config->kurtosisEnabled) {
            this->connect(
                this->pipeline,
                pipeline_config,
                {
                    .dut = blockDut1,
                    .julianDate = blockJulianDate,
                    .buffer = this->kurtosis->getOutputBuffer()
                }
            );
        }
        else {
            this->connect(
                this->pipeline,
                pipeline_config,
                {
                    .dut = blockDut1,
                    .julianDate = blockJulianDate,
                    .buffer = this->inputCaster->getOutputBuffer()
                }
            );

        }
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
    // std::shared_ptr<Runner> pipelineRunner;
    std::shared_ptr<BladePipeline> pipeline;
    using InputCaster = typename Modules::Caster<IT, CF32>;
    std::shared_ptr<InputCaster> inputCaster;
    using Kurtosis = typename Modules::Kurtosis<CF32, CF32>;
    std::shared_ptr<Kurtosis> kurtosis;
    
    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
    Tensor<Device::CPU, F64> blockJulianDate;
    Tensor<Device::CPU, F64> blockDut1;
};

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
    using ModeARunner = ModeARunner<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;

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
        ata_a_config.beamformerBeams + (BLADE_ATA_MODE_A_OUTPUT_INCOHERENT_BEAM ? 1 : 0),
        N_INPUT_CHANNELS * ata_a_config.channelizerRate,
        N_INPUT_BLOCK_TIME / (ata_a_config.channelizerRate * ata_a_config.integrationSize),
        ata_a_config.numberOfOutputPolarizations, // detector enabled
    });

    ModeARunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape
    };
    State.pipelineRunner = std::make_shared<ModeARunner>(
        config,
        &ata_a_config,
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

void blade_ata_a_set_block_time_mjd(double mjd) {
    using ModeARunner = ModeARunner<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeARunner>(State.pipelineRunner)->setJulianDate(mjd);
}

void blade_ata_a_set_block_dut1(double dut1) {
    using ModeARunner = ModeARunner<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;
    std::static_pointer_cast<ModeARunner>(State.pipelineRunner)->setDut1(dut1);
}
