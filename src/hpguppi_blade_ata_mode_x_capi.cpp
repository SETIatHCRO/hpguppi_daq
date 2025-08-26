#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/bundles/generic/mode_x.hh"
#include "blade/modules/kurtosis.hh"

#include "hpguppi_blade_runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_x_capi.h"
}

using namespace Blade;

template<typename IT, typename OT>
class ModeXRunner : public ModeRunner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeXRunner(
        const Config& config,
        U64 channelizationRate,
        U64 integrationRate,
        U64 frequencyIntegrationRate,
        bool kurtosisEnabled,
        int kurtosisSigma,
        std::string kurtosisMaskOutputFilepath
    )
        : inputBuffer(config.inputShape),
          outputBuffer(config.outputShape)
    {
        if (channelizationRate > 1 ) {
            // Spectral line mode (16 MHz firmware)
            if (channelizationRate < config.inputShape.numberOfTimeSamples()) {
                BL_FATAL("Channelizer rate must be a multiple of NTIME: {} < {}.", channelizationRate, config.inputShape.numberOfTimeSamples());
                throw Result::ASSERTION_ERROR;
            }
            if (channelizationRate%config.inputShape.numberOfTimeSamples() != 0) {
                BL_FATAL("Channelizer rate must be a multiple of NTIME: {}%{} != 0.", channelizationRate, config.inputShape.numberOfTimeSamples());
                throw Result::ASSERTION_ERROR;
            }
        }

        bool contiuumNotSpectral = channelizationRate == 1;
        // modeX will always produce T=1 output shapes:
        //   if channelizer is not bypassed, modeX will upchannelize all time first
        //   if channelizer is bypassed, modeX will integrate all input time
        // further integration happens thereafter to reach integrationRate
        U64 preCorrelatorStackerRate = contiuumNotSpectral ? 1 : channelizationRate/config.inputShape.numberOfTimeSamples();
        U64 correlatorIntegrationRate = contiuumNotSpectral ? integrationRate/config.inputShape.numberOfTimeSamples() : integrationRate;
        if (correlatorIntegrationRate % BLADE_ATA_MODE_X_INTEGRATION_FACTOR != 0) {
            BL_FATAL("Correlator integration rate must be a multiple of INTEGRATION_FACTOR: {}%{} != 0.", correlatorIntegrationRate, BLADE_ATA_MODE_X_INTEGRATION_FACTOR);
            throw Result::ASSERTION_ERROR;
        }

        if (kurtosisEnabled) {
            BL_DEBUG("Instantiating input caster module.");
            this->connect(
                inputCaster,
                {},
                {
                    .buf = inputBuffer
                }
            );
            Kurtosis::Config kurtosis_cfg = {
                .debugMode = false,
                .kurtosisChannelLength = 256,
                .numberOfKurtosisStddev = kurtosisSigma,
                .numberOfMaskRuns = 128,
                .maskFilePath = kurtosisMaskOutputFilepath,
            };

            if (kurtosisSigma > 5 || kurtosisSigma < 3) {
                BL_FATAL("Kurtosis sigma ({}) must be in range: [3, 5].", kurtosisSigma);
                throw Result::ASSERTION_ERROR;
            }
            BL_DEBUG("Instantiating Kurtosis module.");
            this->connect(
                kurtosis,
                kurtosis_cfg,
                {
                    .buf = inputCaster->getOutputBuffer()
                }
            );
            
            BL_DEBUG("Instantiating ModeX bundle.");
            this->connect(
                kurtosisPipeline,
                {
                    .inputShape = config.inputShape,
                    .outputShape = config.outputShape,
                    .preChannelizerStackerMultiplier = 1,
                    .channelizerBypass = contiuumNotSpectral,
                    
                    .preCorrelatorStackerMultiplier = BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
                    .correlatorIntegrationRate = correlatorIntegrationRate/BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
                    .correlatorConjugateAntennaIndex = BLADE_ATA_MODE_X_CONJUGATION_INDEX,

                    .correlatorUseSharedMemory = false, //contiuumNotSpectral,
                    .correlatorCalculationMode = CALC_MODE::INTEGER, // contiuumNotSpectral ? CALC_MODE::INTEGER : CALC_MODE::DOUBLE_PRECISION_FP,
                    
                    .postCorrelatorFrequencyIntegrationRate = frequencyIntegrationRate, // TODO support repeated integrations...

                    .correlatorBlockSize = contiuumNotSpectral ? (U64) 64 : (U64) 32
                },
                {
                    // Kurtosis is an inplace operation: its output buffer is its input buffer
                    .buffer = kurtosis->getOutputBuffer()
                }
            );
        }
        else {
            BL_DEBUG("Instantiating ModeX bundle.");
            this->connect(
                pipeline,
                {
                    .inputShape = config.inputShape,
                    .outputShape = config.outputShape,
                    .preChannelizerStackerMultiplier = 1,
                    .channelizerBypass = contiuumNotSpectral,
                    
                    .preCorrelatorStackerMultiplier = BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
                    .correlatorIntegrationRate = correlatorIntegrationRate/BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
                    .correlatorConjugateAntennaIndex = BLADE_ATA_MODE_X_CONJUGATION_INDEX,

                    .correlatorUseSharedMemory = contiuumNotSpectral,
                    .correlatorCalculationMode = CALC_MODE::INTEGER, // contiuumNotSpectral ? CALC_MODE::INTEGER : CALC_MODE::DOUBLE_PRECISION_FP,
                    
                    .postCorrelatorFrequencyIntegrationRate = frequencyIntegrationRate, // TODO support repeated integrations...

                    .correlatorBlockSize = contiuumNotSpectral ? (U64) 64 : (U64) 32
                },
                {
                    // Kurtosis is an inplace operation: its output buffer is its input buffer
                    .buffer = inputBuffer
                }
            );
        }

        this->compile();
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

 private:
 
    using InputCaster = typename Modules::Caster<IT, CF32>;
    std::shared_ptr<InputCaster> inputCaster;
    using Kurtosis = typename Modules::Kurtosis<CF32, CF32>;
    std::shared_ptr<Kurtosis> kurtosis;
    using KurtosisModeX = Bundles::Generic::ModeX<CF32, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>;
    std::shared_ptr<KurtosisModeX> kurtosisPipeline;

    using ModeX = Bundles::Generic::ModeX<IT, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>;
    std::shared_ptr<ModeX> pipeline;

    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
};

bool blade_ata_x_initialize(
    struct blade_ata_mode_x_config ata_x_config
) {
    using namespace std::complex_literals;
    using ModeXRunner = ModeXRunner<CI8, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>;

    BL_INFO("Initializing Mode X...");
    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    State.inputShape = ArrayShape({
        N_INPUT_ASPECTS,
        N_INPUT_CHANNELS,
        N_INPUT_BLOCK_TIME,
        N_INPUT_POL,
    });

    State.outputShape = ArrayShape({
        N_INPUT_ASPECTS*(N_INPUT_ASPECTS+1)/2,
        N_INPUT_CHANNELS*ata_x_config.channelizerRate/ata_x_config.frequencyIntegrationSize,
        1,
        N_INPUT_POL * N_INPUT_POL,
    });

    ModeXRunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape,
    };
    std::string kurtosisMaskOutputFilepath(
        ata_x_config.kurtosisMaskOutputFilepath, sizeof(ata_x_config.kurtosisMaskOutputFilepath)
    );
    State.pipelineRunner = std::make_shared<ModeXRunner>(
        config,
        ata_x_config.channelizerRate,
        ata_x_config.integrationSize,
        ata_x_config.frequencyIntegrationSize,
        ata_x_config.kurtosisEnabled,
        ata_x_config.kurtosisSigma,
        kurtosisMaskOutputFilepath
    );

    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    // State.debugInput = ArrayTensor<Device::CPU, CI8>(State.inputShape);
    // size_t index = 0;
    // for (int a = 0; a < N_INPUT_ASPECTS; a++) {
    //     for (int c = 0; c < N_INPUT_CHANNELS; c++) {
    //         for (int t = 0; t < N_INPUT_BLOCK_TIME; t++) {
    //             for (int p = 0; p < ata_x_config.inputDims.NPOLS; p++) {
    //                 if (a == 2) {
    //                     if (p==0)
    //                         State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                     else
    //                         State.debugInput[index++] = std::complex<int8_t>(0, 0);
    //                 } else if (a == 3) {
    //                     if (p==0)
    //                         State.debugInput[index++] = std::complex<int8_t>(0, 0);
    //                     else
    //                         State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                 } else if (a == 1) {
    //                     State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                 } else {
    //                     State.debugInput[index++] = std::complex<int8_t>(a*3, (t%250-125)+p+1);
    //                 }
    //             }
    //         }
    //     }
    // }
    return true;
}
