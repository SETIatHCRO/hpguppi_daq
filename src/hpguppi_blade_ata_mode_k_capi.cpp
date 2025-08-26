#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/modules/caster.hh"
#include "blade/modules/kurtosis.hh"

#include "hpguppi_blade_runner.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_k_capi.h"
}

using namespace Blade;

template<typename IT, typename OT>
class ModeKRunner : public ModeRunner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeKRunner(
        const Config& config,
        struct blade_ata_mode_k_config* ata_k_config
    )
        : inputBuffer(config.inputShape),
          outputBuffer(config.outputShape)
    {
        BL_DEBUG("Instantiating input caster module.");
        this->connect(
            inputCaster,
            {},
            {
                .buf = inputBuffer
            }
        );
        
        std::string kurtosisMaskOutputFilepath(
            ata_k_config->kurtosisMaskOutputFilepath
        );
        Kurtosis::Config kurtosis_cfg = {
            .debugMode = false,
            .kurtosisChannelLength = (int) ata_k_config->kurtosisChannelLength,
            .numberOfKurtosisStddev = (int) ata_k_config->kurtosisSigma,
            .numberOfMaskRuns = (int) ata_k_config->kurtosisNumberOfMaskRuns,
            .maskFilePath = kurtosisMaskOutputFilepath,
        };

        if (ata_k_config->kurtosisSigma > 5 || ata_k_config->kurtosisSigma < 3) {
            BL_FATAL("Kurtosis sigma ({}) must be in range: [3, 5].", ata_k_config->kurtosisSigma);
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

        this->compile();
    }

    Result transferIn(const ArrayTensor<Device::CPU, IT>& cpuInputBuffer) override {
        BL_CHECK(this->copy(inputBuffer, cpuInputBuffer));
        return Result::SUCCESS;
    }

    Result transferResult() override {
        BL_CHECK(this->copy(outputBuffer, kurtosis->getOutputBuffer()));
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

    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
};

bool blade_ata_k_initialize(
    struct blade_ata_mode_k_config ata_k_config
) {
    using namespace std::complex_literals;
    using ModeKRunner = ModeKRunner<CI8, BLADE_ATA_MODE_K_OUTPUT_ELEMENT_T>;

    BL_INFO("Initializing Mode K...");
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
        N_INPUT_ASPECTS,
        N_INPUT_CHANNELS,
        N_INPUT_BLOCK_TIME,
        N_INPUT_POL,
    });

    ModeKRunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape,
    };
    State.pipelineRunner = std::make_shared<ModeKRunner>(
        config,
        &ata_k_config
    );

    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    return true;
}
