#include <JuceHeader.h>
#include <cmath>
#include <iostream>

const double twoPI = 2 * acos(-1);
const double SAMPLE_RATE = 48000.0;
const double FREQ1 = 1000.0;
const double FREQ10 = 10000.0;
const double AMPLI = 0.5;

class audioCallback final : public juce::AudioIODeviceCallback
{
public:
    void audioDeviceAboutToStart(juce::AudioIODevice* device) override 
    {
        sampleRate = device->getCurrentSampleRate();
        phase1 = 0.0;
        phase10 = 0.0;
    }

    void audioDeviceStopped() override {}

    void audioDeviceIOCallbackWithContext(
        const float* const* inputChannelData,
        int numInputChannels,
        float* const* outputChannelData,
        int numOutputChannels,
        int numSamples,
        const juce::AudioIODeviceCallbackContext& context) override
    {
        const double phaseIncre1 = twoPI * FREQ1 / sampleRate;
        const double phaseIncre10 = twoPI * FREQ10 / sampleRate;
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++)
        {
            const double tone1 = AMPLI * std::sin(phase1);
            const double tone10 = AMPLI * std::sin(phase10);
            const double sample = tone1 + tone10;

            for (int channel = 0; channel < numOutputChannels; channel++)
            {
                outputChannelData[channel][sampleIndex] = static_cast<float>(sample);
            }

            phase1 += phaseIncre1;
            phase10 += phaseIncre10;
            if (phase1 >= twoPI)
            {
                phase1 -= twoPI;
            }
            if (phase10 >= twoPI)
            {
                phase10 -= twoPI;
            }
        }
    }
private:
    double sampleRate = 0.0;
    double phase1 = 0.0;
    double phase10 = 0.0;
};

int main (int argc, char* argv[])
{
    juce::AudioDeviceManager deviceManager;
    deviceManager.initialiseWithDefaultDevices(0, 1);
    auto setup = deviceManager.getAudioDeviceSetup();
    setup.sampleRate = SAMPLE_RATE;
    deviceManager.setAudioDeviceSetup(setup, false);

    audioCallback callback;
    deviceManager.addAudioCallback(&callback);

    std::cout << "Audio callback is running.\n";
    std::cout << "Press ENTER to stop.\n";
    std::cin.get();

    deviceManager.removeAudioCallback(&callback);

    return 0;
}
