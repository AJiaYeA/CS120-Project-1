#include <JuceHeader.h>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <thread>

const double TWOPI = 2 * acos(-1);

const std::size_t EXPECTEDBITCOUNT = 10000;
const std::size_t BITSPERFRAME = 100;
const std::size_t EXPECTEDFRAMECOUNT = EXPECTEDBITCOUNT / BITSPERFRAME;

const std::size_t CRCBITCOUNT = 8;

const std::size_t ENCODEDBITSPERFRAME = BITSPERFRAME + CRCBITCOUNT;
const double SAMPLERATE = 48000.0;

const std::size_t PRESAMPLECOUNT = 480;
const std::size_t PREHALFSAMPLECOUNT = PRESAMPLECOUNT / 2;
const double PRESTARTFREQ = 2000.0;
const double PREPEAKFREQ = 10000.0;
const float PREAMPLI = 0.5f;

const double CARRIERFREQ = 10000.0;
const std::size_t SAMPLESPERBIT = 48;
const float BPSKAMPLI = 0.5f;
const std::size_t MODULATEDSAMPLESPERFRAME = ENCODEDBITSPERFRAME * SAMPLESPERBIT;

const std::size_t GAPSAMPLECOUNT = 96;
const std::size_t AUDIOSAMPLESPERFRAME = PRESAMPLECOUNT + MODULATEDSAMPLESPERFRAME + GAPSAMPLECOUNT;
const std::size_t TOTALSAMPLES = EXPECTEDFRAMECOUNT * AUDIOSAMPLESPERFRAME;

const std::array<std::uint8_t, CRCBITCOUNT + 1> CRCPOLY
{
    1, 1, 0, 1, 0, 0, 1, 1, 1
};

struct PayloadFrame
{
    std::array<std::uint8_t, BITSPERFRAME> payload{};
    std::array<std::uint8_t, CRCBITCOUNT> crc{};
};

struct EncodedFrame
{
    std::array<std::uint8_t, ENCODEDBITSPERFRAME> bits{};
};

struct ModulatedFrame
{
    std::array<float, MODULATEDSAMPLESPERFRAME> samples{};
};

bool loadBitsFromFile(const std::filesystem::path& file_path, std::vector<std::uint8_t>& bits)
{
    std::ifstream input(file_path, std::ios::binary);
    if (!input.is_open())
    {
        std::cerr << "Cannot open input file:\n" << std::filesystem::absolute(file_path) << '\n';
        return false;
    }

    char bit{};
    for (std::size_t byte_pos = 1; input.get(bit); byte_pos++)
    {
        switch (bit)
        {
        case '0':
            bits.push_back(0);
            break;
        case '1':
            bits.push_back(1);
            break;
        }
    }

    return true;
}

void splitIntoFrames(const std::vector<std::uint8_t>& bits, std::vector<PayloadFrame>& frames)
{
    const auto bit_count = bits.size();

    const std::size_t frame_count = bit_count / BITSPERFRAME;

    frames.clear();
    frames.resize(frame_count);

    for (std::size_t frame_index = 0; frame_index < frame_count; frame_index++)
    {
        std::copy_n(
            bits.begin() + static_cast<std::ptrdiff_t>(frame_index * BITSPERFRAME), 
            BITSPERFRAME, 
            frames[frame_index].payload.begin()
        );
    }
}

template <std::size_t data_bit_count> 
std::array<std::uint8_t, CRCBITCOUNT> calculateCrc(const std::array<std::uint8_t, data_bit_count>& data)
{
    std::array<std::uint8_t, data_bit_count + CRCBITCOUNT> working_bits{};
    std::copy(data.begin(), data.end(), working_bits.begin());

    for (std::size_t start = 0; start < data_bit_count; start++)
    {
        if (working_bits[start] == 0)
        {
            continue;
        }

        for (std::size_t poly_index = 0; poly_index < CRCPOLY.size(); poly_index++)
        {
            working_bits[start + poly_index] ^= CRCPOLY[poly_index];
        }
    }

    std::array<std::uint8_t, CRCBITCOUNT> remainder{};
    std::copy_n(
        working_bits.begin() + static_cast<std::ptrdiff_t>(data_bit_count), 
        CRCBITCOUNT, 
        remainder.begin()
    );

    return remainder;
}

void calculateFrameCrcs(std::vector<PayloadFrame>& frames)
{
    for (auto& frame : frames)
    {
        frame.crc = calculateCrc(frame.payload);
    }
}

std::vector<EncodedFrame> buildEncodedFrames(const std::vector<PayloadFrame>& payloadFrames)
{
    const auto frame_count = payloadFrames.size();

    std::vector<EncodedFrame> encoded_frames;
    encoded_frames.resize(frame_count);

    for (std::size_t frame_index = 0; frame_index < frame_count; frame_index++)
    {
        std::copy(
            payloadFrames[frame_index].payload.begin(), 
            payloadFrames[frame_index].payload.end(), 
            encoded_frames[frame_index].bits.begin()
        );
        std::copy(
            payloadFrames[frame_index].crc.begin(), 
            payloadFrames[frame_index].crc.end(), 
            encoded_frames[frame_index].bits.begin() + static_cast<ptrdiff_t>(BITSPERFRAME)
        );
    }

    return encoded_frames;
}

double getPreambleFrequency(std::size_t sample_index)
{
    if (sample_index < PREHALFSAMPLECOUNT)
    {
        const double ratio = static_cast<double>(sample_index) / static_cast<double>(PREHALFSAMPLECOUNT - 1);
        return PRESTARTFREQ + ratio * (PREPEAKFREQ - PRESTARTFREQ);
    }
    else
    {
        sample_index -= PREHALFSAMPLECOUNT;
        const double ratio = static_cast<double>(sample_index) / static_cast<double>(PREHALFSAMPLECOUNT - 1);
        return PREPEAKFREQ + ratio * (PRESTARTFREQ - PREPEAKFREQ);
    }
}

std::array<float, PRESAMPLECOUNT> generatePreamble()
{
    std::array<float, PRESAMPLECOUNT> preamble{};
    preamble[0] = 0.0f;

    double phase = 0.0;

    for (std::size_t sample_index = 1; sample_index < PRESAMPLECOUNT; sample_index++)
    {
        const auto previous_freq = getPreambleFrequency(sample_index - 1);
        const auto current_freq = getPreambleFrequency(sample_index);
        const double average_freq = 0.5 * (previous_freq + current_freq);

        phase += TWOPI * average_freq / SAMPLERATE;
        if (phase >= TWOPI)
        {
            phase -= TWOPI;
        }

        preamble[sample_index] = PREAMPLI * static_cast<float>(std::sin(phase));
    }

    return preamble;
}

std::array<float, SAMPLESPERBIT> generateCarrierTemplate()
{
    std::array<float, SAMPLESPERBIT> carrier{};

    const double phase_incre = TWOPI * CARRIERFREQ / SAMPLERATE;

    for (std::size_t sample_index = 0; sample_index < SAMPLESPERBIT; sample_index++)
    {
        const double phase = phase_incre * static_cast<double>(sample_index);
        carrier[sample_index] = static_cast<float>(std::sin(phase));
    }

    return carrier;
}

ModulatedFrame modulateFrame(const EncodedFrame& encoded_frame, const std::array<float, SAMPLESPERBIT>& carrier)
{
    ModulatedFrame modulated_frame;

    for (std::size_t bit_index = 0; bit_index < encoded_frame.bits.size(); bit_index++)
    {
        const float symbol = encoded_frame.bits[bit_index] ? 1.0f : -1.0f;

        for (std::size_t sample_index = 0; sample_index < SAMPLESPERBIT; sample_index++)
        {
            modulated_frame.samples[bit_index * SAMPLESPERBIT + sample_index] = BPSKAMPLI * symbol * carrier[sample_index];
        }
    }
    
    return modulated_frame;
}

std::vector<ModulatedFrame> modulateAllFrames(const std::vector<EncodedFrame>& encoded_frames, const std::array<float, SAMPLESPERBIT>& carrier)
{
    std::vector<ModulatedFrame> modulated_frames;
    modulated_frames.reserve(encoded_frames.size());

    for (const auto& frame : encoded_frames)
    {
        modulated_frames.push_back(modulateFrame(frame, carrier));
    }

    return modulated_frames;
}

void buildTransmissionSamples(const std::array<float, PRESAMPLECOUNT>& preamble, const std::vector<ModulatedFrame>& modulated_frames, std::vector<float>& transmission_samples)
{
    transmission_samples.clear();
    transmission_samples.reserve(TOTALSAMPLES);

    for (const auto& frame : modulated_frames)
    {
        transmission_samples.insert(
            transmission_samples.end(),
            preamble.begin(),
            preamble.end()
        );
        transmission_samples.insert(
            transmission_samples.end(),
            frame.samples.begin(),
            frame.samples.end()
        );
        transmission_samples.insert(
            transmission_samples.end(),
            GAPSAMPLECOUNT,
            0.0f
        );
    }
}

class TransmissionPlayer final : public juce::AudioIODeviceCallback
{
public:
    explicit TransmissionPlayer(const std::vector<float>& samples) : transmission_samples(samples)
    {
    }

    void audioDeviceAboutToStart(juce::AudioIODevice* device) override
    {
        playback_pos = 0;
        device_stopped = false;
        playback_finished = false;
    }

    void audioDeviceStopped() override
    {
        device_stopped = true;
        playback_finished = true;
    }

    void audioDeviceIOCallbackWithContext(
        const float* const* inputChannelData,
        int numInputChannels,
        float* const* outputChannelData,
        int numOutputChannels,
        int numSamples,
        const juce::AudioIODeviceCallbackContext& context) override
    {
        if (playback_finished)
        {
            return;
        }

        const std::size_t remaining_samples = transmission_samples.size() - playback_pos;
        const std::size_t samples_to_copy = std::min(remaining_samples, static_cast<std::size_t>(numSamples));
        for (int channel = 0; channel < numOutputChannels; channel++)
        {
            std::copy_n(
                transmission_samples.data() + playback_pos,
                samples_to_copy,
                outputChannelData[channel]
            );
        }
        playback_pos += samples_to_copy;

        if (playback_pos == transmission_samples.size())
        {
            playback_finished = true;
        }
    }

    bool isPlaybackFinished() const
    {
        return playback_finished.load();
    }

    bool wasDeviceStopped() const
    {
        return device_stopped.load();
    }
private:
    const std::vector<float>& transmission_samples;
    std::size_t playback_pos = 0;

    std::atomic<bool> playback_finished{ false };
    std::atomic<bool> device_stopped{ false };
};

int main (int argc, char* argv[])
{
    const std::filesystem::path input_path = argc >= 2 ? std::filesystem::path(argv[1]) : std::filesystem::path("INPUT.txt");
    std::vector<std::uint8_t> bits;

    if (!loadBitsFromFile(input_path, bits))
    {
        return 1;
    }

    std::cout << "Successfully loaded " << bits.size() << " bits.\n";

    std::vector<PayloadFrame> frames;
    
    splitIntoFrames(bits, frames);

    calculateFrameCrcs(frames);

    const auto encoded_frames = buildEncodedFrames(frames);

    const auto preamble = generatePreamble();

    const auto carrier = generateCarrierTemplate();

    const auto modulated_frames = modulateAllFrames(encoded_frames, carrier);

    std::vector<float> transmission_samples;

    buildTransmissionSamples(preamble, modulated_frames, transmission_samples);

    juce::AudioDeviceManager device_manager;
    device_manager.initialise(0, 1, nullptr, true);
    auto setup = device_manager.getAudioDeviceSetup();
    setup.sampleRate = SAMPLERATE;
    device_manager.setAudioDeviceSetup(setup, true);
    auto* current_device = device_manager.getCurrentAudioDevice();

    std::cout << "Press ENTER to begin transmission.\n";
    std::cin.get();

    TransmissionPlayer player(transmission_samples);

    const auto start_time = std::chrono::steady_clock::now();
    device_manager.addAudioCallback(&player);

    while (!player.isPlaybackFinished())
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    const auto finish_time = std::chrono::steady_clock::now();
    device_manager.removeAudioCallback(&player);
    const std::chrono::duration<double> elapsed = finish_time - start_time;

    std::cout << "Transmission completed.\n";
    std::cout << "Actual playback time: " << elapsed.count() << " seconds\n";

    return 0;
}
