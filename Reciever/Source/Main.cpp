#include <JuceHeader.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace {
// ---------------------------------------------------------------------------
// Protocol constants (must stay in sync with the transmitter).
// ---------------------------------------------------------------------------
constexpr double TWO_PI = 6.283185307179586476925286766559;
constexpr double SAMPLE_RATE = 48000.0;

constexpr std::size_t BITS_PER_FRAME = 100;
constexpr std::size_t CRC_BIT_COUNT = 8;
constexpr std::size_t ENCODED_BITS_PER_FRAME = BITS_PER_FRAME + CRC_BIT_COUNT; // 108

constexpr double CARRIER_FREQ = 10000.0;
constexpr std::size_t SAMPLES_PER_BIT = 48; // 48000 / 48 = 1000 baud

constexpr std::size_t PRE_SAMPLE_COUNT = 480;
constexpr std::size_t PRE_HALF_SAMPLE_COUNT = PRE_SAMPLE_COUNT / 2;
constexpr double PRE_START_FREQ = 2000.0;
constexpr double PRE_PEAK_FREQ = 10000.0;
constexpr float PRE_AMPLI = 0.5f;

const std::array<std::uint8_t, CRC_BIT_COUNT + 1> CRC_POLY{1, 1, 0, 1, 0, 0, 1, 1, 1};

// ---------------------------------------------------------------------------
// Preamble template: identical chirp to the transmitter
// (2 kHz -> 10 kHz on the first half, 10 kHz -> 2 kHz on the second half).
// ---------------------------------------------------------------------------
double preambleFrequency(std::size_t sampleIndex) {
    if (sampleIndex < PRE_HALF_SAMPLE_COUNT) {
        const double ratio = static_cast<double>(sampleIndex) / static_cast<double>(PRE_HALF_SAMPLE_COUNT - 1);
        return PRE_START_FREQ + ratio * (PRE_PEAK_FREQ - PRE_START_FREQ);
    }

    sampleIndex -= PRE_HALF_SAMPLE_COUNT;
    const double ratio = static_cast<double>(sampleIndex) / static_cast<double>(PRE_HALF_SAMPLE_COUNT - 1);
    return PRE_PEAK_FREQ + ratio * (PRE_START_FREQ - PRE_PEAK_FREQ);
}

std::array<float, PRE_SAMPLE_COUNT> buildPreambleTemplate() {
    std::array<float, PRE_SAMPLE_COUNT> preamble{};
    preamble[0] = 0.0f;

    double phase = 0.0;
    for (std::size_t i = 1; i < PRE_SAMPLE_COUNT; ++i) {
        const double previous = preambleFrequency(i - 1);
        const double current = preambleFrequency(i);
        phase += TWO_PI * 0.5 * (previous + current) / SAMPLE_RATE;
        if (phase >= TWO_PI) phase -= TWO_PI;

        preamble[i] = PRE_AMPLI * static_cast<float>(std::sin(phase));
    }
    return preamble;
}

std::array<double, SAMPLES_PER_BIT> buildCarrierCos() {
    std::array<double, SAMPLES_PER_BIT> c{};
    for (std::size_t i = 0; i < SAMPLES_PER_BIT; ++i)
        c[i] = std::cos(TWO_PI * CARRIER_FREQ * static_cast<double>(i) / SAMPLE_RATE);
    return c;
}

std::array<double, SAMPLES_PER_BIT> buildCarrierSin() {
    std::array<double, SAMPLES_PER_BIT> s{};
    for (std::size_t i = 0; i < SAMPLES_PER_BIT; ++i)
        s[i] = std::sin(TWO_PI * CARRIER_FREQ * static_cast<double>(i) / SAMPLE_RATE);
    return s;
}

// ---------------------------------------------------------------------------
// CRC check: identical algorithm to the transmitter (poly 110100111, MSB first).
// Returns true when the 8 CRC bits match the polynomial remainder.
// ---------------------------------------------------------------------------
bool crcMatches(const std::array<std::uint8_t, ENCODED_BITS_PER_FRAME>& bits) {
    std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> working{};
    std::copy_n(bits.begin(), BITS_PER_FRAME, working.begin());

    for (std::size_t start = 0; start < BITS_PER_FRAME; ++start) {
        if (working[start] == 0) continue;

        for (std::size_t k = 0; k < CRC_POLY.size(); ++k) working[start + k] ^= CRC_POLY[k];
    }

    for (std::size_t k = 0; k < CRC_BIT_COUNT; ++k)
        if (working[BITS_PER_FRAME + k] != bits[BITS_PER_FRAME + k]) return false;

    return true;
}

// ---------------------------------------------------------------------------
// BPSK demodulation with blind carrier-phase estimation.
// ---------------------------------------------------------------------------
struct DemodResult {
    std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> bits{};
    double confidence = 0.0;
    bool crcOk = false;
    // Set when the CRC only passed after the 180-degree ambiguity was resolved
    // by inverting the frame. Used later to infer the polarity of the frames
    // whose own CRC failed on both branches.
    bool flipped = false;
};

DemodResult demodulateFrame(const std::vector<float>& rx,
                            std::size_t dataStart,
                            const std::array<double, SAMPLES_PER_BIT>& cosRef,
                            const std::array<double, SAMPLES_PER_BIT>& sinRef) {
    std::array<double, ENCODED_BITS_PER_FRAME> re{};
    std::array<double, ENCODED_BITS_PER_FRAME> im{};

    for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
        double a = 0.0;
        double b = 0.0;
        for (std::size_t n = 0; n < SAMPLES_PER_BIT; ++n) {
            const double sample = rx[dataStart + j * SAMPLES_PER_BIT + n];
            a += sample * cosRef[n];
            b += sample * sinRef[n];
        }
        re[j] = a;
        im[j] = -b;
    }

    // Build two candidate carrier-phase models:
    //  (1) a single global phase (robust when the two clocks agree);
    //  (2) a per-block linear ramp (tracks a residual frequency offset).
    // We demodulate with both and keep the one that passes the CRC, otherwise
    // the one with the higher soft-decision confidence. The global model is
    // always a candidate, so this can never be worse than using it alone.
    constexpr std::size_t BLOCK = 12;
    constexpr std::size_t NUM_BLOCKS = (ENCODED_BITS_PER_FRAME + BLOCK - 1) / BLOCK;
    constexpr double PI = 3.14159265358979323846;

    // (1) global doubled phase.
    double sumRe2 = 0.0;
    double sumIm2 = 0.0;
    for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
        sumRe2 += re[j] * re[j] - im[j] * im[j];
        sumIm2 += 2.0 * re[j] * im[j];
    }
    const double globalTwoPhi = std::atan2(sumIm2, sumRe2);

    // (2) block-wise doubled phase, unwrapped and linearly interpolated.
    std::array<double, NUM_BLOCKS> twoPhi{};
    for (std::size_t b = 0; b < NUM_BLOCKS; ++b) {
        double blockRe2 = 0.0;
        double blockIm2 = 0.0;
        const std::size_t begin = b * BLOCK;
        const std::size_t end = std::min(begin + BLOCK, ENCODED_BITS_PER_FRAME);
        for (std::size_t j = begin; j < end; ++j) {
            blockRe2 += re[j] * re[j] - im[j] * im[j];
            blockIm2 += 2.0 * re[j] * im[j];
        }
        twoPhi[b] = std::atan2(blockIm2, blockRe2);
    }

    for (std::size_t b = 1; b < NUM_BLOCKS; ++b) {
        while (twoPhi[b] - twoPhi[b - 1] > PI) twoPhi[b] -= TWO_PI;
        while (twoPhi[b] - twoPhi[b - 1] < -PI) twoPhi[b] += TWO_PI;
    }

    const auto demodWith = [&](bool useRamp) {
        DemodResult candidate;
        for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
            double twoPhiAtBit = globalTwoPhi;
            if (useRamp) {
                const double pos = (static_cast<double>(j) + 0.5) / static_cast<double>(BLOCK) - 0.5;
                if (pos <= 0.0) {
                    twoPhiAtBit = twoPhi[0];
                } else if (pos >= static_cast<double>(NUM_BLOCKS - 1)) {
                    twoPhiAtBit = twoPhi[NUM_BLOCKS - 1];
                } else {
                    const auto bLow = static_cast<std::size_t>(pos);
                    const double frac = pos - static_cast<double>(bLow);
                    twoPhiAtBit = (1.0 - frac) * twoPhi[bLow] + frac * twoPhi[bLow + 1];
                }
            }

            const double phi = 0.5 * twoPhiAtBit;
            const double projection = re[j] * std::cos(phi) + im[j] * std::sin(phi);
            candidate.bits[j] = projection > 0.0 ? 1 : 0;
            candidate.confidence += std::abs(projection);
        }

        candidate.crcOk = crcMatches(candidate.bits);

        // Resolve a possible global 180 degree ambiguity.
        if (!candidate.crcOk) {
            std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> flipped{};
            for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) flipped[j] = candidate.bits[j] ? 0 : 1;

            if (crcMatches(flipped)) {
                candidate.bits = flipped;
                candidate.crcOk = true;
                candidate.flipped = true;
            }
        }
        return candidate;
    };

    const DemodResult globalResult = demodWith(false);
    const DemodResult rampResult = demodWith(true);

    // Prefer the global-phase result (stable baseline). The ramp model is only
    // adopted when it actually makes the CRC pass, so it can only help.
    if (globalResult.crcOk) return globalResult;
    if (rampResult.crcOk) return rampResult;
    return globalResult;
}

// ---------------------------------------------------------------------------
// Complex baseband helpers.
//
// The matched filter below can only pin the carrier phase down modulo 180
// degrees (BPSK carries an inherent sign ambiguity). Once a frame is known to
// be correct - because its CRC passed - its bits are known too, so the sign
// ambiguity disappears and the exact carrier phase can be recovered. That
// phase is reused to demodulate the frames whose CRC failed.
// ---------------------------------------------------------------------------
struct Baseband {
    std::array<double, ENCODED_BITS_PER_FRAME> re{};
    std::array<double, ENCODED_BITS_PER_FRAME> im{};
};

Baseband computeBaseband(const std::vector<float>& rx,
                         std::size_t dataStart,
                         const std::array<double, SAMPLES_PER_BIT>& cosRef,
                         const std::array<double, SAMPLES_PER_BIT>& sinRef) {
    Baseband bb;
    for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
        double a = 0.0;
        double b = 0.0;
        for (std::size_t n = 0; n < SAMPLES_PER_BIT; ++n) {
            const double sample = rx[dataStart + j * SAMPLES_PER_BIT + n];
            a += sample * cosRef[n];
            b += sample * sinRef[n];
        }
        bb.re[j] = a;
        bb.im[j] = -b;
    }
    return bb;
}

void projectWithPhase(const Baseband& bb,
                      double phi,
                      std::array<std::uint8_t, ENCODED_BITS_PER_FRAME>& bits,
                      double& confidence) {
    const double c = std::cos(phi);
    const double s = std::sin(phi);
    confidence = 0.0;
    for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
        const double projection = bb.re[j] * c + bb.im[j] * s;
        bits[j] = projection > 0.0 ? 1 : 0;
        confidence += std::abs(projection);
    }
}

double phaseFromKnownBits(const Baseband& bb, const std::array<std::uint8_t, ENCODED_BITS_PER_FRAME>& bits) {
    double sr = 0.0;
    double si = 0.0;
    for (std::size_t j = 0; j < ENCODED_BITS_PER_FRAME; ++j) {
        const double sign = bits[j] ? 1.0 : -1.0;
        sr += sign * bb.re[j];
        si += sign * bb.im[j];
    }
    return std::atan2(si, sr);
}

// Re-demodulate a frame using a phase borrowed from a nearby CRC-valid frame,
// fine-searching the sample timing and refining the phase in a decision-directed
// loop. A frame that failed its own CRC was not necessarily sampled at the
// optimal position, so the timing is re-searched around the borrowed estimate.
struct FineResult {
    std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> bits{};
    double confidence = 0.0;
    bool ok = false;
};

FineResult demodulateWithBorrowedPhase(const std::vector<float>& rx,
                                       std::size_t center,
                                       int searchRadius,
                                       double seedPhase,
                                       std::size_t frameDataSamples,
                                       const std::array<double, SAMPLES_PER_BIT>& cosRef,
                                       const std::array<double, SAMPLES_PER_BIT>& sinRef) {
    FineResult best;
    for (int offset = -searchRadius; offset <= searchRadius; ++offset) {
        const std::ptrdiff_t signedStart = static_cast<std::ptrdiff_t>(center) + offset;
        if (signedStart < 0) continue;

        const auto start = static_cast<std::size_t>(signedStart);
        if (start + frameDataSamples > rx.size()) continue;

        const Baseband bb = computeBaseband(rx, start, cosRef, sinRef);

        std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> bits{};
        double confidence = 0.0;
        double phase = seedPhase;
        for (int iteration = 0; iteration < 4; ++iteration) {
            projectWithPhase(bb, phase, bits, confidence);
            phase = phaseFromKnownBits(bb, bits);
        }

        if (!best.ok || confidence > best.confidence) {
            best.bits = bits;
            best.confidence = confidence;
            best.ok = true;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// Frame synchronisation: normalised sliding cross-correlation between the
// received signal and the preamble template. Peaks mark frame starts.
// ---------------------------------------------------------------------------
std::vector<std::size_t> findPreamblePeaks(const std::vector<float>& rx,
                                           const std::array<float, PRE_SAMPLE_COUNT>& preamble) {
    std::vector<std::size_t> peaks;
    const std::size_t n = PRE_SAMPLE_COUNT;
    if (rx.size() < n) return peaks;

    const std::size_t m = rx.size() - n + 1;

    double preambleEnergy = 0.0;
    for (auto v : preamble) preambleEnergy += static_cast<double>(v) * v;

    std::vector<double> normCorr(m, 0.0);

    double windowEnergy = 0.0;
    for (std::size_t i = 0; i < n; ++i) windowEnergy += static_cast<double>(rx[i]) * rx[i];

    double maxNorm = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        if (i > 0)
            windowEnergy +=
                static_cast<double>(rx[i + n - 1]) * rx[i + n - 1] - static_cast<double>(rx[i - 1]) * rx[i - 1];

        double corr = 0.0;
        for (std::size_t k = 0; k < n; ++k) corr += static_cast<double>(rx[i + k]) * preamble[k];

        const double denom = std::sqrt(std::max(windowEnergy, 1e-9) * preambleEnergy);
        const double value = corr / denom;
        normCorr[i] = value;
        maxNorm = std::max(maxNorm, value);
    }

    const double threshold = std::max(0.25, 0.5 * maxNorm);

    struct Candidate {
        double value;
        std::size_t index;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < m; ++i)
        if (normCorr[i] > threshold) candidates.push_back({normCorr[i], i});

    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.value > b.value; });

    const std::size_t minSeparation = 2000;
    for (const auto& candidate : candidates) {
        bool tooClose = false;
        for (auto idx : peaks) {
            const auto diff = candidate.index > idx ? candidate.index - idx : idx - candidate.index;
            if (diff < minSeparation) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) peaks.push_back(candidate.index);
    }

    std::sort(peaks.begin(), peaks.end());
    return peaks;
}

std::vector<std::uint8_t> decodeTrack(const std::vector<float>& rx, int& detectedFrames, int& crcPassFrames) {
    detectedFrames = 0;
    crcPassFrames = 0;

    const auto preamble = buildPreambleTemplate();
    const auto cosRef = buildCarrierCos();
    const auto sinRef = buildCarrierSin();

    const auto peaks = findPreamblePeaks(rx, preamble);
    std::cout << "Detected " << peaks.size() << " preamble(s).\n";

    const std::size_t frameDataSamples = ENCODED_BITS_PER_FRAME * SAMPLES_PER_BIT;

    struct FrameInfo {
        std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> bits{};
        std::size_t start = 0;
        bool crcOk = false;
        bool hasPhase = false;
        double phase = 0.0;
    };

    std::vector<FrameInfo> frameInfos;
    frameInfos.reserve(peaks.size());

    for (auto peak : peaks) {
        const std::size_t base = peak + PRE_SAMPLE_COUNT;
        if (base + frameDataSamples > rx.size()) continue;

        DemodResult best;
        double bestConfidence = -1.0;
        std::size_t bestStart = base;

        // Try small offsets around the correlation peak to absorb timing error.
        for (int offset = -4; offset <= 4; ++offset) {
            const std::ptrdiff_t signedStart = static_cast<std::ptrdiff_t>(base) + offset;
            if (signedStart < 0) continue;

            const auto start = static_cast<std::size_t>(signedStart);
            if (start + frameDataSamples > rx.size()) continue;

            auto result = demodulateFrame(rx, start, cosRef, sinRef);
            if (result.crcOk) {
                best = result;
                bestStart = start;
                break;
            }
            if (result.confidence > bestConfidence) {
                bestConfidence = result.confidence;
                best = result;
                bestStart = start;
            }
        }

        FrameInfo info;
        info.bits = best.bits;
        info.start = bestStart;
        info.crcOk = best.crcOk;

        if (best.crcOk) {
            // Recover the unambiguous carrier phase of this (known-good) frame.
            const Baseband bb = computeBaseband(rx, bestStart, cosRef, sinRef);
            info.phase = phaseFromKnownBits(bb, best.bits);
            info.hasPhase = true;
        }

        ++detectedFrames;
        if (best.crcOk) ++crcPassFrames;
        frameInfos.push_back(info);
    }

    std::vector<std::uint8_t> payloadBits;

    for (std::size_t i = 0; i < frameInfos.size(); ++i) {
        std::array<std::uint8_t, ENCODED_BITS_PER_FRAME> bits = frameInfos[i].bits;

        // A frame whose CRC failed on both branches has an unresolved 180
        // degree ambiguity, so the raw decision may be entirely inverted. Each
        // frame's data carrier restarts at phase zero, so the true carrier
        // phase is essentially identical from frame to frame; borrow it from
        // the nearest CRC-valid frame and re-demodulate this one.
        if (!frameInfos[i].crcOk) {
            std::size_t nearest = frameInfos.size();
            std::size_t nearestDist = 0;
            for (std::size_t k = 0; k < frameInfos.size(); ++k) {
                if (!frameInfos[k].hasPhase) continue;
                const std::size_t dist = k > i ? k - i : i - k;
                if (nearest == frameInfos.size() || dist < nearestDist) {
                    nearest = k;
                    nearestDist = dist;
                }
            }

            if (nearest != frameInfos.size()) {
                const FineResult refined = demodulateWithBorrowedPhase(
                    rx, frameInfos[i].start, 24, frameInfos[nearest].phase, frameDataSamples, cosRef, sinRef);
                if (refined.ok) bits = refined.bits;
            }
        }

        // INPUT.txt holds exactly 10000 bits = 100 frames; ignore extras.
        if (i < 100) std::copy_n(bits.begin(), BITS_PER_FRAME, std::back_inserter(payloadBits));
    }

    return payloadBits;
}

// ---------------------------------------------------------------------------
// Live microphone capture.
// ---------------------------------------------------------------------------
class AudioRecorder final : public juce::AudioIODeviceCallback {
public:
    void audioDeviceAboutToStart(juce::AudioIODevice *) override {}
    void audioDeviceStopped() override {}

    void audioDeviceIOCallbackWithContext(const float *const *inputChannelData,
                                          int numInputChannels,
                                          float *const *outputChannelData,
                                          int numOutputChannels,
                                          int numSamples,
                                          const juce::AudioIODeviceCallbackContext&) override {
        for (int channel = 0; channel < numOutputChannels; ++channel)
            if (outputChannelData[channel] != nullptr) std::fill_n(outputChannelData[channel], numSamples, 0.0f);

        if (numInputChannels > 0 && inputChannelData[0] != nullptr) {
            const std::lock_guard<std::mutex> lock(bufferMutex);
            buffer.insert(buffer.end(), inputChannelData[0], inputChannelData[0] + numSamples);
        }
    }

    std::vector<float> getBuffer() {
        const std::lock_guard<std::mutex> lock(bufferMutex);
        return buffer;
    }

private:
    std::vector<float> buffer;
    std::mutex bufferMutex;
};

std::vector<float> recordFromMicrophone() {
    juce::AudioDeviceManager deviceManager;
    deviceManager.initialise(1, 0, nullptr, true);

    auto setup = deviceManager.getAudioDeviceSetup();
    setup.sampleRate = SAMPLE_RATE;
    deviceManager.setAudioDeviceSetup(setup, true);

    AudioRecorder recorder;

    std::cout << "Press ENTER to start recording.\n";
    std::cin.get();

    deviceManager.addAudioCallback(&recorder);
    std::cout << "Recording... Press ENTER to stop.\n";
    std::cin.get();

    deviceManager.removeAudioCallback(&recorder);
    return recorder.getBuffer();
}

// ---------------------------------------------------------------------------
// Offline WAV loading (for reproducible debugging).
// ---------------------------------------------------------------------------
std::vector<float> loadWaveFile(const std::filesystem::path& path) {
    juce::AudioFormatManager formatManager;
    formatManager.registerBasicFormats();

    juce::File file(path.string());
    std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(file));
    if (reader == nullptr) {
        std::cerr << "Cannot read wave file: " << path << '\n';
        return {};
    }

    std::cout << "Loaded wave file, sample rate = " << reader->sampleRate << '\n';
    if (std::abs(reader->sampleRate - SAMPLE_RATE) > 1.0)
        std::cerr << "Warning: sample rate is not 48000 Hz; decoding may fail.\n";

    const auto length = static_cast<int>(reader->lengthInSamples);
    juce::AudioBuffer<float> buffer(1, length);
    reader->read(&buffer, 0, length, 0, true, false);

    std::vector<float> samples(static_cast<std::size_t>(length));
    std::copy_n(buffer.getReadPointer(0), length, samples.begin());
    return samples;
}

void writeOutputFile(const std::filesystem::path& path, const std::vector<std::uint8_t>& bits) {
    std::ofstream output(path, std::ios::binary);
    for (auto bit : bits) output << (bit ? '1' : '0');

    std::cout << "Wrote " << bits.size() << " bits to " << path << '\n';
}
} // namespace

//==============================================================================
int main(int argc, char *argv[]) {
    std::vector<float> rx;

    // Usage:
    //   Reciever.exe                 -> record from microphone
    //   Reciever.exe TX.wav          -> decode a wave file, write OUTPUT.txt
    //   Reciever.exe TX.wav OUT.txt  -> decode a wave file, write OUT.txt
    if (argc >= 2)
        rx = loadWaveFile(std::filesystem::path(argv[1]));
    else
        rx = recordFromMicrophone();

    if (rx.empty()) {
        std::cerr << "No audio samples available.\n";
        return 1;
    }

    int detectedFrames = 0;
    int crcPassFrames = 0;
    const auto payloadBits = decodeTrack(rx, detectedFrames, crcPassFrames);

    std::cout << "Frames detected: " << detectedFrames << ", CRC passed: " << crcPassFrames << '\n';

    const std::filesystem::path outputPath =
        argc >= 3 ? std::filesystem::path(argv[2]) : std::filesystem::path("OUTPUT.txt");

    writeOutputFile(outputPath, payloadBits);

    return 0;
}
