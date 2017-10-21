#pragma once

#include <vector>
#include <xmmintrin.h>

#include "helpers/Common.h"
#include "windows/timer.h"

#include "windows/WaveOut_Simple.h"
#include "sfx/Oscillators.h"
#include "sfx/Filters.h"
#include "sfx/Reverb.h"
#include "sfx/ADSR.h"
#include "sfx/Delay.h"
#include "sfx/FiltersDelay.h"
#include "sfx/PlotFilter.h"

#include "Base.h"

#define PARAM_DECLARE(param, val, min, max) \
	scalar param = val; \
	const scalar param##Min = min; \
	const scalar param##Max = max;

namespace KarplusStrongSettings
{
	bool MouseOver = false;
	bool ClassicFilter = false;
	bool PlayMusic = false;

	PARAM_DECLARE(Damping, 0.0f, 0.0f, 1.0f);
	PARAM_DECLARE(Richness, 0.45f, 0.0f, 1.0f);
	PARAM_DECLARE(NoiseMix, 0.5f, 0.0f, 1.0f);
	PARAM_DECLARE(PluckStr, 0.1f, 0.0f, 1.0f);

	PARAM_DECLARE(TuneGain, 0.0f, -1.0f, 1.0f);
	PARAM_DECLARE(TuneDelay, 0.0f, 0.0f, 10.0f);

	PARAM_DECLARE(PickPos, 0.5f, 0.0f, 1.0f);
	PARAM_DECLARE(PickUpPos, 0.6f, 0.0f, 1.0f);
	PARAM_DECLARE(PickStrength, 0.0f, 0.0f, 1.0f);

	PARAM_DECLARE(VibrFreq, 30.0f, 0.0f, 100.0f);
	PARAM_DECLARE(VibrAmp, 0.0f, 0.0f, 0.25f);

	PARAM_DECLARE(S06_Freq, 82.4f, 20.0f, 1000.0f);
	PARAM_DECLARE(S06_LowPassFreq, 1000.0f, 10.0f, 5000.0f);

	PARAM_DECLARE(S05_Freq, 110.0f, 20.0f, 1000.0f);
	PARAM_DECLARE(S05_LowPassFreq, 1000.0f, 10.0f, 5000.0f);

	PARAM_DECLARE(S04_Freq, 146.8f, 20.0f, 1000.0f);
	PARAM_DECLARE(S04_LowPassFreq, 1000.0f, 10.0f, 5000.0f);

	PARAM_DECLARE(S03_Freq, 196.0f, 20.0f, 1000.0f);
	PARAM_DECLARE(S03_LowPassFreq, 1000.0f, 10.0f, 5000.0f);

	PARAM_DECLARE(S02_Freq, 246.9f, 20.0f, 1000.0f);
	PARAM_DECLARE(S02_LowPassFreq, 1000.0f, 10.0f, 5000.0f);

	PARAM_DECLARE(S01_Freq, 329.6f, 20.0f, 1000.0f);
	PARAM_DECLARE(S01_LowPassFreq, 1000.0f, 10.0f, 5000.0f);
}

inline scalar calcFretFreq(scalar baseSemitone, int idxFret = 1)
{
	return powf(2.0f, (baseSemitone + idxFret - 49) / 12.0f) * 440.0f;
}

inline scalar calcFretFreq_BaseFreq(scalar baseFreq, int idxFret = 1)
{
	// First, calculate the base semitone (when the string is open):
	//		baseN = 12 * log2( F(n) / 440.0 ) + 49
	scalar baseSemitone = 12 * m_log2(baseFreq / 440.0f) + 49;
	// One fret moves 1 note down (12 frets: 12 semitone down = full octave)
	// hence we can next calculate the fret freq by:
	//		F(fret) = 2^((baseN + fretNum - 49) / 12) * 440.0
	return calcFretFreq(baseSemitone, idxFret);
}

inline scalar calcFreqStringFret(int idxString = 1, int idxFret = 1)
{
	const scalar stBaseFreqs[6] =
	{
		329.6f,
		246.9f,
		196.0f,
		146.8f,
		110.0f,
		82.4f
	};
	return calcFretFreq_BaseFreq(stBaseFreqs[idxString-1], idxFret);
}

class SoundKarplusStrong : public SoundPlay
{
protected:

	static const int c_numDevices = 6;			// 6 strings
	static const int c_samplesPerSec = 44100;
	static const int c_numChannels = 2;

	uint m_numMarginSamples;
	size_t m_bufSize;
	float * m_buffer;
	short * m_bufInt;

	sfx::LowPassIIR m_lowPass;
	sfx::LowPassIIR m_dynamicsLowPass;
	sfx::LowPassIIR m_finalLowPass;

	sfx::DelayFrac m_KSDelay;
	sfx::CombDelayFF m_pickPosComb;
	sfx::CombDelayFF m_pickupPosComb;

	// Delays for constructing the feed-back/feed-forward all-pass filter
	sfx::Delay m_allPassDelay_Tune;

	bool m_isGenerated = false;

public:

	int genType;

	SoundKarplusStrong(scalar numSeconds = 3.0f):
		SoundPlay(c_numDevices, c_samplesPerSec, c_numChannels),
		genType(1)
	{
		m_numSeconds = numSeconds;

		m_KSDelay.setParams(4096, 0);
		m_allPassDelay_Tune.setParams(4096, 0);

		m_numMarginSamples = (uint)(1.0f * m_samplesPerSec * m_numChannels); 
		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec * m_numChannels);

		m_bufSize = toneBlockSize + m_numMarginSamples;
		m_buffer = new float [m_bufSize];
		m_bufInt = new short [m_bufSize];
	}

	~SoundKarplusStrong()
	{
		delete [] m_buffer;
		delete [] m_bufInt;
	}

	/*
	Table of string freqs by fret:
	
	Fret	S6		S5		S4		S3		S2		S1

	 0		82.4	110.0	146.8	196.0	246.9	329.6
	 1		87.3	116.5	155.6	207.6	261.6	349.2
	 2		92.5	123.5	164.8	220.0	277.2	370.0
	 3		98.0	130.8	174.6	233.1	293.6	392.0
	 4		103.8	138.6	185.0	246.9	311.1	415.3
	//*/

	void KS_GenPluck(
		scalar freq,
		uint sampleLength,
		float noiseMix
		)
	{
		float initialSample;
		for (uint i = 0; i < sampleLength; ++i)
		{
			// Noise sample
			uint rnd30 = rand30() & 65535;
			float noise = (rnd30 - 32768.0f) / 32768.0f;

			float triangleAdd = sfx::oscTriangle(freq * i / (float)m_samplesPerSec, 0.99f);
			noise += (triangleAdd - noise) * (1.0f - noiseMix);

			noise = m_pickPosComb.update(noise);

			m_dynamicsLowPass.applyCont(&initialSample, noise);
			//bandPassTemp.ApplyCont(&initialSample, noise);

#if 0
			mKSDelay.SetSample(initialSample);
			mKSDelay.GetSample();
#else
			m_buffer[i*2 + 0] = initialSample;
			m_buffer[i*2 + 1] = initialSample;
#endif
		}
	}

	void KS_Resonate(
		float baseDelay,
		scalar dampingGain,
		scalar filterRichness,
		scalar vibratoFreq,
		scalar vibratoAmp,
		bool lowOrderFilter,
		uint blockSize,
		scalar tuneGain,
		scalar tuneDelay
		)
	{
		uint outputIdx = 0;
		float prevSample = 0.0f;
		while (outputIdx < blockSize)
		{
			float sample, out;

			// Decay coeff. is [0.7..0.865] (lower is more tension)
			float stringT = 0.9f;
			float stringBase = 0.9f;	// 0.7
			//		sample = (stringBase + stringT * 0.16545f) * ksDelay.GetSample();

			scalar phase = outputIdx / (scalar)m_samplesPerSec;
			scalar vibratoCoeff = vibratoAmp * sinf(vibratoFreq * phase);
#if 0
			// Integer delay line
			mKSDelay.Resize((float)sampleLength + vibratoCoeff);
#else
			// Fractional delay line
			m_KSDelay.resize(baseDelay + vibratoCoeff);
#endif

			sample = m_KSDelay.getSample();

			if (!lowOrderFilter)
			{
				// Higher order low-pass
				sample = dampingGain * sample;
				m_lowPass.applyCont(&out, sample);
				out += m_buffer[outputIdx*2 + 0];
				//out *= dampingGain;
			}
			else
			{
				// KS classic first order low-pass
				sample = dampingGain * sample;
				const scalar lowPassAlpha = filterRichness;//0.55f;
				out = m_buffer[outputIdx*2 + 0] + prevSample * (1.0f - lowPassAlpha) + sample * lowPassAlpha;
				prevSample = out;
			}

			// Tune allpass filter
			//////////////////////////////////////////////////////////////////////////
			if (tuneDelay > 0.01f)
			{
				float tuneDelaySample = m_allPassDelay_Tune.getSample();
				float tuneAllpassFeedback = tuneDelaySample * tuneGain;

				float tuneDelayInput = out + tuneAllpassFeedback;
				scalar tuneAllpassFeedforward = tuneDelayInput * -tuneGain;

				float allPassOutput = tuneDelaySample + tuneAllpassFeedforward;

				m_allPassDelay_Tune.setSample(tuneDelayInput);

				out = allPassOutput;
			}
			//////////////////////////////////////////////////////////////////////////

			m_KSDelay.setSample(out);

#if 0
			// Removed as it is almost uneffective

			// Inharmonicity allpass filter
			//////////////////////////////////////////////////////////////////////////
			if (inharmDelay > 0.01f)
			{
				float inharmDelaySample = mAllPassDelay_Inharm.GetSample();
				float inharmAllpassFeedback = inharmDelaySample * inharmGain;

				float inharmDelayInput = out + inharmAllpassFeedback;
				scalar inharmAllpassFeedforward = inharmDelayInput * -inharmGain;

				float allPassOutput = inharmDelaySample + inharmAllpassFeedforward;

				mAllPassDelay_Inharm.SetSample(inharmDelayInput);

				out = allPassOutput;
			}
			//////////////////////////////////////////////////////////////////////////
#endif

			out = m_pickupPosComb.update(out);

			m_buffer[outputIdx*2 + 0] = out;
			m_buffer[outputIdx*2 + 1] = out;

			++outputIdx;
		}
	}

	void KS_MakeBuffer(
		uint sampleLength,
		uint blockSize,
		uint startMargin,
		scalar maxVolume
		)
	{
		//mFinalLowPass.setParams(1000.0f / m_samplesPerSec, 0.4f);
		//mFinalLowPass.AutoNormalize(0.0f / m_samplesPerSec);
		//mFinalLowPass.ResetHistory();

		for (uint j = 0; j < startMargin * m_numChannels; ++j)
		{
			m_bufInt[j] = 0;
		}

		const uint intBufShift = startMargin * m_numChannels;

		// Convert float buf into shorts
		// UPD: we need to skip some samples (initial noise)
		const scalar volume = 0.9f / maxVolume;
		const uint skipSamples = sampleLength * m_numChannels;
		for (uint j = 0; j < blockSize * m_numChannels; ++j)
		{
			float curBufferElement = 0.0f;
			if (j < blockSize * m_numChannels - skipSamples)
			{
				curBufferElement = m_buffer[j + skipSamples];
			}
			//mFinalLowPass.ApplyCont(&curBufferElement, curBufferElement);
			m_bufInt[j + intBufShift] = int(volume * clamp(curBufferElement, -1.0f, 1.0f ) * 32767);
		}
	}

	void KS_PlayBuffer(
		uint blockSize,
		uint startMargin,
		uint idxString
		)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		// Output converted buffer
		if (m_waveOut[idxString].getBufferNum() != 0)
		{
			m_waveOut[idxString].reset();
		}
		m_waveOut[idxString].play((char *)m_bufInt, (blockSize + startMargin) * m_numChannels * sizeof(short));

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
	}

	void generate(
		int idxString,
		scalar stringFreq,
		scalar lowPassFreq,
		scalar pluckStrength,
		bool lowOrderFilter,
		scalar damping,
		scalar filterRichness,
		scalar noiseMix,
		scalar vibrAmp,
		scalar vibrFreq,

		scalar tuneGain,
		scalar tuneDelay,

		scalar startMarginMS = 0
		)
	{
		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);
		const uint startMargin = (uint)(startMarginMS * m_samplesPerSec);

		const scalar dutyTime = 0.5f;

		const scalar freq = stringFreq;

		m_perfTimer.start();

		scalar maxVolume = 0.0f;

		uint sampleLength = (uint)(m_samplesPerSec / freq + 1);
		// Better sound quality with this one (!)
		sampleLength *= 2;

#if 0
		// Integer delay line
		mKSDelay.init(sampleLength);
#else
		// Calculate phase delay to compensate for tuning all-pass
		scalar tuneAllPassCompensation = 0;

		if (tuneDelay > 0.01f)
		{
			// Order of the all-pass filter
			scalar M = floorf(tuneDelay);

			scalar omegaBase = freq * _2PI / m_samplesPerSec;
			scalar tuneGain2 = tuneGain*tuneGain;
			scalar phase_delay = 1/omegaBase * atan2f( (1-tuneGain2) * sinf(M*omegaBase), (1+tuneGain2) * cosf(M*omegaBase) - 2*tuneGain );

			// Not needed
			scalar group_delay = (1-tuneGain2) / ( 1 - 2*tuneGain*cosf(omegaBase) + tuneGain2 );

			//tuneAllPassCompensation = sqrtf( floorf(tuneDelay) ) * phase_delay;
			tuneAllPassCompensation = phase_delay;
		}

		// Fractional delay line - No need to round off
		scalar baseDelay = m_samplesPerSec / freq + 1 - tuneAllPassCompensation;
		m_KSDelay.init(baseDelay);
#endif
		for (uint i = 0; i < blockSize * m_numChannels; ++i)
		{
			m_buffer[i] = 0.0f;
		}

		// Set up poick position feed-forward comb filter
		scalar pickPos_uN = KarplusStrongSettings::PickPos * sampleLength;
		m_pickPosComb.init(floorf(pickPos_uN));
#if 1
		// Feed-forward comb
		m_pickPosComb.setGain(1.0f - KarplusStrongSettings::PickStrength, 1.0f);
#else
		// Feed-back comb
		mPickPosComb.SetGain(KarplusStrongSettings::PickStrength);
#endif

		scalar pickUpPos_uN = KarplusStrongSettings::PickUpPos * sampleLength;
		m_pickupPosComb.init(floorf(pickUpPos_uN));
#if 1
		// Feed-forward comb
		m_pickupPosComb.setGain(1.0f - KarplusStrongSettings::PickStrength, 1.0f);
#else
		// Feed-back comb
		mPickupPosComb.SetGain(KarplusStrongSettings::PickStrength);
#endif

		// Set up delay for the inharmonicity all-pass filters, the number of samples describes the all-pass order
		//mAllPassDelay_Inharm.init((int)inharmDelay);
		m_allPassDelay_Tune.init((int)tuneDelay);

		// For slightly plucked string, it should be  1k / 0.8, since we need to filter almost any higher partial
		// For a hardly plucked string, it should be 16k / 0.6, to leave most higher partials
#if 0
#if 1
		// Slight
		mDynamicsLowPass.SetParams(1000.0f / m_samplesPerSec, 0.8f);
		mDynamicsLowPass.AutoNormalize(0.0f / m_samplesPerSec);
#else
		// Hard
		mDynamicsLowPass.SetParams(16000.0f / m_samplesPerSec, 0.6f);
		mDynamicsLowPass.AutoNormalize(0.0f / m_samplesPerSec);
#endif
#else
		// From 1k to 16k
		// TODO: Parameters of the filter should also comply to the condition that
		// fundamental frequency should have same amplitude across alll strings,
		// with static filter higher freq strings will have less amplitude
		const scalar dynamicsFilterFreq = 500.0f + 15500.0f * pluckStrength;
		const scalar dynamicsFilterQ = 0.6f + 0.2f * (1.0f - pluckStrength);
		m_dynamicsLowPass.setParams(dynamicsFilterFreq / m_samplesPerSec, dynamicsFilterQ);
		m_dynamicsLowPass.autoNormalize(0.0f / m_samplesPerSec);
#endif
		m_dynamicsLowPass.resetHistory();

		//PlotFreqResponse("ksDynamicsFilterPlot.tga", dynamicsLowPass, m_samplesPerSec * 0.5f, 1024, 256);

		KS_GenPluck(
			freq,
			sampleLength,
			noiseMix
			);

#if 0
		mLowPass.SetParams(2000.0f / m_samplesPerSec, 0.8f);
		mLowPass.AutoNormalize(1200.0f / m_samplesPerSec);
#else
		// From 2k to 16k
		const scalar loopFilterFreq = 2000.0f + 14000.0f * filterRichness;
		const scalar loopFilterQ = 0.2f + 0.4f * (1.0f - filterRichness);
		m_lowPass.setParams(loopFilterFreq / m_samplesPerSec, loopFilterQ);
		m_lowPass.autoNormalize(0.0f / m_samplesPerSec);
#endif
		m_lowPass.resetHistory();

		//PlotFreqResponse("ksFilterPlot.tga", lowPass, m_samplesPerSec * 0.5f, 1024, 256);

		// Nonlinear damping gain
		scalar dampingLerp = clamp<scalar>(freq / 390.0f, 0.0f, 0.9f);
		dampingLerp *= dampingLerp;
		const scalar dampingGain = (1.0f - damping) * (0.98f + 0.02f * dampingLerp);

		const scalar vibratoFreq = vibrFreq;
		const scalar vibratoAmp = sampleLength * vibrAmp;

		KS_Resonate(
			baseDelay,
			dampingGain,
			filterRichness,
			vibratoFreq,
			vibratoAmp,
			lowOrderFilter,
			blockSize,
			tuneGain,
			tuneDelay
			);

		maxVolume = 1.0f;

		KS_MakeBuffer(
			sampleLength,
			blockSize,
			startMargin,
			maxVolume
			);

		m_isGenerated = true;
		m_execTime = (scalar)m_perfTimer.time();
	}
	void generate() override {  }
	void reset()
	{
		m_isGenerated = false;
	}

	virtual short * getBufferGenerated()
	{
		if (m_isGenerated)
		{
			return m_bufInt;
		}
		else
		{
			return 0;
		}
	}
	virtual uint getBufferSizeGenerated()
	{
		if (m_isGenerated)
		{
			return (uint)(m_numSeconds * m_samplesPerSec * m_numChannels);
		}
		else
		{
			return 0;
		}
	}

	void playSound() override
	{
		const scalar startMarginMS = 0.0f;
		const uint idxString = 6;

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);
		const uint startMargin = (uint)(startMarginMS * m_samplesPerSec);

		generate(
			idxString-1,
			calcFreqStringFret(idxString, 3),
			1000.0f,
			KarplusStrongSettings::PluckStr,
			KarplusStrongSettings::ClassicFilter,
			KarplusStrongSettings::Damping,
			KarplusStrongSettings::Richness,
			KarplusStrongSettings::NoiseMix,
			KarplusStrongSettings::VibrAmp,
			KarplusStrongSettings::VibrFreq,
			KarplusStrongSettings::TuneGain,
			KarplusStrongSettings::TuneDelay,

			startMarginMS
			);

		KS_PlayBuffer(
			blockSize,
			startMargin,
			idxString-1
			);
	}
	virtual void storeSound(const char * filename) override
	{
		const uint idxString = 6;
		generate(
			idxString-1,
			calcFreqStringFret(idxString, 3),
			1000.0f,
			KarplusStrongSettings::PluckStr,
			KarplusStrongSettings::ClassicFilter,
			KarplusStrongSettings::Damping,
			KarplusStrongSettings::Richness,
			KarplusStrongSettings::NoiseMix,
			KarplusStrongSettings::VibrAmp,
			KarplusStrongSettings::VibrFreq,
			KarplusStrongSettings::TuneGain,
			KarplusStrongSettings::TuneDelay,

			0.0f
			);

		const uint blockSizeChannel = (uint)(m_numSeconds * m_samplesPerSec);
		formats::writeWAV(filename, m_bufInt, blockSizeChannel, m_samplesPerSec, m_numChannels);
	}

	void playString(
		int idxString,
		scalar stringFreq,
		scalar lowPassFreq,
		scalar pluckStrength,
		bool lowOrderFilter,
		scalar damping,
		scalar filterRichness,
		scalar noiseMix,
		scalar vibrAmp,
		scalar vibrFreq,

		scalar tuneGain,
		scalar tuneDelay,

		scalar startMarginMS = 0
		)
	{
		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);
		const uint startMargin = (uint)(startMarginMS * m_samplesPerSec);

		generate(
			idxString,
			stringFreq,
			lowPassFreq,
			pluckStrength,
			lowOrderFilter,
			damping,
			filterRichness,
			noiseMix,
			vibrAmp,
			vibrFreq,

			tuneGain,
			tuneDelay,

			startMarginMS
			);

		KS_PlayBuffer(
			blockSize,
			startMargin,
			idxString
			);
	}

	void playTone(
		int idxString,
		scalar stringFreq
		)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);

		const scalar dutyTime = 0.5f;

		const scalar freq = stringFreq;

		m_perfTimer.start();

		for (uint i = 0; i < blockSize; ++i)
		{
			scalar sineTune = sinf(_2PI * stringFreq * i / (scalar)m_samplesPerSec);
			m_buffer[i*2 + 0] = sineTune;
			m_buffer[i*2 + 1] = sineTune;
		}

		KS_MakeBuffer(
			0,
			blockSize,
			0,
			1.0f
			);

		KS_PlayBuffer(
			blockSize,
			0,
			idxString
			);

		m_execTime = (scalar)m_perfTimer.time();

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
	}
};

//////////////////////////////////////////////////////////////////////////

class KSSynthSequencer
{
protected:

	double m_seqTime;
	uint m_frame;

	SoundKarplusStrong * m_ksGenerator;

	scalar m_recordingSpeed;
	bool m_recordSamples;
	std::vector<float> m_samplesInterleavedRaw;
	std::vector<short> m_samplesInterleaved;

public:

	struct Event
	{
		uint frame;

		uint idxString;
		double time;
		scalar freq;
		scalar richness;
		scalar damp;
		scalar pluckStr;

		bool useGlobaSettings;

		Event(uint idxString_, double time_, scalar freq_, scalar richness_, scalar damp_, scalar pluckStr_):
			idxString(idxString_),
			time(time_),
			freq(freq_),
			richness(richness_),
			damp(damp_),
			pluckStr(pluckStr_),
			frame(0),
			useGlobaSettings(true)
		{
		}
	};

	double m_totalTime;
	std::vector<Event> Events;

	void setSampleRecordingState(bool recordSamples) { m_recordSamples = recordSamples; }
	bool getSampleRecordingState() const { return m_recordSamples; }
	void setSampleRecordingSpeed(scalar recordingSpeed) { m_recordingSpeed = recordingSpeed; }
	scalar getSampleRecordingSpeed() const { return m_recordingSpeed; }

	const short * getBufferRecorded()
	{
		return m_samplesInterleaved.data();
	}
	uint getBufferSizeRecorded()
	{
		return (uint)m_samplesInterleaved.size();
	}


	// 2 channels only
	void recordSamples(short * samples, uint startSample, uint numSamples)
	{
		uint currentSamplesNum = (uint)m_samplesInterleavedRaw.size();
		uint lastRecordedSamplePos = startSample + numSamples;

		if (lastRecordedSamplePos > currentSamplesNum)
		{
			m_samplesInterleavedRaw.resize(lastRecordedSamplePos);
			for (uint i = 0, iEnd = lastRecordedSamplePos - currentSamplesNum; i < iEnd; ++i)
			{
				m_samplesInterleavedRaw[currentSamplesNum + i] = 0;
			}
		}

		// Mixing in samples, to be normalized later
		const scalar mixGain = 1.0f / 32767;
		for (uint i = 0; i < numSamples; ++i)
		{
			m_samplesInterleavedRaw[startSample + i] += mixGain*samples[i];
		}
	}
	void resetSamplesRecording()
	{
		m_samplesInterleavedRaw.resize(0);
	}
	// Normalize and convert to short
	void prepareRecording()
	{
		const size_t numSamplesRecorded = m_samplesInterleavedRaw.size();

		scalar maxVolume = 1.0f;
		for (size_t i = 0; i < numSamplesRecorded; ++i)
		{
			float absSample = fast_abs(m_samplesInterleavedRaw[i]);
			if (maxVolume < absSample)
			{
				maxVolume = absSample;
			}
		}
		scalar invMaxVolume = 0.95f / maxVolume;
		m_samplesInterleaved.resize(numSamplesRecorded);
		for (size_t i = 0; i < numSamplesRecorded; ++i)
		{
			m_samplesInterleaved[i] = (short)(m_samplesInterleavedRaw[i] * invMaxVolume * 32767);
		}
	}

	KSSynthSequencer(SoundKarplusStrong * ksGenerator):
		m_ksGenerator(ksGenerator),
		m_frame(1),
		m_seqTime(0.0),
		m_totalTime(10.0f),
		m_recordingSpeed(1.0f),
		m_recordSamples(false)
	{
		Events.reserve(1000);
	}

	void update(double dt)
	{
		m_seqTime += dt;
		if (m_seqTime > m_totalTime)
		{
			m_seqTime -= m_totalTime;
			++m_frame;
		}

		const scalar lookAheadTime = 0.0f;

		for (uint i = 0, iend = (uint)Events.size(); i < iend; ++i)
		{
			Event & curEvent = Events[i];
			if (curEvent.time < (m_seqTime + lookAheadTime) && curEvent.frame != m_frame)
			{
				curEvent.frame = m_frame;

				scalar soundDelay = 0.0f;//curEvent.time - mSeqTime;
				//if (soundDelay < 0.0f)
				//	soundDelay = 0.0f;

				// Play event
				if (curEvent.useGlobaSettings)
				{
					m_ksGenerator->playString(curEvent.idxString, curEvent.freq, 1000.0f, KarplusStrongSettings::PluckStr,
						KarplusStrongSettings::ClassicFilter, KarplusStrongSettings::Damping, KarplusStrongSettings::Richness, KarplusStrongSettings::NoiseMix,
						KarplusStrongSettings::VibrAmp, KarplusStrongSettings::VibrFreq, KarplusStrongSettings::TuneGain, KarplusStrongSettings::TuneDelay,
						soundDelay
						);
				}
				else
				{
					m_ksGenerator->playString(curEvent.idxString, curEvent.freq, 1000.0f, curEvent.pluckStr,
						KarplusStrongSettings::ClassicFilter, curEvent.damp, curEvent.richness, KarplusStrongSettings::NoiseMix,
						KarplusStrongSettings::VibrAmp, KarplusStrongSettings::VibrFreq, KarplusStrongSettings::TuneGain, KarplusStrongSettings::TuneDelay);
				}

				if (m_recordSamples)
				{
					uint totalSamplesNum = m_ksGenerator->getBufferSizeGenerated();
					short * samples = m_ksGenerator->getBufferGenerated();

					// The sequence could be played more than once, hence the totalTime addition
					// Recording at ideal event time, trather than at current time (m_seqTime)
					// so that recording won't have uneven timings if system hitches happened
					double sampleTimeOffset = curEvent.time + (m_frame-1)*m_totalTime;
					recordSamples(samples, (uint)(sampleTimeOffset/m_recordingSpeed*m_ksGenerator->getSamplesPerSec())*m_ksGenerator->getNumChannels(), totalSamplesNum);
				}
			}
		}
	}
};

void InitKS_Horizons(KSSynthSequencer * sequencer)
{
	double time = 0.5;

#define KS_ADD_STRING_EVENT(idxString, idxFret, timeOffset) \
	time += timeOffset; \
	sequencer->Events.push_back(KSSynthSequencer::Event(idxString-1, time, calcFreqStringFret(idxString, idxFret), 0.5f, 0.001f, 0.01f));

	// Horizons

#if 1
	
	for (uint i = 0; i < 2; ++i)
	{
		KS_ADD_STRING_EVENT(6, 3, 0.3);

		KS_ADD_STRING_EVENT(4, 0, 0.3);
		KS_ADD_STRING_EVENT(2, 0, 0.3);
		KS_ADD_STRING_EVENT(3, 2, 0.3);
		KS_ADD_STRING_EVENT(2, 0, 0.3);

		KS_ADD_STRING_EVENT(4, 0, 0.3);
		KS_ADD_STRING_EVENT(2, 0, 0.3);
		KS_ADD_STRING_EVENT(3, 0, 0.3);
	}
	//*/

	const double simultStringPluck = 0.05;
	//////////////////////////////////////////////////////////////////////////
	
	KS_ADD_STRING_EVENT(5, 0, 0.3);
	KS_ADD_STRING_EVENT(4, 2, 0.3);
	KS_ADD_STRING_EVENT(3, 2, 0.3);
	KS_ADD_STRING_EVENT(2, 1, simultStringPluck);
	KS_ADD_STRING_EVENT(2, 0, 0.3-simultStringPluck);
	KS_ADD_STRING_EVENT(2, 1, 0.3);
	KS_ADD_STRING_EVENT(6, 3, simultStringPluck);

	KS_ADD_STRING_EVENT(4, 2, 0.35-simultStringPluck);
	KS_ADD_STRING_EVENT(2, 1, 0.3);
	KS_ADD_STRING_EVENT(4, 0, 0.3);

	//////////////////////////////////////////////////////////////////////////

	KS_ADD_STRING_EVENT(6, 2, 0.3);
	KS_ADD_STRING_EVENT(4, 0, 0.3);
	KS_ADD_STRING_EVENT(3, 2, 0.3);
	KS_ADD_STRING_EVENT(2, 1, simultStringPluck);
	KS_ADD_STRING_EVENT(2, 0, 0.3-simultStringPluck);
	KS_ADD_STRING_EVENT(2, 0, 0.3);
	KS_ADD_STRING_EVENT(6, 3, simultStringPluck);

	KS_ADD_STRING_EVENT(4, 0, 0.35-simultStringPluck);
	KS_ADD_STRING_EVENT(3, 2, 0.3);
	KS_ADD_STRING_EVENT(4, 0, 0.3);
	KS_ADD_STRING_EVENT(2, 0, 0.3);
	KS_ADD_STRING_EVENT(4, 0, 0.3);
	KS_ADD_STRING_EVENT(1, 0, 0.3);
	KS_ADD_STRING_EVENT(4, 0, 0.3);

	KS_ADD_STRING_EVENT(6, 2, 0.3);
	KS_ADD_STRING_EVENT(2, 3, simultStringPluck);
	//*/

	//////////////////////////////////////////////////////////////////////////

#endif

	KS_ADD_STRING_EVENT(4, 0, 0.3-simultStringPluck);
	KS_ADD_STRING_EVENT(3, 2, 0.3);
	KS_ADD_STRING_EVENT(2, 0, 0.3);

#if 1
	KS_ADD_STRING_EVENT(4, 0, 0.3);
	KS_ADD_STRING_EVENT(3, 2, simultStringPluck);
#else
	KS_ADD_STRING_EVENT(3, 2, 0.3-simultStringPluck);
#endif
	KS_ADD_STRING_EVENT(2, 1, simultStringPluck);

#if 1
	KS_ADD_STRING_EVENT(4, 0, 0.6);
	KS_ADD_STRING_EVENT(3, 0, simultStringPluck);
#else
	KS_ADD_STRING_EVENT(3, 0, 0.6-simultStringPluck);
#endif
	KS_ADD_STRING_EVENT(2, 0, simultStringPluck);

	time += 1.0f;

#undef KS_ADD_STRING_EVENT

	sequencer->m_totalTime = time;
}

void InitKS_Strumming(KSSynthSequencer * sequencer)
{
	//0		82.4	110.0	146.8	196.0	246.9	329.6

	scalar stBaseFreqs[6] =
	{
		329.6f,
		246.9f,
		196.0f,
		146.8f,
		110.0f,
		82.4f
	};

	scalar stBaseSemitones[6];

	for (uint i = 0; i < 6; ++i)
	{
		stBaseSemitones[i] = 12 * m_log2(stBaseFreqs[i] / 440.0f) + 49;
	}

	int fretsAm[6] =
	{
		0,
		1,
		2,
		2,
		0,
		0
	};
	int fretsC[6] =
	{
		0,
		1,
		0,
		2,
		3,
		0
	};
	int fretsEm[6] =
	{
		0,
		0,
		0,
		2,
		2,
		0
	};


	const uint numStrums = 4;

	int * chords[numStrums] =
	{
		fretsAm,
		fretsC,
		fretsEm,
		fretsEm
	};

	double strumStringTime = 0.015, strumStringVariation = 0.005f;
	double time = 0.5;

	// Strumming
	for (uint idxStrum = 0; idxStrum < numStrums; ++idxStrum)
	{
		int * fretsIdxChord = chords[idxStrum];
		for (uint i = 0; i < 6; ++i)
		{
			// DOWN
			int idxString = 5-i;
			scalar freq = calcFretFreq(stBaseSemitones[idxString], fretsIdxChord[idxString]);
			sequencer->Events.push_back(KSSynthSequencer::Event(idxString, time, freq, 0.5f, 0.001f, 0.01f));
			time += strumStringTime + strumStringVariation*(2.0f*randomN());
		}

		time += 0.475;

		for (uint i = 0; i < 6; ++i)
		{
			// DOWN
			int idxString = 5-i;
			scalar freq = calcFretFreq(stBaseSemitones[idxString], fretsIdxChord[idxString]);
			//			sequencer->Events.push_back(KSSynthSequencer::Event(idxString, time, freq, 0.1f, 0.05f, 0.1f));
			KSSynthSequencer::Event downEvent(idxString, time, freq, KarplusStrongSettings::Richness, 0.03f, KarplusStrongSettings::PluckStr);
			downEvent.useGlobaSettings = false;
			sequencer->Events.push_back(downEvent);
			time += strumStringTime + strumStringVariation*(2.0f*randomN());
		}

		//time += 0.17;
		time += 0.05;

		for (uint i = 0; i < 6; ++i)
		{
			// UP
			int idxString = i;
			scalar freq = calcFretFreq(stBaseSemitones[idxString], fretsIdxChord[idxString]);
			sequencer->Events.push_back(KSSynthSequencer::Event(idxString, time, freq, 0.5f, 0.001f, 0.01f));
			time += strumStringTime + strumStringVariation*(2.0f*randomN());
		}

		//time += 0.17;
		time += 0.29;

		for (uint i = 0; i < 6; ++i)
		{
			// DOWN
			int idxString = 5-i;
			scalar freq = calcFretFreq(stBaseSemitones[idxString], fretsIdxChord[idxString]);
			//			sequencer->Events.push_back(KSSynthSequencer::Event(idxString, time, freq, 0.1f, 0.01f, 0.1f));
			sequencer->Events.push_back(KSSynthSequencer::Event(idxString, time, freq, 0.5f, 0.001f, 0.01f));
			time += strumStringTime + strumStringVariation*(2.0f*randomN());
		}

		time += 1.0;
	}

	// Last strum added +1.0, and we have 0.5 at the beginning(+1.5 total), so to have 1.0 we need to do -0.5
	time -= 0.5;

	sequencer->m_totalTime = time;
}

#undef PARAM_DECLARE
