#pragma once

#include <xmmintrin.h>

#include "windows/Timer.h"
#include "formats/wav.h"

#include "windows/WaveOut_Simple.h"
#include "SFX/Oscillators.h"
#include "SFX/Filters.h"
#include "SFX/Reverb.h"

#include "Base.h"

#define FILTER_BREATH_NOISE	1
#define FILTER_FINAL		1

#define PARAM_DECLARE(param, val, min, max) \
	scalar param = val; \
	const scalar param##Min = min; \
	const scalar param##Max = max;

namespace ReverbSettings
{
	int FeedBackMatrix = 6;

	bool Enable = false;

	PARAM_DECLARE(Delay0Len,  773.0f, 100.0f, 10000.0f);
	PARAM_DECLARE(Delay1Len, 1709.0f, 100.0f, 10000.0f);
	PARAM_DECLARE(Delay2Len, 3001.0f, 100.0f, 10000.0f);
	PARAM_DECLARE(Delay3Len, 4367.0f, 100.0f, 10000.0f);

	PARAM_DECLARE(Halflife0, 9000.0f, 100.0f, 50000.0f);
	PARAM_DECLARE(Halflife1, 3000.0f, 100.0f, 50000.0f);

	PARAM_DECLARE(InOutRot0, 0.0f, -300.0f, 300.0f);
	PARAM_DECLARE(InOutRot1, 0.0f, -300.0f, 300.0f);

	PARAM_DECLARE(Modulation0, 20000.0f, 5000.0f, 100000.0f);
	PARAM_DECLARE(Modulation1, 0.3f, 0.0f, 1.0f);

	PARAM_DECLARE(Mix0, 0.4f, 0.0f, 1.0f);
	PARAM_DECLARE(Mix1, 0.6f, 0.0f, 1.0f);
}

namespace FormantSettings
{
	int Times = 4;

	namespace BandPass0Seq
	{
		PARAM_DECLARE(Freq, 270.0f, 0.0f, 7000.0f);
		PARAM_DECLARE(Q, 0.98f, 0.0f, 1.0f);
		PARAM_DECLARE(Vol, 1.0f, 0.0f, 1.0f);
	}
	namespace BandPass1Seq
	{
		PARAM_DECLARE(Freq, 2300.0f, 0.0f, 7000.0f);
		PARAM_DECLARE(Q, 0.96f, 0.0f, 1.0f);
		PARAM_DECLARE(Vol, 0.2f, 0.0f, 1.0f);
	}
	namespace BandPass2Seq
	{
		PARAM_DECLARE(Freq, 3000.0f, 0.0f, 7000.0f);
		PARAM_DECLARE(Q, 0.99f, 0.0f, 1.0f);
		PARAM_DECLARE(Vol, 0.4f, 0.0f, 1.0f);
	}

	PARAM_DECLARE(BaseFreq, 120.0f, 10.0f, 500.0f);
	PARAM_DECLARE(RandAmp, 0.1f, 0.0f, 1.0f);

	PARAM_DECLARE(Duration, 0.8f, 0.0f, 1.0f);

	PARAM_DECLARE(VibrFreq, 4.0f, 0.1f, 50.0f);
	PARAM_DECLARE(VibrAmp, 20.0f, 1.0f, 300.0f);
}

class SoundFormant : public SoundPlay
{
	static const int c_numDevices = 5;
	static const int c_samplesPerSec = 44100;
	static const int c_numChannels = 2;

	size_t m_bufSize;
	float * m_buffer, * m_auxBuf0, * m_auxBuf1, * m_auxBuf2;
	short * m_bufInt;

	sfx::BandPassIIR bandPass0, bandPass1, bandPass2;
	sfx::LowPassIIR lowPass;

	// Reverb data
	sfx::ReverbFDN4 m_reverb;
	uint m_reverbDelaySizes[4];

	bool m_isGenerated = false;

public:

	int genType;

	SoundFormant(scalar numSeconds = 2.0f):
		SoundPlay(c_numDevices, c_samplesPerSec, c_numChannels),
		genType(1)
	{
		m_numSeconds = numSeconds;

		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec * m_numChannels);

		m_bufSize = toneBlockSize;
		m_buffer = new float [m_bufSize];
		m_auxBuf0 = new float [m_bufSize];
		m_auxBuf1 = new float [m_bufSize];
		m_auxBuf2 = new float [m_bufSize];
		m_bufInt = new short [m_bufSize];

		for (uint i = 0; i < 4; ++i)
			m_reverbDelaySizes[i] = 0;
		m_reverb.init();
	}

	~SoundFormant()
	{
		delete [] m_buffer;
		delete [] m_auxBuf0;
		delete [] m_auxBuf1;
		delete [] m_auxBuf2;
		delete [] m_bufInt;
	}

	virtual void generate() override
	{
		using namespace sfx;

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);

		const scalar dutyTime = 0.5f;

		scalar freq = FormantSettings::BaseFreq;

		lowPass.setParams(freq);
		lowPass.autoNormalize(0.0f*freq);

		m_perfTimer.start();

		// Init reverb
		//////////////////////////////////////////////////////////////////////////
		if (ReverbSettings::Enable)
		{
			uint newReverbDelays[4] =
				{
					(uint)ReverbSettings::Delay0Len,
					(uint)ReverbSettings::Delay1Len,
					(uint)ReverbSettings::Delay2Len,
					(uint)ReverbSettings::Delay3Len
				};
			bool delaySizeChanged = false;
			for (uint i = 0; i < 4; ++i)
			{
				if (newReverbDelays[i] != m_reverbDelaySizes[i])
				{
					delaySizeChanged = true;
					m_reverbDelaySizes[i] = newReverbDelays[i];
				}
			}
			if (delaySizeChanged)
			{
				m_reverb.initDLs(m_reverbDelaySizes);
			}

			m_reverb.setHalflives(ReverbSettings::Halflife0, ReverbSettings::Halflife1);
			m_reverb.setFeedbackMatrix((ReverbMatrixTypes::Values)ReverbSettings::FeedBackMatrix);
			m_reverb.setInOutRotations(ReverbSettings::InOutRot0, ReverbSettings::InOutRot1);
			m_reverb.setModulation(ReverbSettings::Modulation0, ReverbSettings::Modulation1);
			m_reverb.setMixCoeffs(ReverbSettings::Mix0, ReverbSettings::Mix1);

			m_reverb.clearDLs();
		}
		//////////////////////////////////////////////////////////////////////////

		scalar phase0 = 0.0f, phase1 = 1.0f, phase2 = 2.0f;

		scalar invSPC = 1.0f / (scalar)m_samplesPerSec;

		scalar noisePhase = 0.0f;
		const int noiseStretchCoeff = 20;
		int noiseCnt = 0;

		uint samplesProcessed = 0;
		for (uint i = 0, iend = (uint)(blockSize * FormantSettings::Duration); i < iend; ++i)
		{
			scalar phase = i / (scalar)m_samplesPerSec;

			scalar vibrato = 0.2f * FormantSettings::VibrAmp * sinf(FormantSettings::VibrFreq * _2PI * phase);

			float snd;
			scalar val;
			switch (genType)
			{
			case 0:
				{
					// Sine
					val = oscSine(freq, 1.0f, 0.0f, phase, dutyTime);
					break;
				}
			case 1:
				{
					// Saw
					val = 0.5f * oscTriangle(phase0, 0.999f) +
						  0.3f * oscTriangle(phase1, 0.001f) +
						  0.2f * oscTriangle(phase2, 0.001f);

					scalar noise = (rand()%16383) / 8192.0f - 1.0f;

#if (FILTER_BREATH_NOISE == 1)
					scalar noiseFiltered;
					lowPass.applyCont(&noiseFiltered, noise);
					val = 0.8f * val + 0.2f * noiseFiltered;
#else
					val = 0.8f * val + 0.2f * noise;
#endif

					break;
				}
			case 2:
				{
					// Sq
					val = oscSquare(freq, 1.0f, 0.0f, phase, dutyTime);
					break;
				}
			case 3:
				{
					// Tri
					val = oscTriangle(freq, 1.0f, 0.0f, phase, dutyTime);
					break;
				}
			case 4:
				{
					// White noise
					val = (rand()%16383) / 8192.0f - 1.0f;
					break;
				}
			};

			if (noiseCnt > noiseStretchCoeff)
				noiseCnt = 0;

			if (noiseCnt == 0)
				noisePhase = 1.0f + FormantSettings::RandAmp * ((rand()%16383) / 8192.0f - 1.0f);
			++noiseCnt;

			phase0 += invSPC * noisePhase * (freq + vibrato);
			phase1 += invSPC * noisePhase * (freq + vibrato) * 1.01f;
			phase2 += invSPC * noisePhase * (freq + vibrato) * 1.05f;

			snd = val;

			m_buffer[i*m_numChannels + 0] = snd;
			m_buffer[i*m_numChannels + 1] = snd;

			samplesProcessed += m_numChannels;
		}
		for (uint i = samplesProcessed; i < blockSize*m_numChannels; ++i)
		{
			m_buffer[i] = 0.0f;
		}

		lowPass.resetHistory();
		lowPass.setParams(0.05f);
		lowPass.autoNormalize(0.0f);

		scalar maxVolume = 1.0f;

		bool soundUseFiltering = true;
		if (soundUseFiltering)
		{
			// Copy the initial signal into aux buffers
			for (uint i = 0; i < blockSize; ++i)
			{
				int idxL = i*m_numChannels + 0, idxR = i*m_numChannels + 1;

				float curSampleL = m_buffer[idxL];
				float curSampleR = m_buffer[idxR];

				m_auxBuf0[idxL] = curSampleL;
				m_auxBuf0[idxR] = curSampleR;

				m_auxBuf1[idxL] = curSampleL;
				m_auxBuf1[idxR] = curSampleR;

				m_auxBuf2[idxL] = curSampleL;
				m_auxBuf2[idxR] = curSampleR;
			}

			// Formant filter the buffer
			scalar bp0Freq = FormantSettings::BandPass0Seq::Freq / (scalar)m_samplesPerSec;
			bandPass0.setParams(bp0Freq, FormantSettings::BandPass0Seq::Q);
			bandPass0.autoNormalize(bp0Freq);

			scalar bp1Freq = FormantSettings::BandPass1Seq::Freq / (scalar)m_samplesPerSec;
			bandPass1.setParams(bp1Freq, FormantSettings::BandPass1Seq::Q);
			bandPass1.autoNormalize(bp1Freq);

			scalar bp2Freq = FormantSettings::BandPass2Seq::Freq / (scalar)m_samplesPerSec;
			bandPass2.setParams(bp2Freq, FormantSettings::BandPass2Seq::Q);
			bandPass2.autoNormalize(bp2Freq);

			const int times = FormantSettings::Times;

			const int numFormants = 3;
			const int bpTapsI = 3;
			const int bpTapsO = 2;

			const int historyBufSize = numFormants * ( times * (bpTapsI + bpTapsO) );

			float * historyBuf = new float [ historyBufSize ];
			float * historyBufPtr = historyBuf;

			float * history0I = historyBufPtr;
			historyBufPtr += bpTapsI * times;
			float * history0O = historyBufPtr;
			historyBufPtr += bpTapsO * times;

			float * history1I = historyBufPtr;
			historyBufPtr += bpTapsI * times;
			float * history1O = historyBufPtr;
			historyBufPtr += bpTapsO * times;

			float * history2I = historyBufPtr;
			historyBufPtr += bpTapsI * times;
			float * history2O = historyBufPtr;
			historyBufPtr += bpTapsO * times;

			int DBGhistoryBufAllocated = (int)(historyBufPtr - historyBuf);

			// Zero history
			for (int hi = 0; hi < historyBufSize; ++hi)
				historyBuf[hi] = 0.0f;

			float * bp0HistI = bandPass0.getHistoryI();
			float * bp0HistO = bandPass0.getHistoryO();
			float * bp1HistI = bandPass1.getHistoryI();
			float * bp1HistO = bandPass1.getHistoryO();
			float * bp2HistI = bandPass2.getHistoryI();
			float * bp2HistO = bandPass2.getHistoryO();

			// For now, assume mono signal (l == r)
			for (uint i = 0; i < blockSize; ++i)
			{
				int idxL = i*m_numChannels + 0, idxR = i*m_numChannels + 1;
				float bpResult;

				scalar vibrato = FormantSettings::VibrAmp * sinf(FormantSettings::VibrFreq * _2PI * i / (scalar)m_samplesPerSec);
				scalar bp0Freq = (FormantSettings::BandPass0Seq::Freq + vibrato) / (scalar)m_samplesPerSec;
				bandPass0.setParams(bp0Freq, FormantSettings::BandPass0Seq::Q);
				bandPass0.autoNormalize(bp0Freq);

				scalar bp1Freq = (FormantSettings::BandPass1Seq::Freq + vibrato) / (scalar)m_samplesPerSec;
				bandPass1.setParams(bp1Freq, FormantSettings::BandPass1Seq::Q);
				bandPass1.autoNormalize(bp1Freq);

				for (int ti = 0; ti < times; ++ti)
				{
					// Load filters history for current filter #times
					for (uint taps = 0; taps < bpTapsI; ++taps)
					{
						bp0HistI[taps] = history0I[ti*bpTapsI + taps];
						bp1HistI[taps] = history1I[ti*bpTapsI + taps];
						bp2HistI[taps] = history2I[ti*bpTapsI + taps];
					}
					for (uint taps = 0; taps < bpTapsO; ++taps)
					{
						bp0HistO[taps] = history0O[ti*bpTapsO + taps];
						bp1HistO[taps] = history1O[ti*bpTapsO + taps];
						bp2HistO[taps] = history2O[ti*bpTapsO + taps];
					}

					bandPass0.applyCont(&bpResult, m_auxBuf0[idxL]);
					m_auxBuf0[idxL] = bpResult;
					m_auxBuf0[idxR] = bpResult;

					bandPass1.applyCont(&bpResult, m_auxBuf1[idxL]);
					m_auxBuf1[idxL] = bpResult;
					m_auxBuf1[idxR] = bpResult;

					bandPass2.applyCont(&bpResult, m_auxBuf2[idxL]);
					m_auxBuf2[idxL] = bpResult;
					m_auxBuf2[idxR] = bpResult;

					// Save filters history for current filter #times
					for (uint taps = 0; taps < bpTapsI; ++taps)
					{
						//if (bp0HistI[taps] < 1e-5f && bp0HistI[taps] > -1e-5f) bp0HistI[taps] = 0.0f;
						//if (bp1HistI[taps] < 1e-5f && bp1HistI[taps] > -1e-5f) bp1HistI[taps] = 0.0f;
						//if (bp2HistI[taps] < 1e-5f && bp2HistI[taps] > -1e-5f) bp2HistI[taps] = 0.0f;

						history0I[ti*bpTapsI + taps] = bp0HistI[taps];
						history1I[ti*bpTapsI + taps] = bp1HistI[taps];
						history2I[ti*bpTapsI + taps] = bp2HistI[taps];
					}
					for (uint taps = 0; taps < bpTapsO; ++taps)
					{
						//if (bp0HistO[taps] < 1e-5f && bp0HistO[taps] > -1e-5f) bp0HistO[taps] = 0.0f;
						//if (bp1HistO[taps] < 1e-5f && bp1HistO[taps] > -1e-5f) bp1HistO[taps] = 0.0f;
						//if (bp2HistO[taps] < 1e-5f && bp2HistO[taps] > -1e-5f) bp2HistO[taps] = 0.0f;

						history0O[ti*bpTapsO + taps] = bp0HistO[taps];
						history1O[ti*bpTapsO + taps] = bp1HistO[taps];
						history2O[ti*bpTapsO + taps] = bp2HistO[taps];
					}
				}
			}

			delete [] historyBuf;

			maxVolume = -FP_MAX;
			// Combine (mix) all the three BP paths together
			for (uint i = 0; i < blockSize; ++i)
			{
				// Assuming mono again
				scalar mix =
					FormantSettings::BandPass0Seq::Vol * m_auxBuf0[i*m_numChannels + 0] +
					FormantSettings::BandPass1Seq::Vol * m_auxBuf1[i*m_numChannels + 0] +
					FormantSettings::BandPass2Seq::Vol * m_auxBuf2[i*m_numChannels + 0];

#if (FILTER_FINAL == 1)
				scalar mixFiltered;
				lowPass.applyCont(&mixFiltered, mix);
				m_buffer[i*m_numChannels + 0] = mixFiltered;
				m_buffer[i*m_numChannels + 1] = mixFiltered;
#else
				m_Buffer[i*m_NumChannels + 0] = mix;
				m_Buffer[i*m_NumChannels + 1] = mix;
#endif

				if (maxVolume < m_buffer[i*m_numChannels + 0])
					maxVolume = m_buffer[i*m_numChannels + 0];
				if (maxVolume < m_buffer[i*m_numChannels + 1])
					maxVolume = m_buffer[i*m_numChannels + 1];
			}
		}

		if (ReverbSettings::Enable)
		{
			m_reverb.processStereoInterleaved(m_buffer, blockSize);
		}

		// Convert float buf into shorts
		const scalar volume = 0.9f / maxVolume;
		for (uint j = 0; j < blockSize * m_numChannels; ++j)
		{
			m_bufInt[j] = int(volume * clamp( m_buffer[j], -1.0f, 1.0f ) * 32767);
		}
	}
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

	virtual void playSound() override
	{
		if (!m_isGenerated)
		{
			generate();
		}

		const uint blockSizeChannel = (uint)(m_numSeconds * m_samplesPerSec);
		playSoundBuffer(m_bufInt, blockSizeChannel);
	}
	virtual void storeSound(const char * filename) override
	{
		if (!m_isGenerated)
		{
			generate();
		}

		const uint blockSizeChannel = (uint)(m_numSeconds * m_samplesPerSec);
		formats::writeWAV(filename, m_bufInt, blockSizeChannel, m_samplesPerSec, m_numChannels);
	}
};

#undef PARAM_DECLARE
