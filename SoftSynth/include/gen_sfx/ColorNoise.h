#pragma once

#include <xmmintrin.h>

#include "math/AuxMath.h"
#include "windows/Timer.h"
#include "formats/wav.h"

#include "windows/WaveOut_Simple.h"
#include "sfx/Oscillators.h"
#include "sfx/Filters.h"
#include "sfx/Reverb.h"
#include "sfx/Phaser.h"
#include "sfx/Flanger.h"

#include "sfx/PinkNoise.h"
#include "Base.h"

#define PARAM_DECLARE(param, val, min, max) \
	scalar param = val; \
	const scalar param##Min = min; \
	const scalar param##Max = max;

namespace ColorNoiseSettings
{
	PARAM_DECLARE(NoiseAlpha, 1.5f, 0.0f, 2.0f);
	PARAM_DECLARE(SawMix, 0.0f, 0.0f, 1.0f);
}

class SoundColorNoise : public SoundPlay
{
	static const int c_numDevices = 2;
	static const int c_samplesPerSec = 44100;
	static const int c_numChannels = 2;

	size_t m_bufSize;
	float * m_buffer;
	short * m_bufInt;

	sfx::PinkNoise m_pinkNoise;

	bool m_isGenerated = false;

public:

	enum ColorNoiseType
	{
		eAlpha,		// Adjustable
		eOctave,
		eCustom
	};
	ColorNoiseType colorNoiseType;

	SoundColorNoise(scalar numSeconds = 2.0f):
		SoundPlay(c_numDevices, c_samplesPerSec, c_numChannels),
		colorNoiseType(ColorNoiseType::eAlpha)
	{
		m_numSeconds = numSeconds;

		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec * m_numChannels);

		m_bufSize = toneBlockSize;
		m_buffer = new float [m_bufSize];
		m_bufInt = new short [m_bufSize];
	}

	~SoundColorNoise()
	{
		delete [] m_buffer;
		delete [] m_bufInt;
	}

	virtual void generate() override
	{
		using namespace sfx;

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);

		m_perfTimer.start();

		scalar invSPC = 1.0f / (scalar)m_samplesPerSec;

		// Comparison
		if (colorNoiseType == ColorNoiseType::eAlpha)
		{
			m_pinkNoise.setAlpha(ColorNoiseSettings::NoiseAlpha);
			m_pinkNoise.init();

			scalar colNoiseMax = -FP_MAX;
			for (uint i = 0, iend = blockSize; i < iend; ++i)
			{
				float noise = m_pinkNoise.generateSample();
				float saw = oscTriangle(80.0f * i / (scalar)m_samplesPerSec, 0.99f);

				noise = noise * (1.0f - ColorNoiseSettings::SawMix) + saw * ColorNoiseSettings::SawMix;

				scalar absNoise = fast_abs(noise);
				if (colNoiseMax < absNoise)
					colNoiseMax = absNoise;

				m_buffer[i*m_numChannels] = noise;
			}

			for (uint i = 0, iend = blockSize; i < iend; ++i)
			{
				uint rnd30 = rand30() & 65535;
				//float noise = (rnd30 - 32768.0f) / 32768.0f;
				float noise = m_buffer[i*m_numChannels] / (colNoiseMax * 1.1f);

				m_buffer[i*m_numChannels + 0] = noise;
				m_buffer[i*m_numChannels + 1] = noise;
			}
		}
		else if (colorNoiseType == ColorNoiseType::eOctave)
		{
			const uint numOctaves = 6;

			LowPassIIR lowPasses[numOctaves];
			scalar freq = 200.0f / m_samplesPerSec;
			for (uint i = 0; i < numOctaves; ++i)
			{
				lowPasses[i].setParams(freq);
				lowPasses[i].autoNormalize();

				freq *= 2.0f;
			}

			scalar max = 0.0f;
			for (uint i = 0; i < blockSize; ++i)
			{
				scalar snd = 0.0f, amp = 0.5f;
				for (uint j = 0; j < numOctaves; ++j)
				{
					uint rnd30 = rand30() & 65535;
					scalar curOctave, whiteNoise = (rnd30 - 32767.5f) / 32767.5f;
					lowPasses[j].applyCont(&curOctave, whiteNoise);
					snd += curOctave * amp;
					amp *= 0.5f;
				}

				// Amplify a bit
				snd /= 0.7f;

				// Simple wavemapper
				//snd = 1.5f * snd - 0.5f * snd * snd * snd;

				scalar gain = 1.0f;
				if (i > blockSize/2)
				{
					gain -= (float)(i - blockSize/2) / (blockSize/2); 
				}

				m_buffer[i*2 + 0] = gain * snd;
				m_buffer[i*2 + 1] = gain * snd;
			}
		}
		else if (colorNoiseType == ColorNoiseType::eCustom)
		{
			const uint numStages = 3;

			const scalar A[] = { 0.02109238f, 0.07113478f, 0.68873558f }; // rescaled by (1+P)/(1-P)
			const scalar P[] = { 0.3190f,  0.7756f,  0.9613f };

			scalar state[numStages];

			for (uint i = 0; i < numStages; ++i)
			{
				state[i] = 0.0f;
			}

			const scalar RMI2 = 7.0f / 65536.0f;

			uint curHistoryIdx = 0;
			for (uint i = 0; i < blockSize; ++i)
			{
				int noise;

				noise = (rand30() & 65535) - 32768;
				state[0] = P[0] * (state[0] - noise) + noise;
				noise = (rand30() & 65535) - 32768;
				state[1] = P[1] * (state[1] - noise) + noise;
				noise = (rand30() & 65535) - 32768;
				state[2] = P[2] * (state[2] - noise) + noise;

				scalar pink_noise = (A[0]*state[0] + A[1]*state[1] + A[2]*state[2]) * RMI2;
				scalar snd = pink_noise;

				m_buffer[i*2 + 0] = snd;
				m_buffer[i*2 + 1] = snd;
			}
		}

		scalar maxVolume = 1.0f;

		// Convert float buf into shorts
		const scalar volume = 0.9f / maxVolume;
		for (uint j = 0; j < blockSize * m_numChannels; ++j)
		{
			m_bufInt[j] = int(volume * clamp(m_buffer[j], -1.0f, 1.0f) * 32767);
		}

		m_isGenerated = true;
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

