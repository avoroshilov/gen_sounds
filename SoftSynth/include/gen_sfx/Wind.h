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

namespace WindSettings
{
	bool isFlanger = false;

	PARAM_DECLARE(NoiseAlpha, 1.4f, 0.0f, 2.0f);
	PARAM_DECLARE(SawMix, 0.0f, 0.0f, 1.0f);

	PARAM_DECLARE(Feedback, 0.9f, 0.0f, 1.0f);
	PARAM_DECLARE(Rate, 0.05f, 0.0f, 5.0f);
	PARAM_DECLARE(Phase, 0.5f, 0.0f, 1.0f);
	PARAM_DECLARE(Depth, 1.0f, 0.0f, 1.0f);
	PARAM_DECLARE(NumAPs, 6.0f, 1.0f, 30.0f);
	PARAM_DECLARE(RangeMin, 440.0f, 0.0f, 5000.0f);
	PARAM_DECLARE(RangeMax, 1600.0f, 0.0f, 5000.0f);
}

class SoundWind : public SoundPlay
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

	int genType;

	SoundWind(scalar numSeconds = 2.0f):
		SoundPlay(c_numDevices, c_samplesPerSec, c_numChannels),
		genType(1)
	{
		m_numSeconds = numSeconds;

		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec * m_numChannels);

		m_bufSize = toneBlockSize;
		m_buffer = new float [m_bufSize];
		m_bufInt = new short [m_bufSize];
	}

	~SoundWind()
	{
		delete [] m_buffer;
		delete [] m_bufInt;
	}

	virtual void generate() override
	{
		using namespace sfx;

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);

		const scalar dutyTime = 0.5f;
		const scalar freq = 440.0f;

		m_perfTimer.start();

		scalar invSPC = 1.0f / (scalar)m_samplesPerSec;

		m_pinkNoise.setAlpha(WindSettings::NoiseAlpha);
		m_pinkNoise.init();

		Phaser tstPhaser;
		Flanger tstFlanger;

		if (WindSettings::isFlanger)
		{
			tstFlanger.setFeedback(WindSettings::Feedback);

			tstFlanger.setRate(WindSettings::Rate);
			tstFlanger.setPhase(WindSettings::Phase);
			tstFlanger.setDepth(WindSettings::Depth);
			tstFlanger.setRange(WindSettings::RangeMin, WindSettings::RangeMax);

			tstFlanger.init();
		}
		else
		{

#if 0
			// VERY interesting effect
			tstPhaser.SetFeedback(0.9f);
			tstPhaser.SetAllpassNum(24);
			tstPhaser.SetRate(0.05f);
#else
			tstPhaser.setFeedback(WindSettings::Feedback);
			tstPhaser.setAllpassNum((int)(WindSettings::NumAPs + 0.5f));
			tstPhaser.setRate(WindSettings::Rate);
#endif
			tstPhaser.setPhase(WindSettings::Phase);
			tstPhaser.setDepth(WindSettings::Depth);
			tstPhaser.setRange(WindSettings::RangeMin, WindSettings::RangeMax);

			tstPhaser.init();
		}

		scalar colNoiseMax = -FP_MAX;
		for (uint i = 0, iend = blockSize; i < iend; ++i)
		{
			float noise = m_pinkNoise.generateSample();
			float saw = oscTriangle(80.0f * i / (scalar)m_samplesPerSec, 0.99f);

			noise = noise * (1.0f - WindSettings::SawMix) + saw * WindSettings::SawMix;

			scalar absNoise = fast_abs(noise);
			if (colNoiseMax < absNoise)
				colNoiseMax = absNoise;

			m_buffer[i*m_numChannels] = noise;
		}

		scalar noiseMax = -FP_MAX;
		for (uint i = 0, iend = blockSize; i < iend; ++i)
		{
			uint rnd30 = rand30() & 65535;
			//float noise = (rnd30 - 32768.0f) / 32768.0f;
			float noise = m_buffer[i*m_numChannels] / (colNoiseMax * 1.1f);

			scalar absNoise = fast_abs(noise);
			if (noiseMax < absNoise)
				noiseMax = absNoise;

			scalar snd = 0.0f;
			
			if (WindSettings::isFlanger)
				snd = tstFlanger.update(noise);
			else
				snd = tstPhaser.update(noise);

			m_buffer[i*m_numChannels + 0] = snd;
			m_buffer[i*m_numChannels + 1] = snd;
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

