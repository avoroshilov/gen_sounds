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

class SoundGunShot : public SoundPlay
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

	SoundGunShot(scalar numSeconds = 2.0f):
		SoundPlay(c_numDevices, c_samplesPerSec, c_numChannels)
	{
		m_numSeconds = numSeconds;

		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec * m_numChannels);

		m_bufSize = toneBlockSize;
		m_buffer = new float [m_bufSize];
		m_bufInt = new short [m_bufSize];
	}

	~SoundGunShot()
	{
		delete [] m_buffer;
		delete [] m_bufInt;
	}

	virtual void generate() override
	{
		using namespace sfx;

		const uint blockSize = (uint)(m_numSeconds * m_samplesPerSec);

		m_perfTimer.start();

		scalar maxVolume = 1.0f;

		LowPassIIR lowPass1, lowPass2, lowPass3;

		ReverbFDN4 testFDNReverb;
		testFDNReverb.init();

		bool applyReverb = false;
		if (1)
		{
			applyReverb = true;
			uint delayLengths[] = { 2111, 6121, 6748, 7553 };
			testFDNReverb.initDLs(delayLengths);
			testFDNReverb.setHalflives(14000.0f, 1000.0f);
			testFDNReverb.setFeedbackMatrix(ReverbMatrixTypes::eMediumRotation);
			testFDNReverb.setInOutRotations(-45.0f, -45.0f);
			testFDNReverb.setModulation(70000.0f, 50.5f);
			testFDNReverb.setMixCoeffs(0.4f, 0.6f);
		}

		if (applyReverb)
			testFDNReverb.clearDLs();

		lowPass1.setParams((7000.0f + (rand()%1000)) / m_samplesPerSec);
		lowPass1.autoNormalize();
		lowPass2.setParams((1000.0f + (rand()%500)) / m_samplesPerSec);
		lowPass2.autoNormalize();
		lowPass3.setParams((500.0f + (rand()%100)) / m_samplesPerSec);
		lowPass3.autoNormalize();

		for (uint i = 0; i < blockSize * m_numChannels; ++i)
			m_buffer[i] = 0.0f;

		lowPass1.resetHistory();
		lowPass2.resetHistory();
		lowPass3.resetHistory();

		PinkNoise pinkNoise;
		pinkNoise.setAlpha(1.3f);
		pinkNoise.init();

		const bool usePinkNoise = true;

		scalar maxVol = 0.0f;
		for (uint i = 0; i < blockSize; ++i)
		{
			scalar time = i / (scalar)m_samplesPerSec;

			uint rnd30;
			float noise;

			// Part 1
			if (usePinkNoise)
			{
				noise = pinkNoise.generateSample();
			}
			else
			{
				rnd30 = rand30() & 65535;
				noise = (rnd30 - 32768.0f) / 32768.0f;
			}

			// Exponential decay AR
			const scalar ARdecay1 = 15.0f;
			scalar AR = 1.0f * expf(-ARdecay1 * time);

			scalar part1;
			lowPass1.applyCont(&part1, AR * noise);
			//////////////////////////////////////////////////////////////////////////

			// Part 2
			if (usePinkNoise)
			{
				noise = pinkNoise.generateSample();
			}
			else
			{
				rnd30 = rand30() & 65535;
				noise = (rnd30 - 32768.0f) / 32768.0f;
			}

			// Exponential decay AR
			const scalar ARdecay2 = 5.0f;
			AR = 1.0f * expf(-ARdecay2 * time);

			scalar part2;
			lowPass2.applyCont(&part2, AR * noise);
			//////////////////////////////////////////////////////////////////////////

			// Part 3
			if (usePinkNoise)
			{
				noise = pinkNoise.generateSample();
			}
			else
			{
				rnd30 = rand30() & 65535;
				noise = (rnd30 - 32768.0f) / 32768.0f;
			}

			// Exponential decay AR
			const scalar ARdecay3 = 2.0f;
			AR = 1.0f * expf(-ARdecay2 * time);

			scalar part3;
			lowPass3.applyCont(&part3, AR * noise);
			//////////////////////////////////////////////////////////////////////////

			scalar snd = 0.6f * (0.4f * part1 + 0.55f * part2 + 1.0f * part3);

			if (!usePinkNoise)
			{
				// Amplify white noise sample
				snd *= 10.0f;
			}

			if (snd > maxVol)
				maxVol = snd;

			m_buffer[i*2 + 0] = snd;
			m_buffer[i*2 + 1] = snd;
		}
		if (applyReverb)
			testFDNReverb.processStereoInterleaved(m_buffer, blockSize);

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

