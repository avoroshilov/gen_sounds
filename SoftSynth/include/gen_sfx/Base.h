#pragma once

#include "helpers/Common.h"
#include "math/AuxMath.h"
#include "windows/Timer.h"

#include "windows/WaveOut_Simple.h"

class SoundPlay
{
protected:

	uint m_maxWaveDev = 0;

	uint m_numWaveDevs = 0;
	windows::WaveOutSimple * m_waveOut = nullptr;

	uint m_numChannels = 0;
	uint m_samplesPerSec = 0;

	scalar m_numSeconds = 0.0f;

	windows::Timer m_perfTimer;
	scalar m_execTime = 0.0f;

public:

	SoundPlay(uint numWaveDevsRequested, uint samplesPerSec, uint numChannels):
		m_samplesPerSec(samplesPerSec),
		m_numChannels(numChannels),
		m_numWaveDevs(numWaveDevsRequested)
	{
		m_waveOut = new windows::WaveOutSimple[m_numWaveDevs];
		for (uint i = 0; i < m_numWaveDevs; ++i)
			m_waveOut[i].startPlay();

		m_execTime = 0.0f;
	}
	~SoundPlay()
	{
		deinit();
	}

	void init()
	{
		m_waveOut = new windows::WaveOutSimple[m_numWaveDevs];
		for (uint i = 0; i < m_numWaveDevs; ++i)
			m_waveOut[i].startPlay();

		m_execTime = 0.0f;
	}
	void deinit()
	{
		for (uint i = 0; i < m_numWaveDevs; ++i)
			m_waveOut[i].stopPlay();
		delete [] m_waveOut;
	}

	scalar getTimings() const { return m_execTime; }
	scalar getNumSeconds() const { return m_numSeconds; }
	uint getNumChannels() const { return m_numChannels; }
	uint getSamplesPerSec() const { return m_samplesPerSec; }

	// Statistics, shows how much devices were actually used during playback
	uint getMaxWaveDevUsed() const { return m_maxWaveDev; }

	virtual void generate() = 0;
	virtual short * getBufferGenerated() = 0;
	virtual uint getBufferSizeGenerated() = 0;
	virtual void playSound() = 0;
	virtual void storeSound(const char * filename) = 0;

	void playSoundBuffer(short * buffer, uint totalNumSamplesPerChannel)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		// Output converted buffer
		bool devFound = false;
		for (uint i = 0; i < m_numWaveDevs; ++i)
		{
			if (m_waveOut[i].getBufferNum() == 0)
			{
				if (i > m_maxWaveDev)
					m_maxWaveDev = i;

				m_waveOut[i].play((char *)buffer, totalNumSamplesPerChannel * m_numChannels * sizeof(short));
				break;
			}
		}
		if (!devFound)
		{
			// If free device wasn't found, dinf device with lower amount of buffers scheduled and replace it
			
			// TODO: This won't work with simple sample players, since sample is 1 buffer;
			//	in that case, search over maximum time played will be required, or minimum time remained
			uint replacementDevice = 0;
			uint minBuffers = m_waveOut[replacementDevice].getBufferNum();
			for (uint i = 0; i < m_numWaveDevs; ++i)
			{
				uint curBuffersNum = m_waveOut[i].getBufferNum();
				if (curBuffersNum < minBuffers)
				{
					minBuffers = curBuffersNum;
					replacementDevice = i;
				}
			}

			m_waveOut[replacementDevice].reset();
			m_waveOut[replacementDevice].play((char *)buffer, totalNumSamplesPerChannel * m_numChannels * sizeof(short));
		}

		m_execTime = (scalar)m_perfTimer.time();

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
	}
};