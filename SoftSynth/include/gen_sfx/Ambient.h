#pragma once

#include <xmmintrin.h>

#include "windows/timer.h"
#include "formats/wav.h"

#include "windows/WaveOut_Simple.h"
#include "sfx/Oscillators.h"
#include "sfx/Phaser.h"
#include "sfx/Chorus.h"
#include "sfx/Filters.h"
#include "sfx/Reverb.h"
#include "sfx/ADSR.h"

#include "Base.h"

class SoundAmbient : public SoundPlay
{
protected:

	static const int c_numDevices = 2;
	static const int c_samplesPerSec = 44100;
	static const int c_numChannels = 2;

	size_t m_bufSize;
	float * m_buffer, * m_auxBuf0, * m_auxBuf1, * m_auxBuf2;
	short * m_bufInt;

	windows::Timer m_auxTimer;

	bool m_isGenerated = false;

	sfx::ADSRModulator m_ADSR;
	scalar applyADSR(scalar time)
	{
		return m_ADSR.getValue(time);
	}

	inline scalar calcFreqSemitone(scalar baseSemitone, scalar detune)
	{
		return powf(2.0f, (baseSemitone + detune - 49) / 12.0f) * 440.0f;
	}

	//inline scalar m_log2( scalar n )  
	//{
	//	// log(n)/log(2) is m_log2.  
	//	return log( n ) / log( (scalar)2 );  
	//}
	inline scalar calcBaseSemitone(scalar baseFreq)
	{
		// First, calculate the base semitone (when the string is open):
		//		baseN = 12 * m_log2( F(n) / 440.0 ) + 49
		scalar baseSemitone = 12 * m_log2(baseFreq / 440.0f) + 49;
		return baseSemitone;
	}

	scalar Osc3(scalar phase, scalar baseSemitone, scalar detune1, scalar detune2, scalar detune3)
	{
		scalar snd = 0.0f;
		snd += 0.33333f * sfx::oscTriangle(calcFreqSemitone(baseSemitone, -0.2f) * phase, 1.0f);
		snd += 0.33333f * sfx::oscTriangle(calcFreqSemitone(baseSemitone, 0.17f) * phase, 1.0f);
		snd += 0.33333f * sfx::oscTriangle(calcFreqSemitone(baseSemitone, -12.f) * phase, 1.0f);
		return snd;
	}

	bool isBasicTriadMajor(int baseSemitone)
	{
		// C = 4, D = 6, E = 8, F = 9, G = 11, A = 1, B = 3
		int noteSemitone = baseSemitone % 12;

		// Major: C, F, G 
		if (noteSemitone == 4 || noteSemitone == 9 || noteSemitone == 11)
			return true;

		// Minor: D, E, A, B [ and else? ]
		return false;
	}
	
public:

	int genType;

	SoundAmbient(scalar numSeconds = 19.0f):
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
	}

	~SoundAmbient()
	{
		delete [] m_buffer;
		delete [] m_auxBuf0;
		delete [] m_auxBuf1;
		delete [] m_auxBuf2;
		delete [] m_bufInt;
	}

	void generate()
	{
		using namespace sfx;

		const size_t toneBlockSize = (size_t)(m_numSeconds * m_samplesPerSec);

#if 1
		m_ADSR.setStart(0.0f, 0.0f);
		m_ADSR.setAttack(3.0f, 1.0f);
		m_ADSR.setDecay(0.0f, 1.0f);
		m_ADSR.setSustain(0.0f, 1.0f);
		m_ADSR.setRelease(6.0f, 0.0f);
		m_ADSR.setFinalAmp(0.0f, 0.0f);
		m_ADSR.sumDurations();
#else
		m_ADSR.SetStart(0.0f, 0.0f);
		m_ADSR.SetAttack(1.0f, 0.33f);
		m_ADSR.SetDecay(0.0f, 0.33f);
		m_ADSR.SetSustain(0.0f, 0.33f);
		m_ADSR.SetRelease(6.0f, 0.0f);
		m_ADSR.SetFinalAmp(0.0f, 0.0f);
		m_ADSR.SumDurations();
#endif

		const scalar adsrSustainTime = m_ADSR.getSustainTime();
		const scalar adsrReleaseTime = m_ADSR.getReleaseTime();

		LowPassIIR lowPassOsc3, generalLowPass;

		lowPassOsc3.setParams(6000.0f / m_samplesPerSec);
		lowPassOsc3.autoNormalize();
		lowPassOsc3.resetHistory();

		generalLowPass.setParams(5000.0f / m_samplesPerSec);
		generalLowPass.autoNormalize();
		generalLowPass.resetHistory();

		Phaser bassPhaser;
		bassPhaser.setAllpassNum(8);
		bassPhaser.setFeedback(0.075f);
		bassPhaser.setPhase(0.5f);
		bassPhaser.setRate(0.15f);
		//bassPhaser.SetRate(2.0f);
		bassPhaser.setDepth(0.7f);
		bassPhaser.setRange(0.0f, 100.0f);

		bassPhaser.init();

	
		const uint numChorusVoices = 4;
		Chorus bassChorus;
		bassChorus.setNumVoices(numChorusVoices);
		bassChorus.setFeedback(0.075f);
		bassChorus.setDepth(0.9f);

		bassChorus.init();


		ReverbFDN4 generalReverb;
		generalReverb.init();
		{
			uint delayLengths[] = { 1439, 2311, 3347, 4973 };

			generalReverb.initDLs(delayLengths);
			generalReverb.setHalflives(4500.0f, 500.0f);
			generalReverb.setFeedbackMatrix(ReverbMatrixTypes::eStautner2);
			generalReverb.setInOutRotations(90.0f, 180.0f);
			generalReverb.setModulation(60000.0f, 0.8f);
			generalReverb.setMixCoeffs(0.1f, 0.9f);
		}
		generalReverb.clearDLs();

		for (uint i = 0; i < numChorusVoices; ++i)
		{
			bassChorus.setPhase(i, i / (scalar)numChorusVoices);
			//bassChorus.SetRate(i, 0.15f);

			if (i == 100)
			{
				// Special case, low-freq shifter
				bassChorus.setRate(i, 0.1f);

				scalar rangeMin = (rand()%1000) / 1000.0f * 0.0009f + 0.0001f;
				scalar rangeMax = (rand()%1000) / 1000.0f * 0.0250f + 0.0250f;
				bassChorus.setRange(i, rangeMin * m_samplesPerSec, rangeMax * m_samplesPerSec);
			}
			else
			{
				bassChorus.setRate(i, 1.0f * ((rand()%1000) / 1000.0f * 5.0f + 0.5f));
#if 1
				//scalar rangeMin = (rand()%1000) / 1000.0f * 0.005f + 0.005f;
				//scalar rangeMax = (rand()%1000) / 1000.0f * 0.025f + 0.025f;
				//scalar rangeMin = (rand()%1000) / 1000.0f * 0.0009f + 0.0001f;
				//scalar rangeMax = (rand()%1000) / 1000.0f * 0.0025f + 0.0025f;
				scalar rangeMin = (rand()%1000) / 1000.0f * 0.0001f + 0.0001f;
				scalar rangeMax = (rand()%1000) / 1000.0f * 0.0010f + 0.0010f;
#else
				scalar rangeMin = 5.0f * (i+1) / (scalar)numChorusVoices;
				scalar rangeMax = 5.0f * (i+1) / (scalar)numChorusVoices;
#endif
				bassChorus.setRange(i, rangeMin * m_samplesPerSec, rangeMax * m_samplesPerSec);
			}
		}
		bassChorus.initPhases();


		Delay bassDelayL, bassDelayR;
		bassDelayL.setParams(2*m_samplesPerSec, (uint)(0.5f * m_samplesPerSec));
		bassDelayR.setParams(2*m_samplesPerSec, (uint)(0.5f * m_samplesPerSec));


		scalar maxVolume = 1.0f;

		const uint numParts = 4;

		bool isFinished[numParts];
		for (uint j = 0; j < numParts; ++j)
			isFinished[j] = false;
		
		scalar startTimes[numParts] =
		{
			0.0f,
			5.0f,
			10.0f,
			15.0f,
		};
		scalar baseSemitones[numParts] =
		{
			28,	// C3
			30,	// D3
			32,	// E3
			30,	// D3
		};
		scalar keyPressTimes[numParts] =
		{
			1.0f,
			1.0f,
			1.0f,
			1.0f,
		};

		// Not used anywhere - only for debugging/info
		scalar baseFreqs[numParts];
		for (uint j = 0; j < numParts; ++j)
			baseFreqs[j] = calcFreqSemitone(baseSemitones[j], 0.0f);

		scalar rndOffsets[32];
		for (uint j = 0; j < 32; ++j)
		{
			rndOffsets[j] = (rand()%1000) / 1000.0f;
		}

		for (uint j = 0; j < toneBlockSize; ++j)
		{
			const scalar phase = j / (scalar)m_samplesPerSec;

			// Freq = 1/20 Hz
			scalar globalCutoffLFO = sinf(0.05f * _2PI * phase/* + 3*PI/2.0f*/);

			scalar snd = 0.0f;
			
			// Mario coin effect
			//scalar freq = 420.0f;
			//if (phase < 0.1f)
			//	snd = OscSquare(freq * phase, 0.5f);
			//else if (phase < 0.5f)
			//{
			//	scalar relPhase = phase - 0.1f;
			//	snd = (1.0f - relPhase / 0.4f) * OscSquare(2.0f * freq * relPhase, 0.5f);
			//}
			////snd = OscSine(_2PI * freq * phase, 0.5f);

#if 1
			for (uint k = 0; k < numParts; ++k)
			{
				if (!isFinished[k] && phase > startTimes[k])
				{
					scalar baseSemitone = baseSemitones[k]; // 40
					bool isThisTriadMajor = isBasicTriadMajor((int)baseSemitone);

					scalar phase0 = phase - startTimes[k];

					scalar phaseOff0 = rndOffsets[(k*3+0)&31] * 0.2f;
					scalar snd0 = Osc3(phase0 + phaseOff0, baseSemitone, -0.2f, 0.17f, -12.0f);

					// CHORDS
#if 1
					scalar phaseOff1 = rndOffsets[(k*3+1)&31] * 0.2f;
					scalar phaseOff2 = rndOffsets[(k*3+2)&31] * 0.2f;

					if (isThisTriadMajor)
					{
						// Major chord:
						snd0 += Osc3(phase0 + phaseOff1, baseSemitone + 4, -0.2f, 0.17f, -12.0f);
						snd0 += Osc3(phase0 + phaseOff2, baseSemitone + 7, -0.2f, 0.17f, -12.0f);
					}
					else
					{
						// Minor chord:
						snd0 += Osc3(phase0 + phaseOff1, baseSemitone + 3, -0.2f, 0.17f, -12.0f);
						snd0 += Osc3(phase0 + phaseOff2, baseSemitone + 7, -0.2f, 0.17f, -12.0f);
					}

					snd0 /= 3.0f;
#endif

					//bool keyPressed = phase0 > keyPressTimes[k];
					//scalar adsrPhase = phase0;
					//if (keyPressed && adsrPhase > adsrReleaseTime)
					//{

					//}

					scalar adsrAmp = applyADSR(phase0);
					snd0 *= adsrAmp;

					if (phase0 > 1.0f && adsrAmp < 0.05f)
					{
						isFinished[k] = true;
					}

					snd += snd0;
				}
			}

			lowPassOsc3.setParams((3500.0f + 2500.0f * globalCutoffLFO) / m_samplesPerSec);
//			lowPassOsc3.setParams(6000.0f / m_SamplesPerSec);
			lowPassOsc3.autoNormalize();

			lowPassOsc3.applyCont(&snd, snd);
#else
			scalar freq = 220.0f;
			scalar modPhase = phase - (int)(phase / _PI2) * _PI2;
			const scalar dur = 0.1f;
			const scalar ampGrowth = 0.005f;
			if (modPhase < dur)
			{
				scalar amp = 1.0f;
				if (modPhase < ampGrowth) amp = modPhase / ampGrowth;
				else if (modPhase > dur - ampGrowth) amp = (dur - modPhase) / ampGrowth;
				snd = amp * OscTriangle(freq * phase, 0.5f);
			}
#endif

			snd = 0.1f * snd + 0.9f * bassChorus.update(snd);
			snd = 0.1f * snd + 0.9f * bassPhaser.update(snd);

			generalLowPass.applyCont(&snd, snd);

			scalar sndR, sndL;

			scalar delayDryWet = 0.1f;
			scalar delayGain = 0.5f;

#if 0
			// Normal
			scalar delayedSampleR = delayDryWet * snd + (1.0f - delayDryWet) * bassDelayR.GetSample();
			scalar delayedSampleL = delayDryWet * snd + (1.0f - delayDryWet) * bassDelayL.GetSample();
#else
			// Ping-pong
			scalar delayedSampleR = delayDryWet * snd + (1.0f - delayDryWet) * bassDelayL.getSample();
			scalar delayedSampleL = delayDryWet * snd + (1.0f - delayDryWet) * bassDelayR.getSample();
#endif

			bassDelayR.setSample(delayGain * (0*snd + delayedSampleR));
			bassDelayL.setSample(delayGain * (1*snd + delayedSampleL));

			sndR = snd + delayedSampleR;
			sndL = snd + delayedSampleL;

			maxVolume = m_max(maxVolume, fast_abs(snd));

			generalReverb.processCont(sndL, sndR, &sndL, &sndR);

			m_buffer[j*2+0] = sndR;
			m_buffer[j*2+1] = sndL;
		}

		// Convert float buf into shorts
		const scalar volume = 0.9f / maxVolume;
		for (uint j = 0; j < toneBlockSize * m_numChannels; ++j)
		{
			m_bufInt[j] = int(volume * clamp( m_buffer[j], -1.0f, 1.0f ) * 32767);
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

