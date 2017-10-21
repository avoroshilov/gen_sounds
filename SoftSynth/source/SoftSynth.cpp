#pragma warning( disable : 4996 )

#include <stdio.h>
#include <tchar.h>
#include <conio.h>
#include <Windows.h>
#include <MMSystem.h>
#include <math.h>
#include <float.h>

#include "helpers/Common.h"
#include "formats/tga.h"
#include "windows/Timer.h"
#include "windows/WaveOut_SFX.h"
#include "math/AuxMath.h"
#include "math/Vec2.h"
#include "math/Vec3.h"
#include "textures/GenTextures.h"

#include "sfx/Phaser.h"
#include "sfx/Oscillators.h"
#include "sfx/PlotFilter.h"

#include "gen_sfx/Formant.h"
#include "gen_sfx/Wind.h"
#include "gen_sfx/Ambient.h"
#include "gen_sfx/KarplusStrong.h"
#include "gen_sfx/ColorNoise.h"
#include "gen_sfx/GunShot.h"

class SineTable
{
public:

#define SINE_TABLE_SIZE			16384
#define SINE_TABLE_BITMASK		16383
#define SINE_TABLE_SIZE_2PI		2607.5945876176131812373915790953f

	static const int size = SINE_TABLE_SIZE;
	static float tbl[size];

	inline static void Initialize()
	{
		static bool initialized = false;

		if (!initialized)
		{
			float mul = _2PI / size;
			for (int i = 0; i < size; ++i)
			{
				tbl[i] = sinf(i * mul);
			}
			initialized = true;
		}
	}

	static float GetVal(float x)
	{
		Initialize();

#if 0
		if (x < 0.0f)
		{
			x -= _2PI * (static_cast<int>(x / _2PI) + 1);
		}
		else
		if (x > _2PI)
		{
			x -= _2PI * static_cast<int>(x / _2PI);
		}
		return tbl[static_cast<int>(x / _2PI * size)];
#else
		return tbl[static_cast<int>(x * SINE_TABLE_SIZE_2PI) & SINE_TABLE_BITMASK];
#endif
	}
};

float SineTable::tbl[SineTable::size];

#include "windows/WaveOut_Simple.h"
#define FILTERS_AUTONORMALIZE_SINETABLE
#include "sfx/Filters.h"
#include "sfx/Delay.h"
#include "sfx/ADSR.h"

#define PHYSICAL_STRINGS 0

#if (PHYSICAL_STRINGS == 1)
#include "PhysicalStrings.h"
#endif

int keyFromChar(char inChar)
{
	int order = -1;

	switch (inChar)
	{
	case '1': order =  0; break;
	case '2': order =  1; break;
	case '3': order =  2; break;
	case '4': order =  3; break;
	case '5': order =  4; break;
	case '6': order =  5; break;
	case '7': order =  6; break;
	case '8': order =  7; break;
	case '9': order =  8; break;
	case '0': order =  9; break;
	case '-': order = 10; break;
	case '=': order = 11; break;
	case '\\':order = 12; break;

	case 'q': order = 13; break;
	case 'w': order = 14; break;
	case 'e': order = 15; break;
	case 'r': order = 16; break;
	case 't': order = 17; break;
	case 'y': order = 18; break;
	case 'u': order = 19; break;
	case 'i': order = 20; break;
	case 'o': order = 21; break;
	case 'p': order = 22; break;
	case '[': order = 23; break;
	case ']': order = 24; break;

	case 'a': order = 25; break;
	case 's': order = 26; break;
	case 'd': order = 27; break;
	case 'f': order = 28; break;
	case 'g': order = 29; break;
	case 'h': order = 30; break;
	case 'j': order = 31; break;
	case 'k': order = 32; break;
	case 'l': order = 33; break;
	case ';': order = 34; break;
	case '\'':order = 35; break;

	case 'z': order = 36; break;
	case 'x': order = 37; break;
	case 'c': order = 38; break;
	case 'v': order = 39; break;
	case 'b': order = 40; break;
	case 'n': order = 41; break;
	case 'm': order = 42; break;
	case ',': order = 43; break;
	case '.': order = 44; break;
	case '/': order = 45; break;
	}

	return order;
}

int freqFromChar(char key)
{
	int numKeys = 46;

	int order = keyFromChar(key);
	if (order == -1)
		return -1;

#if 0
	int begFreq = 30;
	int endFreq = 4000;

	return (order * (endFreq - begFreq)) / (float)numKeys + begFreq;
#else
	int keyNumber = order + 54; /* MIDI-49-keyboard */ // 27; /* real piano */
	float freq = 8.176005152f * expf(0.5776201269e-1f * keyNumber) - 0.5479535474e-3f;
	
	return (int)freq;
#endif
}

//////////////////////////////////////////////////////////////////////////
void LowPassWrap(float * output, float * samples, uint length, float freqNormalized)
{
	//LowPassIIR::getInstance().setParams(freqNormalized);
	//LowPassIIR::getInstance().autoNormalizeS();

	sfx::LowPassIIR::getInstance().resetHistory();
	for (uint i = 0; i < length; ++i)
	{
		sfx::LowPassIIR::getInstance().applyCont(output + i, samples[i]);
	}
}
void LowPass1pWrap(float * output, float * samples, uint length, float freqNormalized)
{
	sfx::LowPass1pIIR::getInstance().setParams(freqNormalized);
	sfx::LowPass1pIIR::getInstance().resetHistory();
	for (uint i = 0; i < length; ++i)
	{
		sfx::LowPass1pIIR::getInstance().applyCont(output + i, samples[i]);
	}
}
void HighPassWrap(float * output, float * samples, uint length, float freqNormalized)
{
	sfx::HighPassIIR::getInstance().setParams(freqNormalized);
	sfx::HighPassIIR::getInstance().autoNormalize(0.49f);

	sfx::HighPassIIR::getInstance().apply(output, samples, length);
}
void BandPassWrap(float * output, float * samples, uint length, float freqNormalized)
{
	sfx::BandPassIIR::getInstance().resetHistory();
	for (uint i = 0; i < length; ++i)
	{
		sfx::BandPassIIR::getInstance().applyCont(output + i, samples[i]);
	}
}
void NotchWrap(float * output, float * samples, uint length, float freqNormalized)
{
	sfx::NotchIIR::getInstance().resetHistory();
	for (uint i = 0; i < length; ++i)
	{
		sfx::NotchIIR::getInstance().applyCont(output + i, samples[i]);
	}
}
//////////////////////////////////////////////////////////////////////////

int _tmain(int argc, _TCHAR* argv[])
{
	using namespace windows;
	using namespace sfx;

	// Example on how to use new classes

	auto playStoreWait = [](SoundPlay * sound, const char * filename = nullptr)
	{
		sound->playSound();
		if (filename)
		{
			sound->storeSound(filename);
		}
		windows::Timer timer;
		timer.start();
		scalar timeToWait = sound->getNumSeconds();
		while (timer.time() < 1000 * timeToWait);
	};

	const bool enableSamples = true;
	const bool enableSequences = true;
	const bool enableInfPlayback = false;
	const bool enableInfSequence = true;

	const bool storeWAVFiles = true;
	if (enableSamples)
	{
		const float durationSec = 1.5f;
		printf("Formant synthesis (%.1fs)\n", durationSec);

		enum FormantVowel
		{
			eI,
			eAE,
			eA,
			eO,
			eU,
			eWeirdI
		};
		FormantVowel vowel = FormantVowel::eO;

		if (vowel == FormantVowel::eI)
		{
			// I (peep)
			FormantSettings::BandPass0Seq::Freq = 270.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 0.9f;
			FormantSettings::BandPass1Seq::Freq = 3190.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.5f;
			FormantSettings::BandPass2Seq::Freq = 3010.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.4f;
		}
		else if (vowel == FormantVowel::eAE)
		{
			// ae (and)
			FormantSettings::BandPass0Seq::Freq = 660.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 0.9f;
			FormantSettings::BandPass1Seq::Freq = 1720.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.25f;
			FormantSettings::BandPass2Seq::Freq = 2410.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.8f;	// 0.08
		}
		else if (vowel == FormantVowel::eA)
		{
			// /\ (up)
			FormantSettings::BandPass0Seq::Freq = 730.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 0.9f;
			FormantSettings::BandPass1Seq::Freq = 1090.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.57f;
			FormantSettings::BandPass2Seq::Freq = 2440.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.04f;
		}
		else if (vowel == FormantVowel::eO)
		{
			// O (old)
			FormantSettings::BandPass0Seq::Freq = 570.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 1.0f;
			FormantSettings::BandPass1Seq::Freq = 840.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.32f;
			FormantSettings::BandPass2Seq::Freq = 2410.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.21f;

		}
		else if (vowel == FormantVowel::eU)
		{
			// U (ooze)
			FormantSettings::BandPass0Seq::Freq = 440.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 0.9f;
			FormantSettings::BandPass1Seq::Freq = 1020.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.26f;
			FormantSettings::BandPass2Seq::Freq = 2240.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.02f;
		}
		else if (vowel == FormantVowel::eWeirdI)
		{
			// weird i
			FormantSettings::BandPass0Seq::Freq = 190.0f;	FormantSettings::BandPass0Seq::Q = 0.98f;	FormantSettings::BandPass0Seq::Vol = 0.9f;
			FormantSettings::BandPass1Seq::Freq = 1330.0f;	FormantSettings::BandPass1Seq::Q = 0.96f;	FormantSettings::BandPass1Seq::Vol = 0.32f;
			FormantSettings::BandPass2Seq::Freq = 3000.0f;	FormantSettings::BandPass2Seq::Q = 0.99f;	FormantSettings::BandPass2Seq::Vol = 0.26f;
		}

		SoundFormant testSoundGenFormant(durationSec);
		playStoreWait(&testSoundGenFormant, storeWAVFiles ? "formant.wav" : nullptr);
	}
	if (enableSamples)
	{
		const float durationSec = 3.0f;
		printf("Phaser/flanger wind (%.1fs)\n", durationSec);
		SoundWind testSoundGenWind(durationSec);
		playStoreWait(&testSoundGenWind, storeWAVFiles ? "wind.wav" : nullptr);
	}
	if (enableSamples)
	{
		const float durationSec = 26.0f;
		printf("Ambient tune (%.1fs)\n", durationSec);
		SoundAmbient testSoundGenAmbient(durationSec);
		playStoreWait(&testSoundGenAmbient, storeWAVFiles ? "ambient.wav" : nullptr);
	}
	if (enableSamples)
	{
		const float durationSec = 3.0f;
		printf("Extended Karplus-Strong (%.1fs)\n", durationSec);
		SoundKarplusStrong testSoundGenKarplusStrong(durationSec);
		playStoreWait(&testSoundGenKarplusStrong, storeWAVFiles ? "karplus_strong.wav" : nullptr);
	}
	if (enableSamples)
	{
		const float durationSec = 3.0f;
		printf("Color noise (%.1fs)\n", durationSec);
		SoundColorNoise testSoundGenColorNoise(durationSec);

		// Noise generator type:
		// eAlpha - adjustable, default
		//	alpha in [0.0, 2.0) where 0.0 is white, and 1.9 is red
		//	(see ColorNoiseSettings::NoiseAlpha)
		//testSoundGenColorNoise.colorNoiseType = SoundColorNoise::eAlpha;

		// Non-adjustable
		//testSoundGenColorNoise.colorNoiseType = SoundColorNoise::eOctave;
		//testSoundGenColorNoise.colorNoiseType = SoundColorNoise::eCustom;

		playStoreWait(&testSoundGenColorNoise, storeWAVFiles ? "color_noise.wav" : nullptr);
	}
	if (enableSamples)
	{
		const float durationSec = 3.0f;
		printf("Gunshot (%.1fs)\n", durationSec);
		SoundGunShot testSoundGenGunShot(durationSec);
		playStoreWait(&testSoundGenGunShot, storeWAVFiles ? "gunshot.wav" : nullptr);
	}
	if (enableSamples)
	{
		using namespace math;
		using namespace sfx;

		{
			class SoundFXEnv
			{
			protected:

				Vec3 m_sourcePos;
				Vec3 m_sourceVel;

				scalar m_azimuth;
				scalar m_relVel;

				windows::WaveOutSfx * m_soundFX;

			public:

				SoundFXEnv():
					m_azimuth(0.0f),
					m_relVel(0.0f),
					m_soundFX(0)
				{
					m_sourcePos = Vec3C(0.0f, 0.0f, 0.0f);
					m_sourceVel = Vec3C(0.0f, 0.0f, 0.0f);
				}

				void SetSoundFX(windows::WaveOutSfx * pSoundFX) { m_soundFX = pSoundFX; }
				scalar GetAzimuth() { return m_azimuth; }
				scalar GetRelativeVel() { return m_relVel; }

				void UpdateSource(const Vec3 & srcPos, const Vec3 & srcVel)
				{
					m_sourcePos = srcPos;
					m_sourceVel = srcVel;
				}

				void UpdateReciever(const Vec3 & recvPos, const Vec3 & recvVel, const Vec3 & recvDir, const Vec3 & recvUp)
				{
					const Vec3 recvToSource = m_sourcePos - recvPos;
					const Vec3 recvToSourceDir = recvToSource.getNormalized();

					const Vec3 recvStrafe = recvDir.cross(recvUp);
					const Vec3 recvToSourceDirHorizontal = recvToSource - recvUp * recvToSource.dot(recvUp);

					// TODO: optimize this
					scalar sign = recvToSourceDirHorizontal.dot(recvStrafe) < 0.0f ? -1.0f : 1.0f;
					m_azimuth = sign * acosf( recvDir.dot(recvToSourceDirHorizontal.getNormalized()) );

					m_relVel = recvToSourceDir.dot(m_sourceVel - recvVel);		// Negative if moving towards

					const scalar SND_SPEED = 331.0f;							// 311.0f is the speed of the sound
					scalar dopplerEffect = 1.0f - m_relVel / SND_SPEED;

					scalar distanceSq = (m_sourcePos - recvPos).sqLen();

					scalar distFade = clamp( (distanceSq > 1e-4f) ? (1.0f / distanceSq) : 1.0f, 0.0f, 1.0f );

					m_soundFX->setDistanceFade(distFade);
					m_soundFX->setTimeScale(dopplerEffect);
					m_soundFX->setAzimuth(m_azimuth);
				}

			};

			const uint numChannels = 2;
			const uint samplesPerSec = 44100;

			int seconds = 5;
			uint blockSize = seconds * samplesPerSec;

			printf("Doppler + sound positioning (%.1fs)\n", (float)seconds);

			PinkNoise pinkNoise;
			pinkNoise.setAlpha(1.3f);
			pinkNoise.init();
			scalar *tstWave = new scalar[numChannels * blockSize];
			for (uint i = 0; i < blockSize; ++i)
			{
				scalar val = oscTriangle(200.0f, 1.0f, 0.0f, i / (scalar)samplesPerSec, 0.95f) + 0.4f*oscTriangle(100.0f, 1.0f, 0.0f, i / (scalar)samplesPerSec, 0.95f) + 0.1f*oscSquare(50.0f, 1.0f, 0.0f, i / (scalar)samplesPerSec, 0.5f) + 0.05f*pinkNoise.generateSample();
				tstWave[i*2 + 0] = val;
				tstWave[i*2 + 1] = val;
			}

			plotWaveformChannel("tstwf.tga", tstWave, blockSize / 20, 2048, 256);

			short *tstWaveShort = new short[numChannels * blockSize];

			scalar volume = 0.1f;
			for (uint j = 0; j < blockSize * numChannels; ++j)
			{
				tstWaveShort[j] = int(clamp( volume * tstWave[j], -1.0f, 1.0f ) * 32767);
			}

#if 1
			WaveOutSfx testWaveOutSfx;
			testWaveOutSfx.startPlay();

			// Hardcoded, precalculated

			testWaveOutSfx.setIntensity(2000.0f);
#if 0
			// 90 km/h
			testWaveOutSfx.setDistanceFade(0.0002f);
			testWaveOutSfx.setTimeScale(1.0754f);
			testWaveOutSfx.setAzimuth(-1.5f);
#else
			testWaveOutSfx.setDistanceFade(0.000123f);
			testWaveOutSfx.setTimeScale(1.1358762f);
			testWaveOutSfx.setAzimuth(-1.54f);
#endif

			testWaveOutSfx.playSoundStereo(tstWaveShort, blockSize);
			testWaveOutSfx.setSampleRecordingState(true);

			Timer testTimer;
			testTimer.start();

			scalar time = 0.0f;

			Vec3 sourcePos = Vec3C(-90.0f, 0.0f, -3.0f);
			Vec3 sourceVel = Vec3C( 45.0f, 0.0f,  0.0f);	// 90 km/h
			//Vec3 sourceVel = Vec3C( 15.0f, 0.0f,  0.0f);	// 90 km/h

			scalar prevDistance = FP_MAX;

			SoundFXEnv testSoundFX;
			testSoundFX.SetSoundFX(&testWaveOutSfx);

			Vec3 recieverPos = Vec3C(0.0f, 0.0f,  0.0f);
			Vec3 recieverVel = Vec3C(0.0f, 0.0f,  0.0f);
			Vec3 recieverDir = Vec3C(0.0f, 0.0f, -1.0f);
			Vec3 recieverUp  = Vec3C(0.0f, 1.0f,  0.0f);

			double prevTime = testTimer.time();
			while (testWaveOutSfx.getBufferNum() != 0)
			{

				double time = testTimer.time();

				double dt = (time - prevTime) / 1000.0f;
				sourcePos += sourceVel * (float)dt;

				scalar distanceSq = (sourcePos - recieverPos).sqLen();

#if 0
				scalar sign = sourcePos.x < 0.0f ? -1.0f : 1.0f;	// hack
				//scalar az = sign * PI * 0.5f * (1.0f - acosf( Vec3(0.0f, 0.0f, -1.0f).dot((sourcePos - recieverPos).getNormalized() )));
				scalar az = sign * acosf( Vec3(0.0f, 0.0f, -1.0f).dot((sourcePos - recieverPos).getNormalized()) );

				Vec3 diffSrcRcv = (sourcePos - recieverPos).getNormalized();
				scalar relativeVel = diffSrcRcv.dot(sourceVel);		// Negative if moving towards

				scalar dopplerEffect = 1.0f - relativeVel / 331.0f;	// 311.0f is the speed of sound

				scalar distFade = clamp( (distanceSq > 1e-4f) ? (1.0f / distanceSq) : 1.0f, 0.0f, 1.0f );
				testWaveOutSfx.SetDistanceFade(distFade);
				testWaveOutSfx.SetTimeScale(dopplerEffect);
				testWaveOutSfx.SetAzimuth(az);
#else
				testSoundFX.UpdateSource(sourcePos, sourceVel);
				testSoundFX.UpdateReciever(recieverPos, recieverVel, recieverDir, recieverUp);
#endif

				static bool printed = false;
				if (distanceSq > prevDistance && !printed)
				{
					printf("reached!\n");
					printed = true;
				}

				static double elTime = 0.0;
				if (elTime > 0.25)
				{
					printf("az:%f, spx:%f\n", testSoundFX.GetAzimuth(), sourcePos.x);
					elTime = 0.0;
				}
				elTime += dt;

				prevTime = time;

				prevDistance = distanceSq;
			}

			printf("total time: %f\n", testTimer.time());

			testWaveOutSfx.stopPlay();
#else
			WaveOut testWaveOut;
			testWaveOut.StartPlay();
			testWaveOut.Play((char *)tstWaveShort, numChannels * blockSize * sizeof(short));

			while (testWaveOut.GetBufferNum() != 0);

			testWaveOut.StopPlay();
#endif

			delete [] tstWaveShort;
			delete [] tstWave;

			formats::writeWAV("doppler.wav", testWaveOutSfx.m_samplesInterleaved.data(), blockSize, samplesPerSec, numChannels);
		}
	}
	if (enableSequences)
	{
		// WARNING: sequence, not a single sample, thus would require additional mixing to save to a file
		const float durationSec = 24.5f;	// one cycle of 1.0 speed is 14.5
		printf("Sequence: horizons (%.1fs)\n", durationSec);

		const float speedMul = 1.2f;
		const float stringVibrationSec = 3.0f;
		SoundKarplusStrong genKarplusStrong(stringVibrationSec);
		KSSynthSequencer testKSSynthSequencer(&genKarplusStrong);

		InitKS_Horizons(&testKSSynthSequencer);
		testKSSynthSequencer.setSampleRecordingState(true);
		testKSSynthSequencer.setSampleRecordingSpeed(speedMul);

		float oldVirbAmp = KarplusStrongSettings::VibrAmp;

		const bool goofySound = false;
		if (goofySound)
		{
			KarplusStrongSettings::VibrAmp = 0.02f;
		}

		Timer ksSeqTimer;
		double dt = 0.0, elapsedTime = 0.0;
		ksSeqTimer.start();
		while (true)
		{
			// Avoid too frequent updates to reduce FP issues
			while (ksSeqTimer.time() < 1);
			testKSSynthSequencer.update(dt * 0.001 * speedMul);
			dt = ksSeqTimer.time();
			ksSeqTimer.start();
			elapsedTime += dt;
			if (elapsedTime > durationSec * 1000)
				break;
		}

		KarplusStrongSettings::VibrAmp = oldVirbAmp;

		if (testKSSynthSequencer.getSampleRecordingState())
		{
			testKSSynthSequencer.prepareRecording();
			formats::writeWAV(
				"ks_horizons.wav",
				testKSSynthSequencer.getBufferRecorded(),
				testKSSynthSequencer.getBufferSizeRecorded()/genKarplusStrong.getNumChannels(),
				genKarplusStrong.getSamplesPerSec(),
				genKarplusStrong.getNumChannels()
				);
		}
	}
	if (enableSequences)
	{
		// WARNING: sequence, not a single sample, thus would require additional mixing to save to a file
		const float durationSec = 14.5f;	// one cycle is 14.5
		printf("Sequence: strum (%.1fs)\n", durationSec);

		const float speedMul = 1.2f;
		const float stringVibrationSec = 3.0f;
		SoundKarplusStrong genKarplusStrong(stringVibrationSec);
		KSSynthSequencer testKSSynthSequencer(&genKarplusStrong);

		InitKS_Strumming(&testKSSynthSequencer);

		Timer ksSeqTimer;
		double dt = 0.0, elapsedTime = 0.0;
		ksSeqTimer.start();
		while (true)
		{
			testKSSynthSequencer.update(dt * 0.001 * speedMul);
			dt = ksSeqTimer.time();
			ksSeqTimer.start();
			elapsedTime += dt;
			if (elapsedTime > durationSec * 1000)
				break;
		}
	}
	if (enableInfPlayback)
	{
		printf("Sequence: waves (inf)\n");

		const uint c_numChannels = 2;
		const uint c_samplesPerSec = 44100;

		uint invWaveBlockSize;
		int infWaveSeconds = 2;

		invWaveBlockSize = infWaveSeconds * c_samplesPerSec;

		float * infWave = new float [invWaveBlockSize * c_numChannels];
		short * infWaveShort = new short [invWaveBlockSize * c_numChannels];

		float infWaveVolume = 0.5f;

		const uint numWaveDevs = 50;
		uint maxWaveDev = 0;
		WaveOutSimple testWaveOut[numWaveDevs];

		for (uint i = 0; i < numWaveDevs; ++i)
			testWaveOut[i].startPlay();

		const uint numTones = 20;
		scalar toneOffset[numTones];
		scalar toneVelocity[numTones];

		const float c0Freq = 16.35f;
		const float c4Freq = 261.63f;
		for (uint i = 0; i < numTones; ++i)
		{
			toneOffset[i] = 0.0f;
			toneVelocity[i] = c0Freq * (i+1);
		}

		int blockNum = 0;
		float time = 0.0f;
		while (1)
		{
			using namespace math;

			auto testTone = [](float t) -> Vec2
			{
				Vec2 out = Vec2C(0.0f, 0.0f);

				{
					float toneAmp = 0.3f * (sinf(t * 1.5f) * 0.4f + 0.6f);
					float phase = _2PI * 440.0f * t;
					float sine = toneAmp * sinf(phase);
					out.x += sine;
					out.y += sine;
				}
				{
					float toneAmp = 0.3f * (cosf(t * 0.9f + 10.0f) * 0.4f + 0.6f);
					float phase = _2PI * 220.0f * t;
					float sine = toneAmp * sinf(phase);
					out.x += sine;
					out.y += sine;
				}

				return out;
			};

			for (uint i = 0; i < invWaveBlockSize; ++i)
			{
				Vec2 infWaveSound = testTone(time + i / (float)c_samplesPerSec);
				infWave[i*2 + 0] = infWaveSound.x;
				infWave[i*2 + 1] = infWaveSound.y;
			}
			time += infWaveSeconds;

			for (uint j = 0; j < invWaveBlockSize * c_numChannels; ++j)
			{
				infWaveShort[j] = int(infWaveVolume * clamp( infWave[j], -1.0f, 1.0f ) * 32767);
			}

			// Spinwait until sane number of buffers queued
			do
			{
			} while (testWaveOut[0].getBufferNum() > 5);

			printf("Enqueue block #%d\n", blockNum);
			testWaveOut[0].play((char *)infWaveShort, invWaveBlockSize * c_numChannels * sizeof(short));
			++blockNum;
		}

		for (uint i = 0; i < numWaveDevs; ++i)
			testWaveOut[i].stopPlay();

		delete [] infWaveShort;
		delete [] infWave;
	}
	if (enableInfSequence)
	{
		printf("Sequence: chord triage (inf)\n");

		const uint c_numChannels = 2;
		const uint c_samplesPerSec = 44100;

		uint blockSize;
		int seconds = 2;

		scalar volume = 1.0f;

		blockSize = seconds * c_samplesPerSec;

		float * wave;
		short * waveShort;

		wave = new float [blockSize * c_numChannels];
		waveShort = new short [blockSize * c_numChannels];

		for (uint i = 0; i < blockSize * c_numChannels; ++i)
			wave[i] = 0.0f;

		const uint numWaveDevs = 15;
		uint maxWaveDev = 0;
		windows::WaveOutSimple testWaveOut[numWaveDevs];

		for (uint i = 0; i < numWaveDevs; ++i)
			testWaveOut[i].startPlay();

		// Prepare reverb
		ReverbFDN4 testFDNReverb;
		testFDNReverb.init();
		{
			uint delayLengths[] = { 601, 919, 1093, 1187 };
			testFDNReverb.initDLs(delayLengths);
			testFDNReverb.setHalflives(2000.0f, 900.0f);
			testFDNReverb.setFeedbackMatrix(ReverbMatrixTypes::eHouseholder6);
			testFDNReverb.setInOutRotations(0.0f, 0.0f);
			testFDNReverb.setModulation(70000.0f, 0.8f);
			testFDNReverb.setMixCoeffs(0.4f, 0.6f);
		}

		// SNARE
		HighPassIIR & filterPass = HighPassIIR::getInstance();

		BandPassIIR bandPass;
		bandPass.setParams(0.25f, 0.5f);
		bandPass.autoNormalize(0.25f);

		filterPass.setParams((500.0f + (rand()%100)) / c_samplesPerSec);
		filterPass.autoNormalize(0.5f);

		filterPass.resetHistory();
		bandPass.resetHistory();

		for (uint i = 0; i < blockSize; ++i)
		{
			uint rnd30 = rand30() & 65535;
			float noise = (rnd30 - 32768.0f) / 32768.0f;

			scalar time = i / (scalar)c_samplesPerSec;

#if 1
			// Exponential decay AR
			const scalar ARdecay = 10.0f;
			scalar AR = 1.0f * expf(-ARdecay * time);
#else
			// Linear decay AR
			const scalar ARdecay = 2.0f;
			scalar AR = clamp(1.0f - ARdecay * time, 0.0f, 1.0f);
#endif

			scalar intermSnd;
			filterPass.applyCont(&intermSnd, AR * noise);

			scalar snd;
			bandPass.applyCont(&snd, intermSnd);

			// Adding some sine bass
			snd += AR * 1.0f * sinf(time * _2PI * 60.0f);

			const scalar eps = 1e-5f;
			if (AR < eps)
			{
				break;
			}

			wave[i*2 + 0] = snd;
			wave[i*2 + 1] = snd;
		}

		for (uint j = 0; j < blockSize * c_numChannels; ++j)
		{
			waveShort[j] = int(volume * clamp( wave[j], -1.0f, 1.0f ) * 32767);
		}

		double elapsedTime = 0.0;
		windows::Timer dtTimer;
		dtTimer.start();

		windows::Timer DBGoverallTimer;
		DBGoverallTimer.start();
		double DBGoverallTime = 0.0;
		double DBGoveralTimeError = 0.0;

		const double metronomeInterval = 1000.0f;
		uint metronomeSounded = 0;
		double metronomeTime = 0.0;

		// Tone setup
		const scalar toneSeconds = 2.0f;
		const uint toneBlockSize = (uint)(toneSeconds * c_samplesPerSec);
		float * toneWave = new float [toneBlockSize * c_numChannels];
		short * toneWaveShort = new short [toneBlockSize * c_numChannels];

		bool toneValid = false;
		const double toneInterval = 250.0f;
		uint toneSounded = 0;
		double toneTime = 0.0;

		double toneFreqTime = 0.0;

		while (true)
		{
			// Regenerate tone if needed
			if (!toneValid)
			{
				testFDNReverb.clearDLs();

#if 0
				// Raw frequencies
				const scalar minFreq = 200;
				const scalar maxFreq = 700;
				const scalar toneChangeSpeed = 0.0005f;
				scalar freq = 0.5f * (sinf(toneChangeSpeed * toneFreqTime) + 1.0f) * (maxFreq - minFreq) + minFreq;
#else
				// Notes
				const scalar minNote = 20;
				const scalar maxNote = 56;
				const scalar period = 2000.0f;
#if 0
				// SINE WAVE
				scalar phase = toneFreqTime / period * _2PI - PI * 0.5f;
				scalar note = (int)(0.5f * (sinf(phase) + 1.0f) * (maxNote - minNote) + minNote);
#else
				// TRIANGLE WAVE
				scalar phase = (scalar)toneFreqTime;
				scalar clampedPhase = (phase - (int)(phase / period) * period) / period;
				scalar trif;
				if (clampedPhase > 0.5f)
				{
					trif = 1.0f - 4.0f * (clampedPhase - 0.5f);
				}
				else
				{
					trif = 4.0f * clampedPhase - 1.0f;
				}

				const scalar rangeNoteShift = 5;
				scalar subPhase = (scalar)toneFreqTime / 1700.0f * _2PI - PI * 0.5f;
				scalar noteShift = sinf(subPhase) * rangeNoteShift;

				// In RELEASE build FP math is slightly less precise, hence we need to add 1e-4 eps
				const scalar RELEASEeps = 1e-4f;
				scalar note = (scalar)(0.5f * (trif + 1.0f) * (maxNote - minNote) + minNote + 0.0f * noteShift + RELEASEeps);
#endif

#if 1
				scalar freq = powf(2.0f, (note - 49) / 12.0f) * 440.0f; // A = 440Hz
#else
				scalar freq = powf(2.0f, (note - 49) / 12.0f) * 432.0f;	// A = 432Hz
#endif

#define TRIAD_NONE			0
#define TRIAD_MAJOR			1
#define TRIAD_MINOR			2
#define TRIAD_DIMINISHED	3
#define TRIAD_AUGMENTED		4

#define OSC_SINE			1
#define OSC_TRIANGLE		2
#define OSC_SQUARE			3

#define TRIAD_TYPE			TRIAD_DIMINISHED
#define OSC_TYPE			OSC_SINE			

#if (OSC_TYPE == OSC_SINE)
#	define OSC_FUNC			oscSine
#endif
#if (OSC_TYPE == OSC_TRIANGLE)
#	define OSC_FUNC			oscTriangle
#endif
#if (OSC_TYPE == OSC_SQUARE)
#	define OSC_FUNC			oscSquare
#endif

#if (TRIAD_TYPE == TRIAD_MAJOR)
				scalar note_third = note + 4;
				scalar note_fifth = note_third + 3;
#endif
#if (TRIAD_TYPE == TRIAD_MINOR)
				scalar note_third = note + 3;
				scalar note_fifth = note_third + 4;
#endif
#if (TRIAD_TYPE == TRIAD_DIMINISHED)
				scalar note_third = note + 3;
				scalar note_fifth = note_third + 3;
#endif
#if (TRIAD_TYPE == TRIAD_AUGMENTED)
				scalar note_third = note + 4;
				scalar note_fifth = note_third + 4;
#endif

#if (TRIAD_TYPE == TRIAD_NONE)
#else
				scalar freqThird = powf(2.0f, (note_third - 49) / 12.0f) * 440.0f;
				scalar phaseThird = (rand()%1000) / 1000.0f;
				scalar freqFifth = powf(2.0f, (note_fifth - 49) / 12.0f) * 440.0f;
				scalar phaseFifth = (rand()%1000) / 1000.0f;
#endif

#endif

				float dutyTimeMin = 0.2f;
				float dutyTimeMax = 0.6f;
				float dutyTimePeriod = 6000.0f;
				float dutyTime = 0.5f * (sinf((float)toneFreqTime / dutyTimePeriod * _2PI) + 1.0f) * (dutyTimeMax - dutyTimeMin) + dutyTimeMin;
				for (uint i = 0; i < toneBlockSize; ++i)
				{
					scalar time = i / (scalar)c_samplesPerSec;

#define	CALC_ENVELOPE(envTime) \
				if (envTime > 0.0f) \
				{ \
					if (envTime > ARattackTime) \
					{ \
						AR = 1.0f * expf(-ARdecay * (envTime - ARattackTime)); \
					} \
					else \
					{ \
						AR = envTime / ARattackTime; \
					} \
				} \
				else \
				{ \
					AR = 0.0f; \
				}

					// Exponential decay AR
#if 1
					const scalar ARdecay = 4.0f, ARattackTime = 0.005f;
#else
					// How it was tuned long ago
					const scalar ARdecay = 6.0f, ARattackTime = 0.0f;
#endif
					scalar AR;

					scalar snd;
#if (TRIAD_TYPE == TRIAD_NONE)
					CALC_ENVELOPE(time);
					scalar val = OSC_FUNC(freq, 1.0f, 0.0f, time, dutyTime);
					snd = AR * val;
#else
					scalar valR = OSC_FUNC(freq, 1.0f, 0.0f, time, dutyTime);
					scalar valT = OSC_FUNC(freqThird, 1.0f, _2PI * phaseThird, time, dutyTime);
					scalar valF = OSC_FUNC(freqFifth, 1.0f, _2PI * phaseFifth, time, dutyTime);

					const scalar noteTimeVar = 0.025f;

					snd = 0.0f;

					scalar noteTime;

					noteTime = time;
					CALC_ENVELOPE(noteTime);
					snd += AR * valR;

					noteTime = time + noteTimeVar * (phaseThird - 1.0f);
					CALC_ENVELOPE(noteTime);
					snd += AR * valT;

					noteTime = noteTime + noteTimeVar * (phaseFifth - 1.0f);
					CALC_ENVELOPE(noteTime);
					snd += AR * valF;

					snd *= 0.334f;
#endif

#undef CALC_ENVELOPE

					toneWave[i*2 + 0] = snd;
					toneWave[i*2 + 1] = snd;
				}

#if (APPLY_TONE_REVERB == 1)
				testFDNReverb.processStereoInterleaved(toneWave, toneBlockSize);
#endif

				for (uint i = 0; i < toneBlockSize; ++i)
				{
					short shortValL = int(clamp( toneWave[i*2 + 0], -1.0f, 1.0f ) * 32767);
					short shortValR = int(clamp( toneWave[i*2 + 1], -1.0f, 1.0f ) * 32767);
					toneWaveShort[i*2 + 0] = shortValL;
					toneWaveShort[i*2 + 1] = shortValR;
				}

				toneValid = true;
			}

			DBGoverallTime = DBGoverallTimer.time();

#if 1
			__int64 dt_cnt = dtTimer.deltaCount();
			dtTimer.start();

			double dt = dtTimer.timeFromDeltaCount(dt_cnt);
#else
			double dt = dtTimer.Time();
			dtTimer.Start();
#endif
			elapsedTime += dt;

			// Correct delta time
			if (DBGoverallTime > elapsedTime)
			{
				DBGoveralTimeError += DBGoverallTime - elapsedTime;
				dt += DBGoverallTime - elapsedTime;
				elapsedTime = DBGoverallTime;
			}

			metronomeTime += dt;
			toneTime += dt;

			if (metronomeTime >= metronomeInterval)
			{
				++metronomeSounded;
				metronomeTime -= metronomeInterval;

				for (uint i = 0; i < numWaveDevs; ++i)
				{
					if (testWaveOut[i].getBufferNum() == 0)
					{
						if (i > maxWaveDev)
							maxWaveDev = i;

						uint rndMetronome = 0;
						testWaveOut[i].play((char *)waveShort, blockSize * c_numChannels * sizeof(short));
						break;
					}
				}
			}

			if (toneTime >= toneInterval)
			{
				++toneSounded;
				toneTime -= toneInterval;

				for (uint i = 0; i < numWaveDevs; ++i)
				{
					if (testWaveOut[i].getBufferNum() == 0)
					{
						if (i > maxWaveDev)
							maxWaveDev = i;

						testWaveOut[i].play((char *)toneWaveShort, toneBlockSize * c_numChannels * sizeof(short));
						break;
					}
				}

				toneFreqTime = (int)(elapsedTime / toneInterval) * toneInterval;
				toneValid = false;
			}
		}

		delete [] wave;
		delete [] waveShort;
	}

	return 0;
}

