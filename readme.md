# Sound synthesis sample
This project uses shared Core project (via git submodules) - in case of offline project assembly, please manually sync the Core repository and paste files into the `core` folder.

## Sample description
This sample shows how to use the sound synthesis framework. Again, modern neural networks make at least some parts of the sound synth obsolete, but the framework could still be used as an example of how the real-time DSP was done back in the day!

Examples include:

[Complex ambient synthesis sample](examples/ambient.mp3) which includes detuned oscillators, reverb, chorus, delay lines, general IIR filtering and a simple sequencer.

[Sound positioning and Doppler effect sample.](examples/doppler.mp3)

[Extended Karplus-Strong synthesis](examples/ks_horizons.mp3), including sequencing via guitar strings and frets setup.

[Formant vowel synthesis.](examples/formant.mp3)

[Phaser/flanger wind sound synth.](examples/wind.mp3)

## Underlying sound synthesis framework
The framework is a product of old demoscene days, hence it is designed to be fast and small, and doesn't work with the sound loading routines, although could be easily modified to support those too.

**Key** sound manipulations provided:
1. ADSR (Attack-Decay-Sustain-Release amplitude modulator)
2. Generic 2nd order IIR filters (and specific pole-zero configurations of lowpass, highpass, bandpass, notch)
3. Specialized feedback-delay filters: allpass, comb
4. Delay lines
5. Derived processing effects: chorus, flanger, phaser
6. Generators: oscillators and pink noise
7. 4-channel FDN reverb
8. WaveOut playback
9. Spatial positioning for interactive environments:
    1. Doppler effect
    2. Stereo positioning: interaural time difference (ITD), interaural intensity difference (IID), slight touch of lowpass to simulate head shadow; didn't quite get to the full 3D spatialization via HRTF approximation that I had planned


## Acknowledgments
Huge thanks to Yehar for his introductory sound synthesis tutorial, to all authors who published in Hugi demoscene diskmag (check out the coding special edition if you're interested in the demoscene-related knowledge!); and to all the demoscene community which is an open and friendly!

Additionally, Julius Orion Smith III for his wonderful sound processing resources, and to Christian Schüler for the popularization of FDN-based reverb algorithms.  

## License
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License
