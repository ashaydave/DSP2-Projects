#ifndef ASHAY_CBDELAY_H
#define ASHAY_CBDELAY_H

/*
  ==============================================================================

    DelayedMehndi.h
    Created: 29 April 2023 8:57:21pm
    Author:  ashay

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include <math.h>

namespace Ashay
{
    class DelayedMehndi
    {
    public:
        DelayedMehndi();
        ~DelayedMehndi();
        void prepareToPlay(double sampleRate, int samplesPerBlock);
        void process(juce::AudioBuffer<float>& buffer, float delayLengthParameter, float dryWetParameter, float feedbackParameter);
        void reset();

    private:    
        juce::AudioBuffer<float> delayBuffer;
        int delayBufferLength, delayReadIndex, delayWriteIndex, counter;
        float delayLength, dryMix, wetMix, feedback;
        float delayLengthForBuffer;
        double delaySampleRate;
        int delaySamplesPerBlock;
    };
}

#endif