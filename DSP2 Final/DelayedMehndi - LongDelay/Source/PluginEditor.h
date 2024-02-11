/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class DelayedMehndiAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    DelayedMehndiAudioProcessorEditor (DelayedMehndiAudioProcessor&);
    ~DelayedMehndiAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    DelayedMehndiAudioProcessor& audioProcessor;

    juce::Label delayTimeLabel, dryWetLabel, feedbackLabel;
    juce::Slider delayTimeSlider, dryWetSlider, feedbackSlider;

    juce::AudioParameterFloat* delayTimeSliderParameter;
    juce::AudioParameterFloat* dryWetSliderParameter;
    juce::AudioParameterFloat* feedbackSliderParameter;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DelayedMehndiAudioProcessorEditor)
};
