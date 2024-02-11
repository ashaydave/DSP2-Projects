/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
DelayedMehndiAudioProcessorEditor::DelayedMehndiAudioProcessorEditor (DelayedMehndiAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (600, 200);

    auto parameterTree = audioProcessor.getParameters();

    addAndMakeVisible(delayTimeLabel);
    delayTimeLabel.setText("Time", juce::dontSendNotification);
    delayTimeLabel.setColour(juce::Label::textColourId, juce::Colours::mediumvioletred);
    delayTimeLabel.setJustificationType(juce::Justification::centred);
    delayTimeLabel.setFont(juce::Font(14.0f, juce::Font::bold));

    addAndMakeVisible(delayTimeSlider);
    delayTimeSlider.setSliderStyle(juce::Slider::SliderStyle::RotaryVerticalDrag);
    delayTimeSlider.setColour(juce::Slider::ColourIds::rotarySliderFillColourId, juce::Colours::mediumvioletred);
    delayTimeSlider.setColour(juce::Slider::ColourIds::backgroundColourId, juce::Colours::black);
    delayTimeSlider.setTextBoxStyle(juce::Slider::NoTextBox, true, 50, 50);

    delayTimeSliderParameter = (juce::AudioParameterFloat*)parameterTree.getUnchecked(0);
    delayTimeSlider.setRange(delayTimeSliderParameter->range.start, delayTimeSliderParameter->range.end, 0.0001);

    delayTimeSlider.setValue(*delayTimeSliderParameter, juce::NotificationType::dontSendNotification);

    delayTimeSlider.onValueChange = [this]
    {
        *delayTimeSliderParameter = delayTimeSlider.getValue();
    };

    addAndMakeVisible(dryWetLabel);
    dryWetLabel.setText("Dry/Wet", juce::dontSendNotification);
    dryWetLabel.setColour(juce::Label::textColourId, juce::Colours::forestgreen);
    dryWetLabel.setJustificationType(juce::Justification::centred);
    dryWetLabel.setFont(juce::Font(14.0f, juce::Font::bold));
   
    addAndMakeVisible(dryWetSlider);
    dryWetSlider.setSliderStyle(juce::Slider::SliderStyle::RotaryVerticalDrag);
    dryWetSlider.setColour(juce::Slider::ColourIds::rotarySliderFillColourId, juce::Colours::forestgreen);
    dryWetSlider.setColour(juce::Slider::ColourIds::backgroundColourId, juce::Colours::black);
    dryWetSlider.setTextBoxStyle(juce::Slider::NoTextBox, true, 50, 50);

    dryWetSliderParameter = (juce::AudioParameterFloat*)parameterTree.getUnchecked(1);
    dryWetSlider.setRange(dryWetSliderParameter->range.start, dryWetSliderParameter->range.end, dryWetSliderParameter->range.interval);

    dryWetSlider.setValue(0, juce::NotificationType::dontSendNotification);

    dryWetSlider.onValueChange = [this]
    {
        *dryWetSliderParameter = dryWetSlider.getValue();
    };

    addAndMakeVisible(feedbackLabel);
    feedbackLabel.setText("Feedback", juce::dontSendNotification);
    feedbackLabel.setColour(juce::Label::textColourId, juce::Colours::deepskyblue);
    feedbackLabel.setJustificationType(juce::Justification::centred);
    feedbackLabel.setFont(juce::Font(14.0f, juce::Font::bold));
    
    addAndMakeVisible(feedbackSlider);
    feedbackSlider.setSliderStyle(juce::Slider::SliderStyle::RotaryVerticalDrag);
    feedbackSlider.setColour(juce::Slider::ColourIds::rotarySliderFillColourId, juce::Colours::deepskyblue);
    feedbackSlider.setColour(juce::Slider::ColourIds::backgroundColourId, juce::Colours::black);
    feedbackSlider.setTextBoxStyle(juce::Slider::NoTextBox, true, 50, 50);

    feedbackSliderParameter = (juce::AudioParameterFloat*)parameterTree.getUnchecked(2);
    feedbackSlider.setRange(feedbackSliderParameter->range.start, feedbackSliderParameter->range.end, feedbackSliderParameter->range.interval);

    feedbackSlider.setValue(0, juce::NotificationType::dontSendNotification);

    feedbackSlider.onValueChange = [this]
    {
        *feedbackSliderParameter = feedbackSlider.getValue();
    };
}

DelayedMehndiAudioProcessorEditor::~DelayedMehndiAudioProcessorEditor()
{
}

//==============================================================================
void DelayedMehndiAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll(juce::Colours::black);

  /*  g.setColour (juce::Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Hello World!", getLocalBounds(), juce::Justification::centred, 1);*/
}

void DelayedMehndiAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..

    auto x = this->getLocalBounds().getWidth() / 10;
    auto y = this->getLocalBounds().getHeight() / 10;

    delayTimeLabel.setBounds(100, 160, 100, 20);
    dryWetLabel.setBounds(250, 160, 100, 20);
    feedbackLabel.setBounds(400, 160, 100, 20);

    delayTimeSlider.setBounds(0, 2.5 * y, 5 * x, 5 * y);
    dryWetSlider.setBounds(2.5 * x, 2.5 * y, 5 * x, 5 * y);
    feedbackSlider.setBounds(5 * x, 2.5 * y, 5 * x, 5 * y);
}
