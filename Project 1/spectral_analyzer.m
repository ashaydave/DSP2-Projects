
function [magnitude_db,freq_range] = spectral_analyzer(x,fs)
%% MMI - 503/603 Project 1
% Assignment: Develop a spectral analyzer that returns the
% magnitude in decibels and frequency bin values 

% Author : [Ashay Dave]
% Email: [apd122@miami.edu]

%Number of samples
N_FFT = length(x);

%FFT of input signal
x_fft = fft(x); 

%Calculate frequency bin values
freq_binWidth = fs/N_FFT; %comes out to be 1 because fs = N_FFT
freq_range = (0:N_FFT - 1) * freq_binWidth; 

% freq_binValueTest = 0:freq_binWidth:fs-freq_binWidth;

%Calculate magnitude in dB of the FFT
magnitude = abs(x_fft)/N_FFT;
magnitude_db = 20*log10(magnitude);

figure;
plot(freq_range, magnitude_db)
grid on;
xlabel ('Freq Bin Value in Hz')
ylabel ('Magnitude in dB')
title('Magnitude of signal')


end

