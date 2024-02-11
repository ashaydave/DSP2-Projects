%% MMI - 503/603 Project 1
% Assignment: Develop a set of audio analysis functions by generating the
% following audio signals: white noise, sine tone @ 1000 Hz, exponential 
% sine sweep and analyze them


% Author : [Ashay Dave]
% Email: [apd122@miami.edu]

%% Guided Portion
% Generate white noise signal
fs = 48000;

sec = 1;
numOfSamples = fs * sec; %duration in samples
noise = 0.2 * randn(numOfSamples, 1);
plot(noise)
title('Noise')
% soundsc(noise)

% noise = wgn(1000,1,0);
% spectrogram(noise, 512, 0, 1024, fs, 'yaxis')


% Generate sine tone @ 1000 Hz
f = 1000;
Ts = 1/fs;
t_vector = [0:Ts:1]; 
y_amp = 1;

y = y_amp * sin(2*pi*f*t_vector);
plot(t_vector, y)
title('Sine Tone @ 1kHz')
xlabel('Time')
ylabel('Amplitude')
% spectrogram(y, 512, 0, 1024, fs, 'yaxis')
% soundsc(y)


% Generate exponential sine sweep (ESS)


f0 = 20; % Start frequency in Hz
f1 = 20000; % End frequency in Hz
Time = 3; % Duration in seconds
t_vector = [0:1/fs:Time]; % Time array

ess = chirp(t_vector,f0,Time,f1,'logarithmic');

plot(t_vector,ess);
xlabel('Time (s)');
ylabel('Amplitude');
title('Exponential Sine Sweep');
filename = 'ess.wav';
audiowrite(filename,y,fs);
% spectrogram(ess, 512, 0, 1024, fs, 'yaxis')
% soundsc(ess)


% Create a separate function file that analyzes the root-mean squared
% loudness in decibels

% RMS of Noise
rms_noise = rms_loudness(noise);
rms_noise 

% RMS of Sine tone
rms_y = rms_loudness(y);
rms_y

%RMS of ESS
rms_ess = rms_loudness(ess);
rms_ess

% Create a separate function file that analyzes the spectral content and
% returns the magnitude in decibels and frequency bin values

% Mag spectrum of Noise
[magnitude_db, freq_range] = spectral_analyzer(noise,fs);

% Mag spectrum  of Sine Tone
[magnitude_db, freq_range] = spectral_analyzer(y,fs);

% Mag spectrum of ESS
[magnitude_db, freq_range] = spectral_analyzer (ess,fs);


%% Analytical directions
% Loop through the signals using a buffer size of 4096 and an overlap of
% 2048 and analyze the loudness and spectral content on each buffer using
% the functions created above. Plot the results over time.

%% SINE TONE @ 1kHz of 2 Seconds

buffer_size = 4096;
overlap = buffer_size/2;

N = length(y);

numOfFrames = floor(2*N/buffer_size-1);

% Store variables for loudness and magnitude
loudness = zeros(1, numOfFrames);
magnitude = zeros(overlap, numOfFrames);
time = [0:numOfFrames - 1] * buffer_size/fs;

% Frame and window audio
for i = 1:numOfFrames

    % Get the Current Frame
    startInd = (i-1) * overlap + 1;
    endInd = startInd + buffer_size - 1;
    currentFrame = y(1, startInd:endInd);
    currentFrame = currentFrame';
    currentFrame = currentFrame.*hann(buffer_size);
    
    % Do the spectral processing of the current frame

    % Calculate loudness
    loudness(1,i) = rms_loudness(currentFrame);

    % Compute FFT of the Current Frame
    fft_CF = fft(currentFrame);

    % Calculate Magnitude
    magnitude(:,i) = 20*log10(abs(fft_CF(1:overlap)));
% 

end

figure;

plot(time,loudness)
xlabel('Time (seconds)')
ylabel('Loudness (dB)')
title('Loudness over time')

figure;

imagesc(time, [0 fs/2], magnitude)
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
title('Spectrum over time')

%% ESS

buffer_size = 4096;
overlap = buffer_size/2;

N = length(ess);

numOfFrames = floor(2*N/buffer_size-1);

loudness = zeros(1, numOfFrames);
magnitude = zeros(overlap, numOfFrames); 
time = [0:numOfFrames - 1] * buffer_size/fs;

% Frame and window audio
for i = 1:numOfFrames

    % Get the Current Frame
    startInd = (i-1) * overlap + 1;
    endInd = startInd + buffer_size - 1;
    currentFrame = ess(1, startInd:endInd);
    currentFrame = currentFrame';
    currentFrame = currentFrame.*hann(buffer_size);
    
    % Do the spectral processing of the current frame

    % Calculate loudness
    loudness(1,i) = rms_loudness(currentFrame);

    % Compute FFT of the Current Frame
    fft_CF = fft(currentFrame);

    % Calculate Magnitude
    magnitude(:,i) = 20*log10(abs(fft_CF(1:overlap)));
% 

end

figure;

plot(time,loudness)
xlabel('Time (seconds)')
ylabel('Loudness (dB)')
title('Loudness over time')

figure;

imagesc(time, [0 fs/2], magnitude)
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
title('Spectrum over time')




%% Creative directions
% Analyze the loudness and sprectrum of a wav file of the your choice both
% the whole signal and per buffer and plot the results over time.

[signal, fs] = audioread("Guitar.wav");
buffer_size = 4096;
overlap = buffer_size/2;

N = length(signal);

numOfFrames = floor(2*N/buffer_size-1);

loudness = zeros(1, numOfFrames);
magnitude = zeros(overlap, numOfFrames); 
time = [0:numOfFrames - 1] * buffer_size/fs;

% Frame and window audio
for i = 1:numOfFrames

    % Get the Current Frame
    startInd = (i-1) * overlap + 1;
    endInd = startInd + buffer_size - 1;
    currentFrame = signal(startInd:endInd);
    currentFrame = currentFrame'; %To perform assignment because the "size of the left side is 1-by-1 and the size of the right side is 1-by-4096."
    currentFrame = currentFrame.*hann(buffer_size);

    % Do the spectral processing of the current frame

    % Calculate loudness
    loudness(1,i) = rms_loudness(currentFrame);

    % Compute FFT of the Current Frame
    fft_CF = fft(currentFrame);

    % Calculate Magnitude
    magnitude(:,i) = 20*log10(abs(fft_CF(1:overlap)));
% 

end
figure;

plot(time,loudness)
xlabel('Time (seconds)')
ylabel('Loudness (dB)')
title('Loudness over time')

figure;

imagesc(time, [0 fs/2], magnitude)
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
title('Spectrum over time')

