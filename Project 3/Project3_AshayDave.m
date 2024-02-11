%% NON CONSTANT Q WITH TRIANGLE WAVE MODULATION %%
%% CONSTANTS & EFFECT PARAMETERS

fs = 44100;
buffer_size = 4096;

% Defining EQ parameters
fc = 1000;
Q = 0.5;
gain = 12.0;
inDelays = [0,0,0,0];

% Modulating Parameters

f0 = 5;
fd = 1000;
lastValue = 0;

% Defining spectrogram parameters
window = 512;
overlap = 256;
nfft = 1024;

%% SINE TONE

f = 1000;
Ts = 1/fs;
t_vector = 0:Ts:3; 
amp = 1;

sine = amp * sin(2*pi*f*t_vector);

sineBuffer = buffer(sine, buffer_size);
sineOutput = zeros(size(sine));

sineEQ = modEQClass(fc, Q, gain, f0, fd, lastValue, fs);

for i = 1:size(sineBuffer, 2)
    [sine_chunk, inDelays] = sineEQ.process(sineBuffer(:, i),inDelays);
    start_index = (i-1)*buffer_size + 1;
    end_index = start_index + buffer_size - 1;
    sineOutput(start_index:end_index) = sine_chunk;
end

soundsc(sineOutput, fs);

figure;
spectrogram(sine,window,overlap,nfft,fs,'yaxis')
figure;
spectrogram(sineOutput,window,overlap,nfft,fs,'yaxis')

% sineFFT = freqresponse(sine,fs);
% sineOutputFFT = freqresponse(sineOutput,fs);

%% ESS

f1 = 20; 
f2 = 20000;
Time = 5;
t_vector_ess = 0:Ts:Time;

ess = chirp(t_vector_ess,f1,Time,f2,'logarithmic');

essBuffer = buffer(ess, buffer_size);
essOutput = zeros(size(ess));

essEQ = modEQClass(fc, Q, gain, f0, fd, lastValue, fs);

for i = 1:size(essBuffer, 2)
    [essChunk, inDelays] = essEQ.process(essBuffer(:, i), inDelays);
    start_index = (i-1)*buffer_size + 1;
    end_index = start_index + buffer_size - 1;
    essOutput(start_index:end_index) = essChunk;
end

soundsc(essOutput, fs);

figure;
spectrogram(ess,window,overlap,nfft,fs,'yaxis')
figure;
spectrogram(essOutput,window,overlap,nfft,fs,'yaxis')

% essFFT = freqresponse(ess,fs);
% essOutputFFT = freqresponse(essOutput,fs);


%% AUDIO

[wav, fs] = audioread("Ritviz Remix.wav");
wav = wav(:,1) + wav(:,2);
wav = wav(1:50*fs,1);

wavBuffers = buffer(wav, buffer_size);
wavOutput = zeros(size(wavBuffers, 1)*size(wavBuffers,2),1);

wavEQ = modEQClass(fc, Q, gain, f0, fd, lastValue, fs);

for i = 1:size(wavBuffers, 2)
    start_ind = (i-1)*4096 + 1; 
    stop_ind = start_ind + 4095;
    [wavOutput(start_ind:stop_ind,1), inDelays] = wavEQ.process(wavBuffers(:,i), inDelays);

end

soundsc(wavOutput, fs);

figure;
spectrogram(wav,window,overlap,nfft,fs,'yaxis')
figure;
spectrogram(wavOutput,window,overlap,nfft,fs,'yaxis')

% wavFFT = freqresponse(wav,fs);
% wavOutputFFT = freqresponse(wavOutput,fs);

%% LONG DELAY %%

%% SINE TONE

sineDelay = longDelayClass(1, fs, 0.5, buffer_size);
sineDelayOutput = sineDelay.process(sine);
soundsc(sineDelayOutput,fs)

figure;
spectrogram(sineDelayOutput,window,overlap,nfft,fs,'yaxis')



%% ESS

essDelay = longDelayClass(1, fs, 0.5, buffer_size);
essDelayOutput = essDelay.process(ess);
soundsc(essDelayOutput,fs)

figure;
spectrogram(essDelayOutput,window,overlap,nfft,fs,'yaxis')


%% AUDIO

wavDelay = longDelayClass(1, fs, 0.5, buffer_size);
wavDelayOutput = wavDelay.process(wav);
soundsc(wavDelayOutput,fs)

figure;
spectrogram(wavDelayOutput,window,overlap,nfft,fs,'yaxis')


%% COMPRESSOR %%

%% SINE TONE

sineCompress = compressBoi(0.5, 1, -10, 10, 50);
sineCompressOutput = sineCompress.process_audio(sine,fs);

sineCompressOutput = fillmissing(sineCompressOutput,'linear'); % to solve pwelch error "expected x to be finite" - replacing NaN values with linearly interpolated values.

plot(sine)
figure;
plot(sineCompressOutput)

figure;
spectrogram(sineCompressOutput,window,overlap,nfft,fs,'yaxis');

soundsc(sineCompressOutput,fs)

%% ESS

essCompress = compressBoi(0.5, 1, -10, 10, 50);
essCompressOutput = essCompress.process_audio(ess,fs);

plot(ess)
hold on
plot(essCompressOutput)

figure;
spectrogram(essCompressOutput,window,overlap,nfft,fs,'yaxis')

soundsc(essCompressOutput,fs)

%% AUDIO

wavCompress = compressBoi(0.5, 1, -10, 10, 50);
wavCompressOutput = wavCompress.process_audio(wav,fs);

plot(wav)
hold on
plot(wavCompressOutput)

figure;
spectrogram(wavCompressOutput,window,overlap,nfft,fs,'yaxis')

soundsc(wavCompressOutput,fs)