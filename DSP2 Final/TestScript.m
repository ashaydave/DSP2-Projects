fs = 44100;
bufferSize = 4096;

[wav, fs] = audioread("Ritviz Remix.wav");
wav = wav(:,1) + wav(:,2);
wav = wav(1:50*fs,1);

wavDelay = longDelayClass(10, fs, 0.8, bufferSize);
wavDelayOutput = wavDelay.process(wav);
soundsc(wavDelayOutput,fs)