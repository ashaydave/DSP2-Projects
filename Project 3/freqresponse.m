function [xfft] = freqresponse(x,fs)

xfft = fft(x,fs);
xfft = xfft(1:end/2);
xfft = db(abs(xfft));
figure;
plot(xfft)
title("Freq Response")
xlabel('Frequency (kHz)')
ylabel('Magnitude(dB)');
end

