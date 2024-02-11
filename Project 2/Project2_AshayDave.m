%% MMI - 503/603 Project 2
% Follow the directions and conduct the MATLAB Examples directed in the
% assignment documet


% Author : Ashay Dave
% Email: [apd122@miami.edu]

%% GUIDED PORTION 
% 1.a.iii - How to apply filters to a signal

%   Create a dirac-delta signal (do not put the impulse and index 1)
fs = 48000;
fn = fs/2; %Nyquist
buffer = 4096;
delta = direc_delta(buffer/2, buffer);

%   Convert the dirac-delta function to the frequency domain using the FFT

delta_fft = fft(delta,buffer);
figure;
plot(abs(delta_fft(1:(buffer/2))))


%   Create a low-pass filter by attenuating the frequencies from 
%   3/4*nyquist to
%   nyquist

slope = linspace(1,0, (buffer - ((buffer/2) * (3/4)) + 1));
slope = complex(slope);
delta_fft((buffer/2) * (3/4):end) = delta_fft((buffer/2) * (3/4):end).*slope';

figure;
plot(abs(delta_fft(1:buffer/2)));
title("Low Pass Filter")

% % length_of_sample = ((buffer/2) - ((3/4) * buffer/2));
% % for i = (((3/4) * buffer/2) : buffer/2)
% % delta_fft(i,1) = abs(delta_fft(i-1,1)) - (1/length_of_sample);
% % end
% % for i = (buffer/2) + 1 :(buffer/2) + length_of_sample
% %     delta_fft(i,1) = abs(delta_fft(i-1,1)) + 1/length_of_sample;
% % end
% % figure;
% % plot(abs(delta_fft(1:(buffer/2))))


%   Convert the new impulse response to the time domain by using IFFT

lowpass_ir = ifft(delta_fft,buffer, 'symmetric');
figure;
plot(lowpass_ir);
xlabel('Amplitude');
ylabel('Frequency (samples)');
title('LPF IR');
%   Use the built in convolution function to time-based convolve the 
%   original dirac-delta signal with the new impulse respone (make sure the
%   output is the same length as the dirac-delta signal, see MATLAB help 
%   for this).

delta_lpf_time = conv(delta,lowpass_ir,"same");
figure;
plot(delta_lpf_time);
ylabel('Magnitude')
xlabel('Frequency (samples)')
title("Time Based Convolution")

%   Take both dirac-delta signal and the impulse response into the 
%   frequency domain and multiply them to do the frequency-based 
%   convolution. Take the IFFT of the resu  lt into the time domain.

delta1 = direc_delta(buffer/2, buffer);
delta_fft_freq = fft(delta1,buffer); %new variable for delta signal in frequency domain
dd_freqconv = delta_fft_freq .* delta_fft; %'delta_fft' here is the lowpass ir
delta_lpf_freq = ifft(dd_freqconv,buffer,'symmetric');

figure;
plot(delta_lpf_freq);
ylabel('Magnitude')
xlabel('Frequency (samples)')
title("Freq Based Convolution")


%   Plot the magnitude response of the time-based convolution and the 
%   frequency-based convolution in the same plot to compare them.

figure;
subplot(2,1,1)
plot(delta_lpf_time);
ylabel('Magnitude')
xlabel('Frequency (samples)')
title("Time Based Convolution")

subplot(2,1,2)
plot(delta_lpf_freq);
ylabel('Magnitude')
xlabel('Frequency (samples)')
title("Freq Based Convolution")



%   1.b.iv - Compute the coefficients for the shelving filter provided. 
%   The angles are 0 and place the pole at 0.3 on the real axis and the zero 
%   at 0.7 on the real axis. 

rs_zero = 0.7;
rs_pole = 0.3;
theta = 0;

a0_s = 1;
a1_s = -rs_zero*cos(theta);
b1_s = -rs_pole*cos(theta);

%   1.b.v - Apply the dirac delta function to the difference equation using
%   the coefficients computed above. Plot the magnitude response of the
%   output of the signal.

buffer = 4096;
x = direc_delta(buffer/2, buffer); %delta signal

y = zeros(length(buffer),1); %output

a1_delay = 0;
b1_delay = 0;

for n = 1:length(x)
    
    y(n) = a0_s*x(n) + a1_s*a1_delay - b1_s*b1_delay;
    a1_delay = x(n);
    b1_delay = y(n);
end


y_l = length(y); % Length of output signal
Y = fft(y); % Compute FFT of output signal

% Magnitude response of output signal
magY = (abs(Y));


% Plot the magnitude response
figure;
plot(magY);
xlim([0 length(magY)/2])
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');



%% Advanced Portion

% 2.d - Compute the coefficients with the following poles and zeros of a
% biquad filter. Use the angles Pi/4 and –Pi/4 for the zeros and ¾*Pi and 
% –¾*Pi for the poles. Use an R value of 0.5 for the Zeros and 0.25 for 
% the Poles. 

rb_zero = 0.5;
rb_pole = 0.25;
theta_zero = pi/4;
theta_pole = (3*pi)/4;

a0_b = 1;
a1_b = -2*rb_zero*cos(theta_zero);
a2_b = rb_zero * rb_zero;
b1_b = -2*rb_pole*cos(theta_pole);
b2_b = rb_pole * rb_pole;


% 2.e - Apply the coefficients and use the dirac-delta function to test the
% filter. Take the FFT of the output and plot the magnitude. 

buffer1 = 4096;
x1 = direc_delta(buffer1/2, buffer1); %delta signal

y1 = zeros((buffer1),1); %output

a1b_delay = 0;
b1b_delay = 0;
a2b_delay = 0;
b2b_delay = 0;

for n = 1:length(y1)
    y1(n) = a0_b*x1(n) + a1_b*a1b_delay + a2_b*a2b_delay - b1_b*b1b_delay - b2_b*b2b_delay;
    a2b_delay = a1b_delay;
    b2b_delay = b1b_delay;
    a1b_delay = x1(n);
    b1b_delay = y1(n);
    
end


y1_l= length(y1); % Length of output signal
Y1 = fft(y1); % Compute FFT of output signal

% Magnitude response of output signal
magY1 = (abs(Y1));


% Plot the magnitude response
figure;
plot(magY1);
xlim([0 length(magY1)/2])
xlabel('Frequency');
ylabel('Magnitude)');
title('Magnitude Response of Biquad');

%% Creative Portion

% 3.a - pick a song or a track, plot the signal, and the spectrogram

[audio, fs] = audioread('Shaanti.wav');
audio = audio(:,1) + audio(:,2);

% Plot audio signal
t = linspace(0, length(audio)/fs, length(audio));
figure;
plot(t, audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Audio Signal');

window = 512;
overlap = 256;
nfft = 1024;
figure;
spectrogram(audio,window,overlap,nfft,fs,'yaxis')

length_audio = length(audio);
frv = (0:length_audio - 1)*(fs/length_audio);
plot(frv, abs(fft(audio)))
title('FFT of Audio File')
xlabel('Frequency in kHz')
ylabel('Magnitude')



% 3.b - Choose a new zero and pole for the shelving filter and pass the
% song through it. Plot the output and the spectrogram (you can use the
% built-in function for this).

rs_z = 0.9; %Low shelf
rs_p = 0.2;
%theta = 0

a0_crs = 0.9;
a1_crs = -rs_z*cos(theta); 
b1_crs = -rs_p*cos(theta);

buffer = 4096;
x1_a = audio;

y1_a = zeros(length(buffer),1);
a1_delay1 = 0;
b1_delay1 = 0;

for n = 1:length(x1_a)
    
    y1_a(n) = a0_crs*x1_a(n) + a1_crs*a1_delay1 - b1_crs*b1_delay1;
    a1_delay1 = x1_a(n);
    b1_delay1 = y1_a(n);
end

window = 512;
overlap = 256;
nfft = 1024;
figure;
spectrogram(y1_a,window,overlap,nfft,fs,'yaxis')


length_y1_a = length(y1_a);
frv1 = (0:length_y1_a - 1)*(fs/length_y1_a);
plot(frv1, abs(fft(y1_a)))
title('FFT of Audio File through Low Shelf')
xlabel('Frequency in kHz')
ylabel('Magnitude')

%soundsc(y1_a,fs)

% 3.c - Choose a new set of zeros and poles for the biquad filter and pass
% the song through it. Plot the output and the spectrogram (you can use the
% built-in function for this).



%reloading audio signal
[audio2, fs] = audioread('Shaanti.wav');
audio2 = audio2(:,1) + audio2(:,2);

buffer2 = 4096;
x2_a = audio2;

rb_z = 0.9;
rb_p = 0.1;
theta_bz = pi/2;
theta_bp = (3*pi)/4;

a0_crb = 1;
a1_crb = -2*rb_z*cos(theta_bz);
a2_crb = rb_z*rb_z;
b1_crb = -2*rb_p*cos(theta_bp);
b2_crb = rb_p * rb_p;

y2_a = zeros(size(x2_a)); %output

a1b_delay2 = 0;
b1b_delay2 = 0;
a2b_delay2 = 0;
b2b_delay2 = 0;

for n = 1:length(y2_a)
    y2_a(n) = a0_crb*x2_a(n) + a1_crb*a1b_delay2 + a2_crb*a2b_delay2 - b1_crb*b1b_delay2 - b2_crb*b2b_delay2;
    a2b_delay2 = a1b_delay2;
    b2b_delay2 = b1b_delay2;
    a1b_delay2 = x2_a(n);
    b1b_delay2 = y2_a(n);

end

length_y2_a = length(y2_a);
frv2 = (0:length_y2_a - 1)*(fs/length_y2_a);
plot(frv2, abs(fft(y2_a)))
title('FFT of Audio File through Biquad')
xlabel('Frequency in kHz')
ylabel('Magnitude')


window = 512;
overlap = 256;
nfft = 1024;
figure;
spectrogram(y2_a,window,overlap,nfft,fs,'yaxis')



