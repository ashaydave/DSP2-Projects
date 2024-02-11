    
function rms_db = rms_loudness(x)
%% MMI - 503/603 Project 1
% Assignment: Develop an rms analyzer that returns the root-mean square in
% decibels
% Author : [Ashay Dave]
% Email: [apd122@miami.edu]

%Calculate RMS of audio
y_rms = rms(x);

%Convert to dB
rms_db = 20*log10(y_rms);

end

