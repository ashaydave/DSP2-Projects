function [output,modValue] = modulatorBoi(f0, fs, fc, lastValue, fd)
%MODULATOR Summary of this function goes here
%   Detailed explanation goes here
inc = f0/fs;
modValue = lastValue + inc;
if modValue > 1
    modValue = modValue - 2;
elseif modValue < -1
    modValue = modValue + 2;
end
saw = modValue * fd;
triangleWave = 2 * abs(saw) - 1;
output = fc + triangleWave;
end
