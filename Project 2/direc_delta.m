function [delta_out] = direc_delta(delta_index, buffer)

%   Creates a direc delta impulse signal at the choice of index and buffer
%   size which can be initialized before calling the function

delta_out = zeros(buffer,1);
delta_out(delta_index) = 1;

figure;
plot(delta_out);
end

