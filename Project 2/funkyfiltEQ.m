function [out_buffer, out_delays] = funkyfiltEQ(in_buffer, f_center, Q, gain, in_delays, fs)
% Constant Q Parametric EQ and filter

% Author: Dr. Justin D. Mathew, Ph.D
% Email: jdm423@miami.edu

K = tan(pi * f_center/fs);
V_out = 10^(gain/20);
d0 = 1+ (1/Q)*K +K^2;
e0 = 1 + 1/(V_out*Q)*K + K^2;

alpha = 1 + (V_out/Q)*K + K^2;
beta = 2 * (K^2 - 1);
gamma =  1 - (V_out/Q)*K + K^2;
delta = 1 - (1/Q)*K + K^2;
eta = 1-(1/(V_out*Q))*K+K^2;

if (gain > 0)
    a0=alpha/d0;
    a1 = beta/d0;
    a2 = gamma/d0;
    b1 = beta/d0;
    b2 = delta/d0;
    c0 = 1.0;
    d0= 0.0;

elseif (gain < 0)
    a0 = d0/e0;
    a1 = beta/e0;
    a2 = delta/e0;
    b1 = beta/e0;
    b2 = eta/e0;
    c0 = 1.0;
    d0 = 0.0;

else
    out_buffer = in_buffer;
    return;
end

out_buffer = zeros(length(in_buffer),1);
a1_delay = in_delays(1);
a2_delay = in_delays(2);
b1_delay = in_delays(3);
b2_delay = in_delays(4);

x = in_buffer;
z = zeros(length(in_buffer),1);

for n = 1:length(out_buffer)
    
    z(n) =  (x(n)*a0 + ...
        a1*a1_delay + a2*a2_delay - b1*b1_delay - b2*b2_delay);
    a2_delay = a1_delay;
    b2_delay = b1_delay;
    a1_delay = x(n);
    b1_delay = z(n);

end

out_buffer = d0*x + c0*z;
out_delays = [a1_delay, a2_delay, b1_delay, b2_delay];

end

