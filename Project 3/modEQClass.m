classdef modEQClass
    %MODEQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fc
        Q
        gain
        fs
        
        % Modulation Parameters
        f0
        fd
        lastValue

        % Co Effcients
        a0
        a1
        a2
        b1
        b2
        c0
        d0

        % Delay Variables
        a1_delay
        a2_delay
        b1_delay
        b2_delay
    end
    
    methods
        function obj = modEQClass(fc, Q, gain, f0, fd, lastValue, fs)
            obj.fc = fc;
            obj.Q = Q;
            obj.gain = gain;
            obj.fs = fs;
            obj.f0 = f0;
            obj.fd = fd;
            obj.lastValue = lastValue;
            
        
           
        end
        
        function [output,modValue] = modeffect(obj)
                  
                    inc = obj.f0/obj.fs;
                    modValue = obj.lastValue + inc;
                    if modValue > 1
                        modValue = modValue - 2;
                    elseif modValue < -1
                        modValue = modValue + 2;
                    end
                    saw = modValue * obj.fd;
                    triangleWave = 2 * abs(saw) - 1;
                    output = obj.fc + triangleWave;
                    end

        
        function [outBuffer,outDelays] = process(obj, in_buffer, in_delays)
            obj.a1_delay = in_delays(1);
            obj.a2_delay = in_delays(2);
            obj.b1_delay = in_delays(3);
            obj.b2_delay = in_delays(4);
            x = in_buffer;
            z = zeros(length(in_buffer),1);

            for n = 1:length(x)

                [output,modValue] = obj.modeffect();
                obj.lastValue = modValue;

                % Compute Coefficients
                theta_c = 2*pi*(output/obj.fs);
                V_out = 10^(obj.gain/20);
                zeta = 4/(1+V_out);
                beta = 0.5*((1 - zeta*tan(theta_c/2*obj.Q))/(1 + zeta * tan(theta_c/2*obj.Q)));
                gamma = (0.5 + beta) * cos(theta_c);
                
                obj.a0 = 0.5 - beta;
                obj.a1 = 0;
                obj.a2 = -(0.5 - beta);
                obj.b1 = -2*gamma;
                obj.b2 = 2*beta;
                obj.c0 = V_out - 1.0;
                obj.d0 = 1.0;

                % Difference equation
                z(n) =  x(n)*obj.a0 + obj.a1*obj.a1_delay + obj.a2*obj.a2_delay - obj.b1*obj.b1_delay - obj.b2*obj.b2_delay;
                obj.a2_delay = obj.a1_delay;
                obj.b2_delay = obj.b1_delay;
                obj.a1_delay = x(n);
                obj.b1_delay = z(n);
            end

            outBuffer = obj.d0*x + obj.c0*z;
            outDelays = [obj.a1_delay, obj.a2_delay, obj.b1_delay, obj.b2_delay];
        end
    end
end

