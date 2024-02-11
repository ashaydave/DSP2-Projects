classdef audio_envelope_detector
    %AUDIO_ENVELOPE_DETECTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rectifier_type; %MS or RMS or normal
        sample_type; %dB or normal
        attack_time; %in ms
        release_time; %in ms
    end
    
    methods
        function obj = audio_envelope_detector(rectifier_type,sample_type, attack_time, release_time)
            %AUDIO_ENVELOPE_DETECTOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.rectifier_type = rectifier_type;
            obj.sample_type = sample_type;
            obj.attack_time = attack_time;
            obj.release_time = release_time;
        end
        
        function [envelope_values] = process_audio(obj,input_buffer, fs)
            %Process Audio Summary of this method goes here
            %   Detailed explanation goes here

            % rectifier
            if strcmp(obj.rectifier_type, "MS") == 1 || strcmp(obj.rectifier_type, "RMS") == 1 
                    rect_input = abs(input_buffer).^2;
            elseif strcmp(obj.rectifier_type, "normal") == 1
                    rect_input = abs(input_buffer);
            else
                    error("rectifier_type must be normal, MS, or RMS");
            end

            % RC simulator
            % compute coefficients of filter based on attack time and
            % release time
            Tc = log(0.368);
            a0 = Tc/(exp(1)^(fs*obj.attack_time*0.001));
            b1 = Tc/(exp(1)^(fs*obj.release_time*0.001));
            rc_buff = zeros(size(input_buffer));
            forward_val = 0;
            delay_val = 0;
            for i = 1:length(input_buffer)
                forward_val = a0*(rect_input(i)+ b1*delay_val);
                rc_buff(i) = rect_input(i) + forward_val;
                delay_val = rc_buff(i);
            end

            %RMS and log conversion
            if strcmp(obj.sample_type, "normal") == 1
                    envelope_values = rc_buff;
            elseif strcmp(obj.sample_type, "dB") == 1
                    envelope_values = 20*log10(sqrt(rc_buff));
            else
                    error("sample_type must be normal or dB");
            end

        end
    end
end

