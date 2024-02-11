classdef compressBoi
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        detector % Audio Env Detector Instance
        threshold % In dB
        ratio 
        knee_width
        y_n % Compression Result
        g_n % Gain Diff
        G_n % Final Gain
    end
    
    methods
        function obj = compressBoi(attack, release, threshold, knee_width, ratio)
            
            obj.detector = audio_envelope_detector("RMS", "dB", attack, release);
            obj.threshold = threshold;
            obj.knee_width = knee_width;
            obj.ratio = ratio;
        end
        
        function output = process_audio(obj,input, fs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d_n = obj.detector.process_audio(input, fs);
            d_n = abs(d_n);

            if d_n <= obj.threshold
               obj.y_n = d_n;
               

            elseif d_n > obj.threshold
               obj.y_n = obj.threshold + (d_n - obj.threshold)/obj.ratio;
                
            end
            obj.g_n = obj.y_n - d_n;
            obj.G_n = 10.^(obj.g_n/20);

            output = input .* obj.G_n.*3.5;

        end
    end
end

