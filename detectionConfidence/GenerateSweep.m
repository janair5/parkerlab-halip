function wave = GenerateSweep(SamplingRate, StartFreq, EndFreq, Duration)
    t=1/SamplingRate:1/SamplingRate:1;    
    k=(EndFreq-StartFreq)/(Duration-t(1));
    wave=cos(2*pi*(k/2*t+StartFreq).*t+(pi/2));
end