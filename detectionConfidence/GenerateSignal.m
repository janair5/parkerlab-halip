function signal=GenerateSignal(StimulusSettings, NoiseAttFactor, EmbedSignal)

global BpodSystem %we need this for volume adjustment

%% abbreviate variable names and clip impossible values for better handling
SamplingRate=StimulusSettings.SamplingRate;
SignalRamp=StimulusSettings.Ramp; %UPDATE HERE IF NO NOISE IS USED
SignalDuration=StimulusSettings.SignalDuration;
SignalForm=StimulusSettings.SignalForm;
SignalMinFreq=StimulusSettings.SignalMinFreq;
SignalMaxFreq=StimulusSettings.SignalMaxFreq;
SignalVolume=max(min(StimulusSettings.SignalVolume,StimulusSettings.MaxVolume),StimulusSettings.MinVolume);%clip signal volume to Min and Max

    %make signalclose all
    
    t=linspace(0,SignalDuration,[SamplingRate*SignalDuration]); %time vector for chirp
    switch SignalForm
        case 'LinearUpsweep'
            signal=chirp(t,SignalMinFreq,SignalDuration,SignalMaxFreq);
            % signal = GenerateSineSweep(SamplingRate, SignalMinFreq, SignalMaxFreq, SignalDuration);
            % disp('### ALERT ### Audio sweep generated using GenerateSineSweep. Return to chirp() before production')
            freqvec=SignalMinFreq+(SignalMaxFreq-SignalMinFreq)*t;
        case 'LinearDownsweep' %gaussian noise from mean 0 std .25
            signal=chirp(t,SignalMaxFreq,SignalDuration,SignalMinFreq);
            freqvec=SignalMaxFreq+(SignalMinFreq-SignalMaxFreq)*t;
        case 'QuadraticConvex'
            tnew=t-mean(t);
            signal=chirp(tnew,SignalMinFreq,SignalDuration./2,SignalMaxFreq,'quadratic',[],'convex'); %make chirp
            freqvec=SignalMinFreq+(SignalMaxFreq-SignalMinFreq)./tnew(1)*tnew.^2;
    end
    
    signal=[signal;signal];
    
    %adjust signal volume
    SoundCal = BpodSystem.CalibrationTables.SoundCal;
    if(isempty(SoundCal))
        disp('Error: no sound calibration file specified');
        return
    end
    if size(SoundCal,2)<2
        disp('Error: no two speaker sound calibration file specified');
        return
    end
    
    for s=1:2 %loop over two speakers
        toneAtt = polyval(SoundCal(1,s).Coefficient,freqvec);%Frequency dependent attenuation factor with less attenuation for higher frequency (based on calibration polynomial)
        %toneAtt = [polyval(SoundCal(1,1).Coefficient,toneFreq)' polyval(SoundCal(1,2).Coefficient,toneFreq)']; in Torben's script
        diffSPL = SignalVolume - [SoundCal(1,s).TargetSPL];
        attFactor = sqrt(10.^(diffSPL./10)); %sqrt(10.^(diffSPL./10)) in Torben's script WHY sqrt?
        att = toneAtt.*attFactor;%this is the value for multiplying signal scaled/clipped to [-1 to 1]
        signal(s,:)=signal(s,:).*att; %should the two speakers dB be added?
    end
    
%put an envelope to avoide clicking sounds at beginning and end
omega=(acos(sqrt(0.1))-acos(sqrt(0.9)))/(SignalRamp/pi*2); % This is for the envelope with Ramp duration duration
t=0 : (1/SamplingRate) : pi/2/omega;
t=t(1:(end-1));
RaiseVec= (cos(omega*t)).^2;

Envelope = ones(length(signal),1); % This is the envelope
Envelope(1:length(RaiseVec)) = fliplr(RaiseVec);
Envelope(end-length(RaiseVec)+1:end) = (RaiseVec);

signal = signal.*Envelope';

% Add noise to match HiFi module's noise generator
BGNoise = (1-(rand(size(signal))*2))*NoiseAttFactor;

if EmbedSignal
    signal = signal + BGNoise;
else
    signal = BGNoise;
end
