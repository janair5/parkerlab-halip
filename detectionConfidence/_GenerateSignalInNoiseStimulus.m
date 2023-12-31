function [stimulus,embed]=GenerateSignalInNoiseStimulus(StimulusSettings)
%makes auditory noise and amplitude modulation envelope (to be implemented) with the following properties
%NoiseSettings.Dur - Duration of Noise in seconds
%NoiseSettings.SamplingRate - Sampling Rate of Sound Card
%NoiseSettings.Color - Type of Noise 'white_uniform', 'white_gaussian' or 'pink_gaussian'


global BpodSystem %we need this for volume adjustment

%% abbreviate variable names and clip impossible values for better handling
SamplingRate=StimulusSettings.SamplingRate;
Ramp=StimulusSettings.Ramp;
EmbedSignal=StimulusSettings.EmbedSignal;
NoiseColor=StimulusSettings.NoiseColor;
NoiseDuration=StimulusSettings.NoiseDuration;
NoiseVolume=StimulusSettings.NoiseVolume;
SignalDuration=StimulusSettings.SignalDuration;
SignalForm=StimulusSettings.SignalForm;
SignalMinFreq=StimulusSettings.SignalMinFreq;
SignalMaxFreq=StimulusSettings.SignalMaxFreq;
SignalVolume=StimulusSettings.SignalVolume;

%check whether Ramp is reasonable
if 2*Ramp>.5*NoiseDuration
   error('Your stimulus ramp is too long! Check your stimulus settings.')
end

%% generate noise
%generate noise vector
samplenum=round(SamplingRate * NoiseDuration);
switch NoiseColor
    case 'WhiteUniform'
        noise = 2 * rand(1, samplenum) - 1;%make white uniform noise -1 to 1
    case 'WhiteGaussian' %gaussian noise from mean 0 std .25
        noise = .25 * randn(1, samplenum);%make white gausian noise and clip to [-1 1]
        noise(noise<-1)=-1;
        noise(noise>1)=1;
    case 'PinkGaussian'
        noise = f_alpha_gaussian (samplenum,.2^2,1)';%make colored gaussian noise with std .2 and 1/f^alpha
        %     case 'naturalistic' %https://www.ncbi.nlm.nih.gov/pubmed/18301738 "Efficient coding of naturalistic stimuli"
        %         noise = randn(1, samplenum);%
        %         noise1 = f_alpha_gaussian ( samplenum,1, 0.1)';%make colored noise  1/f^alpha
        %         noise2 = f_alpha_gaussian ( samplenum,1, 0.1)';%make colored noise 1/f^alpha
        %         noisemodulator=(noise1.^2+noise2.^2).^.1;%make colored noise modulator with amplitude modulation spectrume Raleigh
end

%clip noise to [-1 1]
noise(noise<-1)=-1;
noise(noise>1)=1;

%put in double speaker
noise = [noise;noise];

%adjust noise volume
SoundCal = BpodSystem.CalibrationTables.SoundCal;
if(isempty(SoundCal))
    disp('Error: no sound calibration file specified');
    return
    endc
toneAtt = mean(polyval(SoundCal(1,1).Coefficient,linspace(SignalMinFreq,SignalMaxFreq))); %just take the mean over signal frequencies -
=======
if size(SoundCal,2)<2
   disp('Error: no stereo sound calibration file specified') 
end
toneAtt = [mean(polyval(SoundCal(1,1).Coefficient,linspace(SignalMinFreq,SignalMaxFreq))),...
    mean(polyval(SoundCal(1,2).Coefficient,linspace(SignalMinFreq,SignalMaxFreq)))]; %just take the mean over signal frequencies -
>>>>>>> develop
%toneAtt = [polyval(SoundCal(1,1).Coefficient,toneFreq)' polyval(SoundCal(1,2).Coefficient,toneFreq)']; in Torben's script
diffSPL = NoiseVolume - [SoundCal(1,1).TargetSPL];
attFactor = sqrt(10.^(diffSPL./10)); %sqrt(10.^(diffSPL./10)) in Torben's script WHY sqrt?
att = toneAtt.*attFactor;%this is the value for multiplying signal scaled/clipped to [-1 to 1]
noise(1,:)=noise(1,:)*att(1);
noise(2,:)=noise(2,:)*att(2);
%should the two speakers dB be added?

%% generate signal
if EmbedSignal
    %make signal
    t=linspace(0,SignalDuration,[SamplingRate*SignalDuration]); %time vector for chirp
    switch SignalForm
        case 'LinearUpsweep'
            signal=chirp(t,SignalMinFreq,SignalDuration,SignalMaxFreq);
            freqvec=SignalMinFreq+(SignalMaxFreq-SignalMinFreq)*t;
        case 'LinearDownsweep' %gaussian noise from mean 0 std .25
            signal=chirp(t,SignalMaxFreq,SignalDuration,SignalMinFreq);
            freqvec=SignalMaxFreq+(SignalMinFreq-SignalMaxFreq)*t;
        case 'QuadraticConvex'
            tnew=t-mean(t);
            signal=chirp(tnew,SignalMinFreq,SignalDuration./2,SignalMaxFreq,'quadratic',[],'convex'); %make chirp
            freqvec=SignalMinFreq+(SignalMaxFreq-SignalMinFreq)./tnew(1)*tnew.^2;
    end
    
    %adjust signal volume
    SoundCal = BpodSystem.CalibrationTables.SoundCal;
    if(isempty(SoundCal))
        disp('Error: no sound calibration file specified');
        return
    end
    toneAtt = polyval(SoundCal(1,1).Coefficient,freqvec);%Frequency dependent attenuation factor with less attenuation for higher frequency (based on calibration polynomial)
    %toneAtt = [polyval(SoundCal(1,1).Coefficient,toneFreq)' polyval(SoundCal(1,2).Coefficient,toneFreq)']; in Torben's script
    diffSPL = SignalVolume - [SoundCal.TargetSPL];
    attFactor = sqrt(10.^(diffSPL./10)); %sqrt(10.^(diffSPL./10)) in Torben's script WHY sqrt?
    att = toneAtt.*attFactor;%this is the value for multiplying signal scaled/clipped to [-1 to 1]
    signal=signal.*att;
end


%% embed signal in noise at t with t being drawn from an exponential distribution with a mean such that signal is centered in noise (if you are not me and reading this: sorry for this umstaendlich Beschreibung; if you are not me and not German and reading this: sorry for writing in Germlish)
if EmbedSignal
    maxLat=NoiseDuration-SignalDuration;%maximum latency
    if maxLat<0
        warning('Cannot place signal in noise. Check noise and signal duration. Will produce pure noise stimulus...')
        stimulus=noise;
        embed=nan;
    else
        SignalLatency=(NoiseDuration-SignalDuration)./2;
        lat=min( [exprnd(SignalLatency) maxLat] );        
        stimulus=zeros(size(noise));
        SignalStartSample=floor(lat*SamplingRate)+1;
        SignalEndSample=SignalStartSample+length(signal)-1;
        stimulus(SignalStartSample:SignalEndSample)=signal;
        stimulus=stimulus+noise;
        embed=lat;
    end
else
    stimulus=noise;
    embed=nan;
end


%put an envelope to avoide clicking sounds at beginning and end
omega=(acos(sqrt(0.1))-acos(sqrt(0.9)))/(Ramp/pi*2); % This is for the envelope with Ramp duration duration
t=0 : (1/SamplingRate) : pi/2/omega;
t=t(1:(end-1));
RaiseVec= (cos(omega*t)).^2;

Envelope = ones(length(stimulus),1); % This is the envelope
Envelope(1:length(RaiseVec)) = fliplr(RaiseVec);
Envelope(end-length(RaiseVec)+1:end) = (RaiseVec);

stimulus = stimulus.*Envelope';
end

    function x = f_alpha_gaussian ( n, q_d, alpha )
        
        %*****************************************************************************80
        %
        %% F_ALPHA_GAUSSIAN generates noise using a Gaussian distribuion.
        %
        %  Discussion:
        %
        %    This function generates a discrete colored noise vector X of size N
        %    with a power spectrum distribution that behaves like 1 / f^ALPHA.
        %    The white noise used in the generation is sampled from a Gaussian
        %    distribution with zero mean and a variance of Q_D.
        %
        %  Licensing:
        %
        %    This code is distributed under the GNU LGPL license.
        %
        %  Modified:
        %
        %    30 March 2011
        %
        %  Author:
        %
        %    Miroslav Stoyanov
        %
        %  Reference:
        %
        %    Jeremy Kasdin,
        %    Discrete Simulation of Colored Noise and Stochastic Processes
        %    and 1/f^a Power Law Noise Generation,
        %    Proceedings of the IEEE,
        %    Volume 83, Number 5, 1995, pages 802-827.
        %
        %    Miroslav Stoyanov, Max Gunzburger, John Burkardt,
        %    Pink Noise, 1/F^Alpha Noise, and Their Effect on Solution
        %    of Differential Equations,
        %    submitted, International Journal for Uncertainty Quantification.
        %
        %  Parameters:
        %
        %    Input, integer N, the number of elements in the discrete noise vector.
        %
        %    Input, real Q_D, the variance of the Gaussian distribution.  A standard
        %    Gaussian distribution has a variance of 1.  The variance is the square
        %    of the standard deviation.
        %
        %    Input, real ALPHA, specifies that the computed noise is to have a
        %    power spectrum that behaves like 1/f^alpha.  Normally 0 <= ALPHA <= 2.
        %
        %    Output, real X(N), the computed discrete noise vector.
        %
        
        %
        %  Set the standard deviation of the white noise.
        %
        stdev = sqrt ( abs ( q_d ) );
        %
        %  Generate the coefficients.
        %
        hfa = zeros ( 2 * n, 1 );
        hfa(1) = 1.0;
        for i = 2 : n
            hfa(i) = hfa(i-1) * ( 0.5 * alpha + ( i - 2 ) ) / ( i - 1 );
        end
        hfa(n+1:2*n) = 0.0;
        %
        %  Sample the white noise.
        %
        wfa = stdev * randn ( n, 1 );
        %
        %  Pad the array with zeros in anticipation of the FFT process.
        %
        z = zeros ( n, 1 );
        wfa = [ wfa; z ];
        %
        %  Perform the discrete Fourier transforms.
        %
        fh = fft ( hfa );
        fw = fft ( wfa );
        %
        %  Multiply the two complex vectors.
        %
        fh = fh(1:n+1);
        fw = fw(1:n+1);
        
        fw = fh .* fw;
        %
        %  Scale to match the conventions of the Numerical Recipes FFT code.
        %
        fw(1)   = fw(1)   / 2.0;
        fw(end) = fw(end) / 2.0;
        %
        %  Pad the array with zeros in anticipation of the FFT process.
        %
        z = zeros ( n - 1, 1 );
        fw = [ fw; z ];
        %
        %  Take the inverse Fourier transform.
        %
        x = ifft ( fw );
        %
        %  Only the first half of the inverse Fourier transform is useful.
        %
        x = 2.0 * real ( x(1:n) );
        
        return
    end