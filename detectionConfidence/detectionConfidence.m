function detectionConfidence()
% Bpod script for an auditory detection in noise task with
% time-investment-based confidence reports
% Katharina Schmack, Cold FpharmSpring Harbr Laboratory, March 2018

global BpodSystem
global TaskParameters

%% Resolve HiFi Module USB port
if (isfield(BpodSystem.ModuleUSB, 'HiFi1'))
    %% Create an instance of the HiFi module
    H = BpodHiFi(BpodSystem.ModuleUSB.HiFi1);
else
    error('Error: To run this protocol, you must first pair the HiFi module with its USB port. Click the USB config button on the Bpod console.')
end

%% prepare TaskParameters
TaskParameters = BpodSystem.ProtocolSettings;%get saved TaskParameters
TaskParameters=syncTaskParameters(TaskParameters);%sync saved TaskParameters with DefaultTaskParameters
BpodParameterGUI_KEPECS('init', TaskParameters);
BpodSystem.Status.Pause = 1;
HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
TaskParameters = BpodParameterGUI_KEPECS('sync', TaskParameters); % Sync parameters with BpodParameterGUI plugin
BpodSystem.ProtocolSettings = TaskParameters; % copy settings back prior to saving
SaveBpodProtocolSettings;

%% Setup HiFi Module
H.SynthWaveform = 'WhiteNoise';
H.SynthAmplitudeFade = 1920; % White noise fades in and out over 1920 samples, 0.1s at 192khz (it then remains on for entire experiment)
H.DigitalAttenuation_dB = TaskParameters.GUI.AudioAttenuation;

%% define some pathes
[~,BpodSystem.Data.Custom.Rig] = system('hostname');
BpodSystem.Data.Custom.Subject = BpodSystem.Status.CurrentSubjectName;
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler';
[filepath,filename,~]=fileparts(BpodSystem.Path.CurrentDataFile);
BpodSystem.Data.Custom.StimulusPath=fullfile(strrep(filepath,'Session Data','Session Stimuli'),filename);
if ~exist(BpodSystem.Data.Custom.StimulusPath,'dir')
    mkdir(BpodSystem.Data.Custom.StimulusPath)
end
BpodSystem.Data.Custom.recordBaselineTrial=false;


%% User Input for Photometry & Pharmacoloy & Blocks (to do: add option in GUI to switch on and off)
if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode
    site = questdlg('Where are you recording from?', ...
        'photometry site', ...
        'left','right','right');
    BpodSystem.Data.Custom.PhotometrySite=site;
end
if TaskParameters.GUI.PharmacologyOn
    drug = questdlg('What did you administer?', ...
        'drug condition', ...
        'haloperidol015','vehicle','haloperidol020');
    BpodSystem.Data.Custom.Pharmacology=drug;
end
if TaskParameters.GUI.LaserPercentage>0
    PulsePal('COM3');
    load('ParameterMatrix25Hz5ms200ms.mat','ParameterMatrix');
    ParameterMatrix(2:13,4)={0 5 0 0.004 0 0 1/2-0.004 TaskParameters.GUI.LaserBlockLength*30 0 TaskParameters.GUI.LaserBlockLength*30 0 1};
    ParameterMatrix(2:13,5)={0 5 0 0.004 0 0 1/2-0.004 TaskParameters.GUI.LaserBlockLength*30 0 TaskParameters.GUI.LaserBlockLength*30 0 1};
    ParameterMatrix{2,8}=1;
    ProgramPulsePal(ParameterMatrix);
        site = questdlg('Where are you stimulating?', ...
        'stimulation site', ...
        'left','left10','right','right');
    BpodSystem.Data.Custom.StimulationSite=site;
end

if TaskParameters.GUI.DetermineFirstBlock
    firstblock = questdlg('Which block do you want to start with?', ...
        'first block', ...
        '30','50','70','50');
    BpodSystem.Data.Custom.firstblock=firstblock;
end
%% Initialize photometry if PhotometryOn
if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode
    TaskParameters.nidaq.duration = 120;
    TaskParameters.nidaq.IsContinuous = true;
    TaskParameters.nidaq.updateInterval = 0.1; % save new data every n seconds
    BpodSystem.PluginObjects.Photometry.baselinePeriod = [0 1]; % kludge, FS
    BpodSystem.ProtocolSettings = TaskParameters; % copy settings back because syncPhotometrySettings relies upon BpodSystem.ProtocolSettings
    TaskParameters = initPhotometry(TaskParameters);
    
    %% initialize photometry plot
    updatePhotometryPlotKatharina('init',[0 0],{'Reward','Stimulus'});
    
    %% alternate LED modulation mode
    if TaskParameters.GUI.PhotometryOn==2
        % store initial LED settings
        storedLED1_amp = TaskParameters.GUI.LED1_amp;
        storedLED2_amp = TaskParameters.GUI.LED2_amp;
    end
end

%% Initialize trial specific values and save them in BpodSystem.Data.Custom
updateCustomDataFields(0, H);
BpodSystem.Data.Custom = orderfields(BpodSystem.Data.Custom);

%% Initialize main plot
BpodSystem.ProtocolFigures.SideOutcomePlotFig = figure('Position', [0, 600, 1000, 400],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleOutcome = axes('Position',    [  .05  0.08 0.9 .45]);
BpodSystem.GUIHandles.OutcomePlot.HandlePsycAud = axes('Position',    [2*.05 + 1*.08   .62  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleTrialRate = axes('Position',  [3*.05 + 2*.08   .62  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleFix = axes('Position',        [4*.05 + 3*.08   .62  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleST = axes('Position',         [5*.05 + 4*.08   .62  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleFeedback = axes('Position',   [6*.05 + 5*.08   .62  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleVevaiometric = axes('Position',   [7*.05 + 6*.08   .62  .1  .3], 'Visible', 'off');
MainPlot(BpodSystem.GUIHandles.OutcomePlot,'init');


%% Task loop over trials
RunSession = true;
iTrial = 1;
while RunSession
    
    % update TaskParameters
    if TaskParameters.GUI.PhotometryOn==2
        LEDmode = rem(iTrial, 3);
        switch LEDmode
            case 1
                TaskParameters.GUI.LED1_amp = storedLED1_amp;
                TaskParameters.GUI.LED2_amp = storedLED2_amp;
            case 2
                TaskParameters.GUI.LED1_amp = storedLED1_amp;
                TaskParameters.GUI.LED2_amp = 0;
            case 0
                TaskParameters.GUI.LED1_amp = 0;
                TaskParameters.GUI.LED2_amp = storedLED2_amp;
        end
    end
    TaskParameters = BpodParameterGUI_KEPECS('sync', TaskParameters);
    BpodSystem.ProtocolSettings = TaskParameters; % copy settings back prior to saving
    SaveBpodProtocolSettings;
    
    %% prepare stateMatrix
    %record baseline on first trial 
%     if TaskParameters.GUI.PhotometryOn~=0&&TaskParameters.GUI.BaselineRecording>0&&iTrial==1
%         iTrial=recordBaselineTrial(iTrial);
%     end
    sma = stateMatrix(iTrial);
    SendStateMatrix(sma);
    
    %% prepare photometry
    if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode
        preparePhotometryAcq(TaskParameters);
    end
    
    %% RUN!!!
    RawEvents = RunStateMatrix();
    
    %% stop photometry session
    if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode
        stopPhotometryAcq;
    end
    
    %% process photometry session
    if ~isempty(fieldnames(RawEvents))
        if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode
            processPhotometryAcq(iTrial);
        end
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents);
        SaveBpodSessionData;
    else
        disp([' *** Trial # ' num2str(iTrial) ':  aborted, data not saved ***']); % happens when you abort early (I think), e.g. when you are halting session
    end
    
    %% check for protocol manually ended by user
    if BpodSystem.Status.BeingUsed == 0
        H.SynthAmplitude = 0;
        return
    end
    
    %% analyze behavior and create new trial-specific parameters for next trial
    updateCustomDataFields(iTrial, H)%
    
    %% update main plot
    MainPlot(BpodSystem.GUIHandles.OutcomePlot,'update',iTrial);
    
    %% update photometry plot
    if TaskParameters.GUI.PhotometryOn && ~BpodSystem.EmulatorMode 
        processPhotometryOnline(iTrial);
        if BpodSystem.Data.Custom.RewardReceivedCorrect(iTrial)>0
            rewtime=BpodSystem.Data.Custom.RewardStartTime(iTrial);
        else
            rewtime=nan;
        end
        if (BpodSystem.Data.Custom.CoutEarly(iTrial))~=1
            stimtime=BpodSystem.Data.Custom.StimulusStartTime(iTrial);
        else
            stimtime=nan;
        end
        try %kludge, problem with plotting baseline
        updatePhotometryPlotKatharina('update', [rewtime stimtime],{'reward','stimulus'});
        catch
            fprintf('Trial %d: Error in updatePhotometryPlotKatharina.m\n',iTrial)
        end
    end
    
    %% Go on to next trial
    iTrial = iTrial + 1;
%     %record baseline on 
%     if BpodSystem.Pause==1
%         iTrial=recordBaselineTrial(iTrial);
%     end
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
end %while loop for trials

end



