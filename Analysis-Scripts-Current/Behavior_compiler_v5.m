%% script to combine Bpod HALIP data across multiple sessions with the same treatment
% 
path = {
%    
% '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\Parker_Lab\JustinAnair\Data\HALIP\DATCre_TRPV1KO_AAVTRPV1\Bilateral Expression\Veh_Veh'
'\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\Parker_Lab\JustinAnair\Data\HALIP\DATCre_TRPV1KO_AAVTRPV1\Bilateral Expression\Veh_Cap'
% % '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\Parker_Lab\JustinAnair\Data\HALIP\B6'
% '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\Parker_Lab\JustinAnair\Data\HALIP\TRPV1Ko'
};
warning('off', 'MATLAB:MKDIR:DirectoryExists');

for x = 1:length(path)
    
    cd(path{x}); 
    fileNames = dir;
    fileNames = struct2cell(fileNames); 
    
    fIndices = find(contains(fileNames(1,:), 'f')); 
    mIndices = find(contains(fileNames(1,:), 'm'));
    mouseIndices = [fIndices mIndices];
    
    
    for y = 1:length(mouseIndices); 
        
        cd(path{x}); 
        cd(fileNames{1, mouseIndices(y)}); 
        fileNames2 = dir;
        fileNames2 = struct2cell(fileNames2); 
        fileIndices = find(contains(fileNames2(1,:), 'detectionConfidence'));
    
    
%options
easyTrials = 26; %the first trial we want to consider in terms of performance
FAcutoff = 1.5; %for some reasons these there are about 0.5 s shorter than what they seeem. Schmack paper uses 2s cutoff
% which corresongs to (1.5s) here. We can play with this number to include/sxclude more or less trials. 
bckg = 40; %dB for subtracting form Signal Volume for signal-to-noise differnece for psychometric curve (note that no signal trials = -40dB)
pctCorrectCutoff = .65; %fraction of high signal-to-noise trials that must be correct to retain the session
SNRthresh = 55; %sound intensity threshold for calculating pctCorrectCutoff (i.e., high SNR trials) 
binsize = 3;
numberbinsper15 = 15/binsize;

%pre-allocation

%data compiled for this mouse
CompiledData.nTrials = 0; 
CompiledData.RewardReceivedTotal = []; 
CompiledData.ResponseCorrect = []; 
CompiledData.ResponseTime = []; 
CompiledData.EmbedSignal = []; 
CompiledData.RewardStartTime = []; 
CompiledData.RewardReceivedCorrect = []; 
CompiledData.CatchTrial = []; 
CompiledData.Trial = {};%the raw event data for every trial compiled 
CompiledData.Ports = 0;
CompiledData.psychData = []; 
CompiledData.SignalVolume = [];
CompiledData.Exclusions.rawTrialTotal = [];
CompiledData.Exclusions.trialsExcluded =[];
CompiledData.Exclusions.nSessions = [];
CompiledData.Exclusions.nSessionsClean = [];
CompiledData.Binned.BinIndices = [];
CompiledData.Binned.sessionFArate = [];
CompiledData.Binned.FAcount = [];
CompiledData.Binned.concatbinFA = [];
CompiledData.Binned.concatbinInvest = [];
CompiledData.Binned.sessionFAinvest = [];

%summary data for compiling across mice (except for the psychData, which must be concatenated across the trials of all mice within a treatment condition. 
OutputData.FalseAlarmRate = 0;
OutputData.FalseAlarmInvest = 0;
OutputData.PercentCorrect = 0;
OutputData.TotalWater = 0;
OutputData.CorrectSignalTrial = 0;
OutputData.MeanCatchTrialInvest = 0;
OutputData.psychData = []; %left column are SPL intensities (minus background 40dB). right column are 1-0's denoting self-report of having heard the stimulus (1).
OutputData.Exclusions.rawTrialTotal = [];
OutputData.Exclusions.nSessionsClean = [];
OutputData.normfalseAlarmInvestments = [];
OutputData.timelineFA = [];
OutputData.Binned.CorrectSignal = [];
OutputData.Binned.FArate = [];
OutputData.Binned.FAcount = [];
OutputData.Binned.concatbinFA = [];
OutputData.Binned.concatbinInvest = [];
OutputData.Binned.FAinvest = [];

%metrics for binning FA rate and bin counts
sessionFABinrate = [];
FABinned = [];
sessionFABincount = [];
FABin_counts = [];
sessionFABin_invest = [];
FABin_invest =[];


for session_idx = 1:length(fileIndices) %calculating metrics for all sessions for one mouse
     
    load(fileNames2{1, fileIndices(session_idx)}); 
    
    
    theseWaitTimes = SessionData.Custom.WaitingTime(easyTrials:end); 
    trials2keep = find(SessionData.Custom.WaitingTime >= FAcutoff);
    trials2keep(trials2keep < easyTrials) = []; 
    
    
    theseSignalVolumes = SessionData.Custom.SignalVolume(trials2keep); 
    theseResponseCorrect = SessionData.Custom.ResponseCorrect(trials2keep); 
    SignalCorrectThreshold = find(theseSignalVolumes >= SNRthresh & SessionData.Custom.EmbedSignal(trials2keep) == 1);
    
    
    if nanmean(theseResponseCorrect(SignalCorrectThreshold)) >= pctCorrectCutoff %cut out sessions less than this performance on signal trials
    CompiledData.Exclusions.nSessionsClean = [sum(CompiledData.Exclusions.nSessionsClean) + 1];
    
    CompiledData.RewardReceivedTotal = [CompiledData.RewardReceivedTotal SessionData.Custom.RewardReceivedTotal]; 
    
    CompiledData.ResponseCorrect = [CompiledData.ResponseCorrect SessionData.Custom.ResponseCorrect(trials2keep)];
    
    CompiledData.ResponseTime = [CompiledData.ResponseTime SessionData.Custom.ResponseTime(trials2keep)];
    
    CompiledData.EmbedSignal = [CompiledData.EmbedSignal SessionData.Custom.EmbedSignal(trials2keep)];
    
    CompiledData.RewardStartTime = [CompiledData.RewardStartTime SessionData.Custom.RewardStartTime(trials2keep)];
    
    CompiledData.RewardReceivedCorrect = [CompiledData.RewardReceivedCorrect SessionData.Custom.RewardReceivedCorrect(trials2keep)]; 
    
    CompiledData.CatchTrial = [CompiledData.CatchTrial SessionData.Custom.CatchTrial(trials2keep)]; 
    
    CompiledData.Trial = [CompiledData.Trial SessionData.RawEvents.Trial(trials2keep)];
    
    CompiledData.Exclusions.rawTrialTotal = [CompiledData.Exclusions.rawTrialTotal numel(SessionData.Custom.WaitingTime)];
    
    CompiledData.Ports(1) = SessionData.SettingsFile.GUI.Ports_LMR;
       
    CompiledData.SignalVolume = [CompiledData.SignalVolume SessionData.Custom.SignalVolume(trials2keep)];

    Bin_Indices = [];
    Bin_Indices = binFA(SessionData, binsize, easyTrials); %adjust bin size by changing value after SessionData
   
    %% pullout useful information to evaluate in the for loop for graph line
        currentSession = SessionData;
    % find incorrect trials to find FA trials
        SessionData4Bin = currentSession.RawEvents.Trial(easyTrials:length(currentSession.RawEvents.Trial));
        incorrectresponse = find(currentSession.Custom.ResponseCorrect == 0);
        incorrectresponse_keep = intersect(incorrectresponse,trials2keep); 

        noSignalTrial = find(currentSession.Custom.EmbedSignal == 0);
        FAresponse = intersect(incorrectresponse_keep, noSignalTrial) - (easyTrials - 1); %find where these overlap and subtract easy trials to start at 1 (matches bin_indices better)
    % find correct trials to measure correct performance
        correctresponse = find(currentSession.Custom.ResponseCorrect == 1);
        correctresponse_keep = intersect(correctresponse,trials2keep);

        SignalTrial = find(currentSession.Custom.EmbedSignal == 1);
        correctsignal = intersect(correctresponse_keep, SignalTrial) - (easyTrials - 1);

%compute the percentage of signal trials that were correct or FA per bin after injection
 if Bin_Indices(1,1) ~= 0
    for idx2 = 1:length(Bin_Indices)
            currentTrials = Bin_Indices(1,idx2):Bin_Indices(2,idx2); %determine trials for first block to evaluate
            if max(currentTrials) - min(currentTrials) <=3 %criteria for excluding bins with little or no trials
                FABin_counts = horzcat(FABin_counts, NaN);
                FABinned = horzcat(FABinned, NaN); %store each FA rate per block in this array
                FABin_invest = horzcat(FABin_invest, NaN);
            else 
                currentFAtrial = intersect(currentTrials,FAresponse); %find where FA trials intersect within the first block
                if length(currentFAtrial) > 0
                for i3 = 1:length(currentFAtrial)
                     if CompiledData.Ports == 321;
                         if isfield(SessionData4Bin{1, currentFAtrial(1,i3)}.Events, 'Port3Out') 
                             FA_invest= SessionData4Bin{1, currentFAtrial(1,i3)}.Events.Port3Out(end) - SessionData4Bin{1, currentFAtrial(1,i3)}.Events.GlobalTimer1_Start;
                         else     
                             FA_invest = nan; 
                         end
                    end

                     if CompiledData.Ports == 123;
                         if isfield(SessionData4Bin{1, currentFAtrial(1,i3)}.Events, 'Port1Out') 
                             FA_invest = SessionData4Bin{1, currentFAtrial(1,i3)}.Events.Port1Out(end) - SessionData4Bin{1, currentFAtrial(1,i3)}.Events.GlobalTimer1_Start;
                         else     
                             FA_invest = nan; 
                         end
                     end
                     Indivinvest(1,i3) = FA_invest;
                end
                else
                    FA_invest = 0;
                    Indivinvest = FA_invest;
                end
                FABin_counts = horzcat(FABin_counts, sum(length(currentFAtrial)));
                rate_onebin = (sum(currentFAtrial) / sum(currentTrials))*100; %determine FA rate for the entire block
                FABinned = horzcat(FABinned, rate_onebin); %store each FA rate per block in this array
                FABin_invest = horzcat(FABin_invest, mean(Indivinvest));
            end 
    end 
 end 
        sessionFABincount = [sessionFABincount; FABin_counts];
        sessionFABinrate = [sessionFABinrate; FABinned];
        sessionFABin_invest = [sessionFABin_invest; FABin_invest];
        FABinned = []; %reset FAbinned for next session
        FABin_invest = [];
        FABin_counts = [];
        CompiledData.Binned.BinIndices{session_idx, 1} = Bin_Indices; %store for troubleshooting
    end

end 
sessionFABin_invest(sessionFABin_invest == 0) = NaN; % omit zeros to avoid weighing the investment time score down.
CompiledData.Binned.sessionFAinvest = sessionFABin_invest;
CompiledData.Binned.FAcount = sessionFABincount;
CompiledData.Binned.concatbinFA = sessionFABin_invest;
CompiledData.Binned.sessionFArate = sessionFABinrate;


OutputData.Binned.FAcount = [OutputData.Binned.FAcount, sum(sessionFABincount)];
OutputData.Binned.FArate = [OutputData.Binned.FArate; mean(sessionFABinrate,'omitmissing')];
OutputData.Binned.FAinvest = [OutputData.Binned.FAinvest; mean(sessionFABin_invest, 'omitmissing')];

for i2 = 1:numberbinsper15
    FABin_Summed(1, i2) = mean(mean(sessionFABinrate(:,i2:numberbinsper15:end), 'omitmissing'));
    FABin_investsummed(1,i2) = mean(mean(sessionFABin_invest(:,i2:numberbinsper15:end), 'omitmissing'));
end 

CompiledData.Binned.concatbinFA = FABin_Summed;
OutputData.Binned.concatbinFA = FABin_Summed;
OutputData.Binned.concatbinInvest = FABin_investsummed;
%


% calculate trials and sessions excluded

CompiledData.nSessions = numel(fileIndices);

OutputData.Exclusions.rawTrialTotal = CompiledData.Exclusions.rawTrialTotal;
OutputData.Exclusions.nSessionsClean = CompiledData.Exclusions.nSessionsClean;
OutputData.nSessions = CompiledData.Exclusions.nSessions;

%basic data for output
OutputData.PercentCorrect = nanmean(CompiledData.ResponseCorrect)*100; %determines correct trials from struct
OutputData.TotalWater = nansum(CompiledData.RewardReceivedTotal); 

        
%identify oddball trials (early withdrawal from center port)
oddErrors = find(isnan(CompiledData.ResponseCorrect)); 

%identify oddball trials (early withdrawal from correct reward port)
oddErrors_2 = find(isnan(CompiledData.RewardStartTime)); 

%pull these data out of the structre
noSignalTrials = find(CompiledData.EmbedSignal == 0);
signalTrials = find(CompiledData.EmbedSignal == 1); 

% correctTrials = find(CompiledData.RewardReceivedCorrect);
% incorrectTrials = find(CompiledData.RewardReceivedCorrect == 0);

correctTrials = find(CompiledData.ResponseCorrect);
incorrectTrials = find(CompiledData.ResponseCorrect == 0);

catchTrials = find(CompiledData.CatchTrial); 
signalVolumes = CompiledData.SignalVolume; 

%remove indices of oddball trials (early withdrawal from center port) from all

noSignalTrials = setdiff(noSignalTrials, oddErrors);
signalTrials = setdiff(signalTrials, oddErrors);

correctTrials = setdiff(correctTrials, oddErrors);
incorrectTrials = setdiff(incorrectTrials, oddErrors); 
catchTrials = setdiff(catchTrials, oddErrors); %ignore catch trials where mouse did not stay in center

signalVolumes(1, noSignalTrials) = 0; %set the pre-determined SPL for no signal trials to zero; 
signalVolumes = signalVolumes - bckg; %subtract the background SPL

%remove indices of oddball trials (early withdrawal from correct reward port) from
%incorrect trials and ignore catchTrials that were incorrect
incorrectTrials = setdiff(incorrectTrials, oddErrors_2);
catchTrials = setdiff(catchTrials, incorrectTrials); %ignore catch trials when incorrect


%compute the percentage of signal trials that were correct overall

correctSignalTrials = intersect(correctTrials, signalTrials); 
missTrials = intersect(incorrectTrials, signalTrials); %incorrect signal trials
falseAlarms = intersect(noSignalTrials, incorrectTrials); 

OutputData.CorrectSignalTrial = nnz(correctSignalTrials)/nnz(signalTrials);





%% compute missTrial time investments

if isempty(missTrials)
    OutputData.missTrialRate = nan;
    OutputData.missTrialInvest = nan;  
    medianMissTrialInvest = nan; 
else
    
OutputData.missTrialRate = length(missTrials)/length(signalTrials)*100; 



%now compute the false alarm time investments
missTrialInvestments = nan(1, length(missTrials)); 

for idx2 = 1:length(missTrials)

    if CompiledData.Ports == 321;
        
        if isfield(CompiledData.Trial{1, missTrials(idx2)}.Events, 'Port1Out') 
       
        missTrialInvestments(idx2) = CompiledData.Trial{1, missTrials(idx2)}.Events.Port1Out(end) - CompiledData.Trial{1, missTrials(idx2)}.Events.GlobalTimer1_Start;
    
        else
            
            missTrialInvestments(idx2) = nan; 
            
        end
               
        
        
   end
    
    
    
    if SessionData.SettingsFile.GUI.Ports_LMR == 123;
        
        if isfield(CompiledData.Trial{1, missTrials(idx2)}.Events, 'Port3Out') 
    
        missTrialInvestments(1,idx2) = CompiledData.Trial{1, missTrials(idx2)}.Events.Port3Out(end) - CompiledData.Trial{1, missTrials(idx2)}.Events.GlobalTimer1_Start;
        
                else
            
            missTrialInvestments(1,idx2) = nan; 
            
        end
               
        
    end
end


medianMissTrialInvest = nanmedian(missTrialInvestments);
% falseAlarmInvestments(falseAlarmInvestments < FAcutoff) = []; %cutoff (2s) below which we don't consider them to be investments. 

OutputData.missTrialInvest = nanmean(missTrialInvestments); 


end


%% compute the false alarm rate and time investmtnets


if isempty(falseAlarms)
    OutputData.FalseAlarmRate = nan;
    OutputData.FalseAlarmInvest = nan;  
    pctHiFAinvestment = nan; 
else
    
OutputData.FalseAlarmRate = length(falseAlarms)/length(noSignalTrials)*100; 



%now compute the false alarm time investments
falseAlarmInvestments = nan(1, length(falseAlarms)); 

for idx2 = 1:length(falseAlarms)

    if CompiledData.Ports == 321;
        
        if isfield(CompiledData.Trial{1, falseAlarms(idx2)}.Events, 'Port3Out') 
       
        falseAlarmInvestments(idx2) = CompiledData.Trial{1, falseAlarms(idx2)}.Events.Port3Out(end) - CompiledData.Trial{1, falseAlarms(idx2)}.Events.GlobalTimer1_Start;
    
        else
            
        falseAlarmInvestments(idx2) = nan; 
            
        end
               
        
        
   end
    
    
    
    if SessionData.SettingsFile.GUI.Ports_LMR == 123;
        
        if isfield(CompiledData.Trial{1, falseAlarms(idx2)}.Events, 'Port1Out') 
    
        falseAlarmInvestments(idx2) = CompiledData.Trial{1, falseAlarms(idx2)}.Events.Port1Out(end) - CompiledData.Trial{1, falseAlarms(idx2)}.Events.GlobalTimer1_Start;
        
                else
            
         falseAlarmInvestments(idx2) = nan; 
            
        end
               
        
    end
end

% falseAlarmInvestments(falseAlarmInvestments < FAcutoff) = []; %cutoff (2s) below which we don't consider them to be investments. 

OutputData.FalseAlarmInvest = nanmean(falseAlarmInvestments); 
OutputData.normfalseAlarmInvestments = zscore(falseAlarmInvestments);
OutputData.pctHiFAinvestment = nnz(falseAlarmInvestments > medianMissTrialInvest)/length(noSignalTrials)*100; 

end




%% compute the catchtrial time investments

% there are both signal and no signal catch trials. 
noSignalCatchTrials = intersect(catchTrials, noSignalTrials);
signalCatchTrials = intersect(catchTrials, correctSignalTrials);


if isempty(noSignalCatchTrials)
    
    noSignalCatchTrialInvestments = []; 
    
  else
    
    noSignalCatchTrials(noSignalCatchTrials > length(CompiledData.Trial(1, :))) = []; %get rid of the last one, since the chamber stops recording so we don't get the timestamp. 
    
end

if ~isempty(noSignalCatchTrials)  
    
    for session_idx = 1:length(noSignalCatchTrials)

    if CompiledData.Ports == 321;
        
        if isfield(CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events, 'Port1Out') 
        
        noSignalCatchTrialInvestments(session_idx) = CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events.Port1Out(end) - CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events.GlobalTimer1_Start;
        
        else
            
        noSignalCatchTrialInvestments(session_idx) = nan; 
        
        end
        
    end
    
    if CompiledData.Ports == 123;
        
        if isfield(CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events, 'Port3Out')

        
        noSignalCatchTrialInvestments(session_idx) = CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events.Port3Out(end) - CompiledData.Trial{1, noSignalCatchTrials(session_idx)}.Events.GlobalTimer1_Start;
        
        else
            
        noSignalCatchTrialInvestments(session_idx) = nan; 
        
        end
    
    end
    end
end


if isempty(signalCatchTrials)
    
    signalCatchTrialInvestments = []; 
    
  else
    
    signalCatchTrials(signalCatchTrials > length(CompiledData.Trial(1, :))) = []; %get rid of the last one, since the chamber stops recording so we don't get the timestamp. 
    
end

if ~isempty(signalCatchTrials)  
    
    for session_idx = 1:length(signalCatchTrials)

    if CompiledData.Ports == 321;
        
        if isfield(CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events, 'Port3Out') 
        
            signalCatchTrialInvestments(session_idx) = CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events.Port3Out(end) - CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events.GlobalTimer1_Start;
        
        else
            
        signalCatchTrialInvestments(session_idx) = nan; 
        
        end
        
    end
    
    if CompiledData.Ports == 123;
        
        if isfield(CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events, 'Port1Out')
        
        signalCatchTrialInvestments(session_idx) = CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events.Port1Out(end) - CompiledData.Trial{1, signalCatchTrials(session_idx)}.Events.GlobalTimer1_Start;
        
        else
            
        signalCatchTrialInvestments(session_idx) = nan; 
        
        end
        
    
    end
    end
end
    
    catchTrialInvestments = [noSignalCatchTrialInvestments signalCatchTrialInvestments]; 
    

    if isempty(catchTrials)  
    
    OutputData.MeanCatchTrialInvest = nan; 
    
    else
        
    OutputData.MeanCatchTrialInvest = nanmean(catchTrialInvestments); 
    
    
    end
    
    
        if isempty(noSignalCatchTrials)
        
        OutputData.MeanNoSignalCatchTrialInvestment = nan; 
        
        else
            
        OutputData.MeanNoSignalCatchTrialInvestment = nanmean(noSignalCatchTrialInvestments); 
        
        end
        
        
        if isempty(noSignalCatchTrials)
            
        OutputData.MeanSignalCatchTrialInvestment = nan; 
        
        else
            
        OutputData.MeanSignalCatchTrialInvestment =nanmean(signalCatchTrialInvestments); 
        
        end
        
            
            

%% generate vector for psychometric curve fitting


choseSignalNoSignal = intersect(incorrectTrials, noSignalTrials); 
choseSignalTrials = sort([correctSignalTrials choseSignalNoSignal]); 

choiceTrace = zeros(1,length(signalVolumes)); 
choiceTrace(1, choseSignalTrials) = 1; 

clean_coiceTrace = choiceTrace; 
clean_coiceTrace(oddErrors) = []; %remove early exit trials

clean_signalVolumes = signalVolumes; 
clean_signalVolumes(oddErrors) = []; %remove early exit trials


CompiledData.psychData = nan(length(clean_signalVolumes), 2); 

CompiledData.psychData(:, 1) = clean_signalVolumes'; 
CompiledData.psychData(:, 2) = clean_coiceTrace'; 

OutputData.psychData = CompiledData.psychData; 

 
    
    %export compiled data for individual mouse
    
    mkdir('Analysis'); 
    
    cd('Analysis')
    
    save('CompiledData.mat', 'CompiledData'); 
    save('OutputData.mat', 'OutputData'); 
     end
    
end

disp('Behavior_Compiler_v5 analysis is complete')
clear all;
   
