%% get folder
myFolder = uigetdir();

%% Create a list of all files in the folder
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for Idx = 1:length(theFiles)
    baseFileName = theFiles(Idx).name;
    str = strcat(fullfile(theFiles(Idx).folder, baseFileName));
    myFilenames{Idx} = str;
end
%% Get Dates
Dates = extractBetween(myFilenames,113,119); %will change based on path
%% Looping Analysis (Extract Water volume, psychometric data, average response time, anything I want to look at really)
DataTable = cell(Idx,9); %if you want to track more variables, add more columns (change "5")
for Idx2 = 1:length(myFilenames)
    BehaviorSource = char(myFilenames(Idx2));
    load(BehaviorSource);
%% Get the session dates
    GetSessionDate = SessionData.Info.SessionDate; %loads session date from struct
        DateStr = GetSessionDate;
        formatIn = "dd-mmm-yyyy"
        dt = datenum(DateStr,formatIn)
        if Idx2 <= 1
            dt1 = dt;
            timediff = dt1 - dt
        else
            timediff = dt - dt1
        end
        %DataTable(Idx2,1) = cellstr(GetSessionDate);
        DataTable(Idx2, 1) = num2cell(timediff);
%% determine some basic metrics
    TotalTrials = SessionData.nTrials;
        DataTable(Idx2,2) = num2cell(TotalTrials);
    TotalWater = sum(SessionData.Custom.RewardReceivedTotal, 'omitnan'); %sums water volume from struct
        DataTable(Idx2,3) = num2cell(TotalWater);
    CorrectTrials = mean(SessionData.Custom.ResponseCorrect(26:end), 'omitnan'); %determines correct trials from struct
        DataTable(Idx2,4) = num2cell(CorrectTrials);
    AverageResponseTime = mean(SessionData.Custom.ResponseTime(26:end), 'omitnan'); %determines avg response from struct
        DataTable(Idx2,5) = num2cell(AverageResponseTime);

%% Determine False Alarms
%identify oddball trials (early withdrawal from center port)
oddErrors = find(isnan(SessionData.Custom.ResponseCorrect)); 

%identify oddball trials (early withdrawal from correct reward port)
oddErrors_2 = find(isnan(SessionData.Custom.RewardStartTime)); 

%pull these data out of the structre
noSignalTrials = SessionData.Custom.EmbedSignal; 
incorrectTrials = SessionData.Custom.ResponseCorrect; 
catchTrials = SessionData.Custom.CatchTrial; 

CorrectTrials = find(SessionData.Custom.ResponseCorrect); 
CorrectTrials = setdiff(CorrectTrials, oddErrors); 
SignalTrials = find(SessionData.Custom.EmbedSignal); 
SignalTrials = setdiff(SignalTrials, oddErrors); 
allCorrectSignalTrials = intersect(CorrectTrials,SignalTrials);  

tempErrors_1 = find(isnan(SessionData.Custom.ResponseCorrect(26:end)));
tempCorrectTrials = SessionData.Custom.ResponseCorrect(26:end); 
tempCorrectTrials(tempErrors_1) = 0; 

tempErrors_2 = find(isnan(SessionData.Custom.ResponseCorrect(26:end)));
tempSignalTrials = SessionData.Custom.EmbedSignal(26:end-1);
tempSignalTrials(tempErrors_2) = 0; 

correctSignalTrials = intersect(find(tempCorrectTrials), find(tempSignalTrials)); 
pctCorrectSignalTrials = nnz(correctSignalTrials)/nnz(tempSignalTrials); 

%find true incorrect and no signal trials in cleaned data
incorrectTrials = find(incorrectTrials == 0); 
noSignalTrials = find(noSignalTrials == 0);
catchTrials = find(catchTrials);

%remove indices of oddball trials
noSignalTrials = setdiff(noSignalTrials, oddErrors); 
incorrectTrials = setdiff(incorrectTrials, oddErrors); 
catchTrials = setdiff(catchTrials, oddErrors); %ignore catch trials where mouse did not stay in center


incorrectTrials = setdiff(incorrectTrials, oddErrors_2);

catchTrials = setdiff(catchTrials, incorrectTrials); %ignore catch trials when incorrect



% incorrectTrials = setdiff(incorrectTrials, oddErrors_2); 

falseAlarms = intersect(noSignalTrials, incorrectTrials); 

% there are both signal and no signal catch trials. 
noSignalCatchTrials = intersect(catchTrials, noSignalTrials);
signalCatchTrials = intersect(catchTrials, allCorrectSignalTrials);


falseAlarms(falseAlarms < 26) = []; 
noSignalTrials(noSignalTrials < 26) = []; 

if isempty(falseAlarms)
    falseAlarmPct = nan;
    meanfalseAlarmInvestment = nan;  
    pctHiFAinvestment = nan; 
else
    
falseAlarmPct = length(falseAlarms)/length(noSignalTrials)*100; 
falseAlarmInvestments = nan(1, length(falseAlarms)); 

for idx = 1:length(falseAlarms)

    if SessionData.SettingsFile.GUI.Ports_LMR == 321;
        
        if isfield(SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events, 'Port3Out') 
       
        falseAlarmInvestments(idx) = SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events.Port3Out(end) - SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        falseAlarmInvestments(idx) = nan; 
        
        end
        
    end
    
    if SessionData.SettingsFile.GUI.Ports_LMR == 123;
        
        if isfield(SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events, 'Port1Out') 
    
        falseAlarmInvestments(idx) = SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events.Port1Out(end) - SessionData.RawEvents.Trial{1, falseAlarms(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        falseAlarmInvestments(idx) = nan; 
        
        end
            
            
    end
end

falseAlarmInvestments(falseAlarmInvestments < 2) = []; 

meanfalseAlarmInvestment = nanmean(falseAlarmInvestments); 

pctHiFAinvestment = nnz(falseAlarmInvestments > nanmedian(falseAlarmInvestments))/length(falseAlarmInvestments); 

end

% noSignalCatchTrials 
% 
% signalCatchTrials 


if isempty(noSignalCatchTrials)
    
    noSignalCatchTrialInvestments = []; 
    
  else
    
    noSignalCatchTrials(noSignalCatchTrials > length(SessionData.RawEvents.Trial(1, :))) = []; %get rid of the last one, since the chamber stops recording so we don't get the timestamp. 
    
end

if ~isempty(noSignalCatchTrials)  
    
    for idx = 1:length(noSignalCatchTrials)

    if SessionData.SettingsFile.GUI.Ports_LMR == 321;
        
        if isfield(SessionData.RawEvents.Trial{1, noSignalCatchTrials(idx)}.Events, 'Port1Out') 
        
            noSignalCatchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, noSignalCatchTrials(idx)}.Events.Port1Out(end) - SessionData.RawEvents.Trial{1, noSignalCatchTrials(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        noSignalCatchTrialInvestments(idx) = nan; 
        
        end
        
    end
    
    if SessionData.SettingsFile.GUI.Ports_LMR == 123;
        
        if isfield(SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events, 'Port3Out')

        
        noSignalCatchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, noSignalCatchTrials(idx)}.Events.Port3Out(end) - SessionData.RawEvents.Trial{1, noSignalCatchTrials(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        noSignalCatchTrialInvestments(idx) = nan; 
        
        end
        
    
    end
    end
end


if isempty(signalCatchTrials)
    
    signalCatchTrialInvestments = []; 
    
  else
    
    signalCatchTrials(signalCatchTrials > length(SessionData.RawEvents.Trial(1, :))) = []; %get rid of the last one, since the chamber stops recording so we don't get the timestamp. 
    
end

if ~isempty(signalCatchTrials)  
    
    for idx = 1:length(signalCatchTrials)

    if SessionData.SettingsFile.GUI.Ports_LMR == 321;
        
        if isfield(SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events, 'Port3Out') 
        
            signalCatchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events.Port3Out(end) - SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        signalCatchTrialInvestments(idx) = nan; 
        
        end
        
    end
    
    if SessionData.SettingsFile.GUI.Ports_LMR == 123;
        
        if isfield(SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events, 'Port1Out')

        
        signalCatchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events.Port1Out(end) - SessionData.RawEvents.Trial{1, signalCatchTrials(idx)}.Events.GlobalTimer1_Start;
        
        else
            
        signalCatchTrialInvestments(idx) = nan; 
        
        end
        
    
    end
    end
end


    
    catchTrialInvestments = [noSignalCatchTrialInvestments signalCatchTrialInvestments]; 
    

    if isempty(catchTrials)  
    
    meanCatchTrialInvestment = nan; 
    
    else
        
    meanCatchTrialInvestment = nanmean(catchTrialInvestments); 
    
    end
    
    
    
    
    
    
    
    
% if isempty(catchTrials)  
%     
%     meanCatchTrialInvestment = nan; 
%     
% else
%     
%     catchTrials(catchTrials > length(SessionData.RawEvents.Trial(1, :))) = []; %get rid of the last one, since the chamber stops recording so we don't get the timestamp. 
%     
% for idx = 1:length(catchTrials)
% 
%     if SessionData.SettingsFile.GUI.Ports_LMR == 321;
%         
%         if isfield(SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events, 'Port3Out') 
%         
%         catchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events.Port3Out(end) - SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events.GlobalTimer1_Start;
%         
%         else
%             
%         catchTrialInvestments(idx) = nan; 
%         
%         end
%         
%     end
%     
%     if SessionData.SettingsFile.GUI.Ports_LMR == 123;
%         
%         if isfield(SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events, 'Port1Out')
% 
%         
%         catchTrialInvestments(idx) = SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events.Port1Out(end) - SessionData.RawEvents.Trial{1, catchTrials(idx)}.Events.GlobalTimer1_Start;
%         
%         else
%             
%         catchTrialInvestments(idx) = nan; 
%         
%         end
%         
%     
%     end
%     
% end
% 
% 
% meanCatchTrialInvestment = nanmean(catchTrialInvestments); 

% end
        DataTable(Idx2,6) = num2cell(falseAlarmPct);
        DataTable(Idx2,7) = num2cell(meanfalseAlarmInvestment);

        DataTable(Idx2,8) = num2cell(pctHiFAinvestment); 
        DataTable(Idx2,9) = num2cell(pctCorrectSignalTrials);
        DataTable(Idx2,10) = num2cell(meanCatchTrialInvestment);
        
%% Determine which phase the animal is on
    if SessionData.SettingsFile.GUI.ContinuousTable.SignalLimits(1) >=64
        PhaseType = "Phase 1 (one light) or 2 (two light)";
    elseif SessionData.SettingsFile.GUI.ContinuousTable.SignalLimits(1) >=35 & SessionData.SettingsFile.GUI.FeedbackDelay == 0
        PhaseType = "Phase 3 (fine discrim)";
    elseif SessionData.SettingsFile.GUI.FeedbackDelay > 0
        PhaseType = "Phase 4 (No Catch)";
    else SessionData.SettingsFile.GUI.FeedbackDelay > 0 & SessionData.SettingsFile.GUI.PercentCatch > 0
        PhaseType = "Phase 5 (Catch Trials)";
end 
        DataTable(Idx2,11) = cellstr(PhaseType);
end 
%% Generate GUI Panels
f = figure;
  p1 = uipanel('Title','Main Table','BackgroundColor','white','Position', ...
   [0 .5 .5 .5], 'TitlePosition','centertop','BorderType','line', 'BorderWidth', 2, 'FontSize', 24, 'Scrollable',"on");
  p2 = uipanel('Title','Plots','BackgroundColor','white','Position', ...
   [0 0 1 .5], 'TitlePosition','centertop','BorderType','line', 'BorderWidth', 2, 'FontSize', 24);
  
%% Generate Table for panel 1
  uit = uitable('Parent',p1,'Units','normalized','Position', [0 0 1 1],...
'ColumnName', {'DaysSinceStart','TotalTrials','WaterVolume', 'PercentCorrect', 'AverageResponse', 'FalseAlarmPct', 'meanFalseAlarmInvest','pctHiFAinvestment','pctCorrectSignalTrials', 'meanCatchTrialInvestment', 'Which phase was the mouse on'}, ...
    'RowName', Dates, 'ColumnWidth', 'auto');
uit.Data = DataTable;

%% Generate a set of summary graphs for panel 2
X = cell2mat(DataTable(:,1));
y1 = cell2mat(DataTable(:,2));
y2 = cell2mat(DataTable(:,3));
y3 = cell2mat(DataTable(:,4));
y4 = cell2mat(DataTable(:,5));
t = tiledlayout(p2,2,2)

% Top left plot
nexttile(t); 
    plot(X,y1,'o','MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black'), 
    title('Trial Number'), 
    ylabel('# of trials')
    xlabel('days since start')
% Bottom left plot
nexttile(t); 
    plot(X,y2,'o','MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black'), 
    title('Total Water'), 
    ylabel('WaterVolume'),
    xlabel('days since start')
% Top right plot
nexttile(t); 
    plot(X,y3,'o','MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black'), 
    title('SuccessRate'), 
    ylabel('Percent Correct'),
    xlabel('days since start'),
    yline(.65, 'LineStyle','--','Color','#FFA500', 'LineWidth', 1.5),
    yline(.7,'LineStyle','-','Color','#FFFF00', 'LineWidth', 1.5)
% Bottom right plot
nexttile(t); 
    plot(X,y4,'o','MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black') 
    title('Response Rate'), 
    ylabel('Average Response (s)')
    xlabel('days since start')
   
