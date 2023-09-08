%% script to combine Bpod HALIP data across multiple mice within a given treatment condition


path = {

'R:\PBS\Parker_Lab\JustinAnair\Data\HALIP\DATCre_TRPV1KO_AAVTRPV1\Bilateral Expression\Veh_Veh'
'R:\PBS\Parker_Lab\JustinAnair\Data\HALIP\DATCre_TRPV1KO_AAVTRPV1\Bilateral Expression\Veh_Cap'
% 'R:\PBS\Parker_Lab\JustinAnair\Data\HALIP\B6'
% 'R:\PBS\Parker_Lab\JustinAnair\Data\HALIP\TRPV1Ko'
};

deltaBins = [-40, -5, 5,15,20, 25];
color = {'k' 'm' 'r' 'c' 'b' 'g'};
line = {'-k' '-m' '-r' '-c' '-b' '-g'};
legends = {'Veh_Veh','Veh_Cap','Hal_Cap','VEH + CAP', 'HALDOL + CAP','HALDOL + VEH'}; 


for x = 1:length(path)
    
    cd(path{x}); 
    fileNames = dir;
    fileNames = struct2cell(fileNames); 
    
    
    
    fIndices = find(contains(fileNames(1,:), 'f')); 
    mIndices = find(contains(fileNames(1,:), 'm'));
    mouseIndices = [fIndices mIndices];
    
    combinedStruct.FalseAlarmRate = nan(1,length(mouseIndices)); 
    combinedStruct.FalseAlarmInvest = nan(1,length(mouseIndices));
    combinedStruct.normFalseAlarmInvestments = []; 
    combinedStruct.PercentCorrect = nan(1,length(mouseIndices)); 
    combinedStruct.TotalWater = nan(1,length(mouseIndices)); 
    combinedStruct.MeanCatchTrialInvest = nan(1,length(mouseIndices)); 
    combinedStruct.slope = nan(1,length(mouseIndices)); 
    combinedStruct.bias = nan(1,length(mouseIndices)); 
    combinedStruct.pctHiFAinvestment = nan(1,length(mouseIndices)); 
    combinedStruct.missTrialInvest = nan(1,length(mouseIndices));  
    combinedStruct.psychData = []; 
    combinedStruct.Exclusions.nMice = []; 
    combinedStruct.nTrials_psych = []; 
    combinedStruct.fit_results = [];
    combinedStruct.Exclusions.rawTrialTotal = [];
    combinedStruct.Exclusions.nSessions = [];
    combinedStruct.Exclusions.nSessionsClean = [];
    combinedStruct.Exclusions.SessionsExcludedPct = []; 
    combinedStruct.Binned.FArate = [];
    combinedStruct.Binned.FAcount = [];
    combinedStruct.Binned.concatFArate = [];
    combinedStruct.Binned.concatFAinvest = [];

    
    for y = 1:length(mouseIndices)  
        
        cd(path{x}); 
        cd(fileNames{1, mouseIndices(y)}); 
        cd('Analysis'); 
        load('OutputData.mat'); 
        
    combinedStruct.FalseAlarmRate(y) = OutputData.FalseAlarmRate; 
    combinedStruct.FalseAlarmInvest(y) = OutputData.FalseAlarmInvest;
    combinedStruct.normFalseAlarmInvestments = [combinedStruct.normFalseAlarmInvestments; OutputData.normfalseAlarmInvestments'];
    combinedStruct.PercentCorrect(y) = OutputData.PercentCorrect;
    combinedStruct.TotalWater(y) = OutputData.TotalWater; 
    combinedStruct.MeanCatchTrialInvest(y) = OutputData.MeanCatchTrialInvest;
    combinedStruct.pctHiFAinvestment(y) = OutputData.pctHiFAinvestment;
    combinedStruct.missTrialInvest(y) = OutputData.missTrialInvest;
    combinedStruct.psychData = [combinedStruct.psychData; OutputData.psychData];
    combinedStruct.Exclusions.nSessions = [combinedStruct.Exclusions.nSessions; OutputData.Exclusions.nSessionsClean];
    combinedStruct.Exclusions.nSessionsClean = [combinedStruct.Exclusions.nSessionsClean; OutputData.Exclusions.nSessionsClean];
    combinedStruct.Exclusions.rawTrialTotal = [combinedStruct.Exclusions.rawTrialTotal; sum(OutputData.Exclusions.rawTrialTotal)];
    combinedStruct.Binned.FArate = [combinedStruct.Binned.FArate; OutputData.Binned.FArate];
    combinedStruct.Binned.FAcount = [combinedStruct.Binned.FAcount; OutputData.Binned.FAcount];
    combinedStruct.Binned.concatFArate = [combinedStruct.Binned.concatFArate; OutputData.Binned.concatbinFA];
    combinedStruct.Binned.concatFAinvest = [combinedStruct.Binned.concatFAinvest; OutputData.Binned.concatbinInvest];
    
    if ~isempty(OutputData.psychData)
    fit_results = psychFit(OutputData.psychData(:, 1), OutputData.psychData(:, 2), deltaBins);
   
    combinedStruct.slope(y) = fit_results.slope;
    combinedStruct.bias(y) = fit_results.bias;
    
    %uncomment to see individual fits
     % plot(fit_results.delta_fit, fit_results.pright_fit, '-k', 'LineWidth', 2); 
%    close all; 
    
    else
        
    combinedStruct.slope(y) = nan;
    combinedStruct.bias(y) = nan;
   
    end
    
    end
    %calculate session number, number of mice, excluded sessions and trials
    combinedStruct.Binned.FArate = mean(combinedStruct.Binned.FArate);
    combinedStruct.Binned.FAcount = sum(combinedStruct.Binned.FAcount, 'omitmissing');
    combinedStruct.nMice = length(mouseIndices); 
    combinedStruct.nTrials_psych = size(combinedStruct.psychData, 1); 
    combinedStruct.trialsExcludedPct = (1 - (combinedStruct.nTrials_psych / sum(combinedStruct.Exclusions.rawTrialTotal))) * 100;
    combinedStruct.Exclusions.SessionsExcludedPct = (1 - sum(combinedStruct.Exclusions.nSessionsClean) / sum(combinedStruct.Exclusions.nSessions)) * 100;
    
    %old psychometric fit
    combinedStruct.fit_results = psychFit(combinedStruct.psychData(:, 1), combinedStruct.psychData(:, 2), deltaBins);
    
    %new psychometric fit
%     combinedStruct.psych = psychometricFit(combinedStruct.psychData(:, 2),combinedStruct.psychData(:, 1), 0,deltaBins,1,0)
    
    %conversion?
%     pright_data = psych.perfPsych_binned   
%     pright_error =  psych.perfPsychJ_binnedSEM
%     pright_fit =  psych.fit.curve
%     delta_fit = psych.fit.curve
%     
%     
%     perfPsych_binned, psych.perfPsychJ_binnedSEM, psych.fit.curve, psych.fit.curve
    
    cd('..\'); 
    cd('..\'); 
    
    mkdir('Analysis')
    
    cd('Analysis'); 
    
    save('combinedStruct.mat', 'combinedStruct'); 

    %plot psychFit
    yneg = combinedStruct.fit_results.pright_data - combinedStruct.fit_results.pright_error(1,:);
    ypos = combinedStruct.fit_results.pright_error(2,:) - combinedStruct.fit_results.pright_data;
    
    
    hold on; 
    errorbar(combinedStruct.fit_results.delta_data,combinedStruct.fit_results.pright_data,yneg,ypos,'o', ...
            'MarkerSize',10,'MarkerEdgeColor',color{x},'MarkerFaceColor',color{x}, 'Color',color{x});
    ylabel('Signal Choice (%)'); xlabel('Signal-to-noise'); ylim([0 100]); 
    
    if x == 1; 
        
    p1 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2); 
    
    end
    
    if x == 2; 
        
    p2 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2); 
    end
    if x == 3; 
        
    p3 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2);     
    end
    if x == 4; 
        
    p4 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2); 
    end
%     if x == 5; 
        
%     p5 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2);
%     end
%     if x == 6; 
%         
%     p6 = plot(combinedStruct.fit_results.delta_fit, combinedStruct.fit_results.pright_fit, line{x}, 'LineWidth', 2);
%     
%     end
    
    allStruct{x} = combinedStruct; %allStruct places combinedStructs into a cell array; struct number corresponds to order of file paths
end

% get plot handles for desired data to have legend
hleglines = [p1(1) p2(1) p3(1)]; %p4(1) p5(1) p6(1)];

% create the legend
legend(hleglines,'Veh+Veh','Veh+Cap','Hal+Cap','VEH + CAP', 'HALDOL + CAP','HALDOL + VEH', 'Location', 'northwest' );

figure(2);
plot([1:length(allStruct{1,1}.FABinned)], allStruct{1, 1}.FABinned,"Color", "k")
hold on;
plot(allStruct{1, 2}.FABinned, "Color", "r")
% hist(allStruct{1,1}.normFalseAlarmInvestments(:,1),25,'black');
% h = findobj(gca,'Type','patch');
% h.CData = [0 0 0];;
% hold on;
% 
% hist(allStruct{1,2}.normFalseAlarmInvestments(:,1),25);
% h = findobj(gca,'Type','patch');
% h.CData = [1 0 0];
% hold on;
% 
% hist(allStruct{1,1}.normFalseAlarmInvestments(:,1),25);
% h.CData = [1 0 1];