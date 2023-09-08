function Bin_Indices = binFA(SessionData, binsize, easyTrials) % x = a session to bin, y = output double (should be 20 values)
    newTime = SessionData.TrialEndTimestamp(easyTrials:end) / 60; %convert seconds to minutes
    timeAdjusted = newTime - (newTime(1,1)); %standardize times to 0 seconds for convenience
    timeblocks = [0:binsize:(60-binsize)]; %set timeblocks for indices
    Bin_Indices = [];
    if max(timeAdjusted) < max(timeblocks)
        Bin_Indices = vertcat(Bin_Indices, [0 0]');
    else
        for idx = timeblocks
            if idx == 0;
                findEnd = find(timeAdjusted == interp1(timeAdjusted, timeAdjusted, 3, 'nearest')); %finds trial nearest to 15 min
                block = [1:findEnd];
                Bin_Indices = [min(block), max(block)]';

            elseif idx == max(timeblocks);
                 %this is built in in case the adjusted time is less than 60 somehow. it basically find the last trial independent of time.
                if max(timeAdjusted) > max(timeblocks) & max(timeAdjusted) < max(timeblocks) + binsize;
                     findBeginning = findEnd + 1; 
                     findEnd = length(SessionData.Custom.ResponseCorrect) - easyTrials;
                     block = [block findBeginning:findEnd];
                     Bin_Indices = [Bin_Indices [findBeginning findEnd]'];

                else floor(interp1(timeAdjusted, timeAdjusted, [idx+binsize], 'nearest')) | ceil(interp1(timeAdjusted, timeAdjusted, [idx+binsize], 'nearest'))  == 60;
                     findBeginning = findEnd + 1; 
                     findEnd = find(timeAdjusted == interp1(timeAdjusted, timeAdjusted, [idx+binsize], 'nearest')); %finds range for 30, 45, 60
                     block = [block findBeginning:findEnd];
                     Bin_Indices = [Bin_Indices [findBeginning findEnd]'];
              
                end 

            else idx ~= 0;
                findBeginning = findEnd + 1; 
                findEnd = find(timeAdjusted == interp1(timeAdjusted, timeAdjusted, [idx+binsize], 'nearest')); %finds range for 30, 45, 60
                block = [block findBeginning:findEnd];
                Bin_Indices = [Bin_Indices [findBeginning findEnd]'];

            end 
        end 
        % timelineFA = intersect(intersect(noSignalTrials, incorrectTrials,'stable'),block, 'stable')
    end 
        