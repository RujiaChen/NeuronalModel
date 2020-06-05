close all;
NumCueDelay=cell(1,2);
for imonkey=1:2    
    if imonkey==1
        dateUsed={'111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
    else
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
    end
    figure;
    for idate=1:numel(dateUsed)
        load(['Z:\RujiaChen\Results\CorrectTrialInfo_' dateUsed{idate} '_new.mat']);
        load(['Z:\RujiaChen\Results\CorrectTrialParam_' dateUsed{idate} '.mat']);
        load(['Z:\RujiaChen\Results\flanker_TrlParam_' dateUsed{idate} '.mat']);
        idx1 = find(TrlParam.trial_response==-105&TrlParam.is_mouse_trial==0);
        idx2 = find(CorrectTrialParam.is_mouse_trial==0);
        if imonkey==1||(imonkey==2&&idate<=7)
            CueDelayDuration=Correct.ArrayOnsetT-Correct.DelayT;
            ArrayDelay=TrlParam.holdshape_duration(idx1)';   %-(Correct.ArrayOnsetT(idx2)-Correct.ArrayT(idx2));
        else
            CueDelayDuration=Correct.ArrayT-Correct.DelayT;
            ArrayDelay=TrlParam.holdshape_duration(idx1);
        end
        subplot(4,4,idate);
        NumCueDelay{imonkey}(idate,1)=min(ArrayDelay);
        NumCueDelay{imonkey}(idate,2)=max(ArrayDelay);
        hist(CueDelayDuration); hold on;
%         hist(ArrayDelay); hold on;
        
    end
end