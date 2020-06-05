%%% plot the behavioral results for each session
% dateUsed={'112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  %'042418','051218','050118','050418','050518','042118','051718',
clear all;
close all;
dateUsed{1}={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
dateUsed{2}={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
meanCR=zeros(2, 2, 3);
eeCR=zeros(2, 2, 3);
HitRate=cell(2);
for imonkey=1:2    
%     figure;
    for icueTy=1:2
        isexo=2-icueTy;
        CorrectRate=zeros(numel(dateUsed{imonkey}), 3);
        for idate=1:numel(dateUsed{imonkey})
            date=dateUsed{imonkey}{idate};
            load(['Z:\RujiaChen\Results\flanker_TrlParam_' date '.mat']);          
            
%             rr=unique(TrlParam.trial_Response);
            %         [hh, pp]=hist(TrlParam.trial_response, rr);
            idxx=TrlParam.isexocue==isexo;
            hit1=sum(TrlParam.trial_response==-105&idxx);
            err1=sum(TrlParam.trial_response==132&idxx);
            err2=sum(ismember(TrlParam.trial_response, [-15, -16])&idxx);
            allTrl=sum(TrlParam.trial_response~=129&idxx);
            
            CorrectRate(idate, 1)=hit1/allTrl;
            CorrectRate(idate, 2)=err1/allTrl;
            CorrectRate(idate, 3)=err2/allTrl;
        end
%         subplot(1,2,icueTy);
%         bar(CorrectRate, 'stacked'); hold on;
%         line([0 numel(dateUsed{imonkey})+1], [0.8 0.8], 'linestyle', '--', 'color', 'k'); hold on;
%         set(gca,'box' ,'off');
%         xlabel('Session Number');
%         ylabel('Behavioral Performance');
        meanCR(imonkey, icueTy,:)=mean(CorrectRate,1);
        eeCR(imonkey, icueTy,:)=std(CorrectRate,[],1);  %/sqrt(size(CorrectRate,1));
        HitRate{imonkey,icueTy}=CorrectRate(:,1);
    end
end

for imonkey=1:2
   [tt(imonkey),pp(imonkey)]=ttest2(HitRate{imonkey,1},HitRate{imonkey,2});
end
figure;
clr2=[1 1 1; 0.5 0.5 0.5];
%%%%% plot the hit rates averaged across sessions
for itype=1:2
    for imonkey=1:2        
        bb(itype)=bar(itype+(imonkey-1)*3, squeeze(meanCR(imonkey,itype,1)),'facecolor',clr2(itype,:),'edgecolor',[0.5 0.5 0.5],'linewidth',1.5); hold on;
        errorbar(itype+(imonkey-1)*3, meanCR(imonkey,itype,1), eeCR(imonkey,itype,1),'color',[0.4 0.4 0.4],'linewidth',1.5); hold on;
    end
end
legend(bb,'Exo','Endo');
xticks([1.5 4.5]);
xticklabels({'Mikey','Vasco'});
ylabel('Hit rate');
set(gca,'box','off');


%% plot the averaged ratio of different reactions: correct/early response/late+error response
close all;
figure;
clr2=[0 0 0; 0.5 0.5 0.5];
for ii=1:2
    errorbar(1:3, meanCR(ii,:), eeCR(ii,:), 'color',clr2(ii,:),'linewidth', 1.5); hold on;
end
legend('Exo', 'Endo','box','off');
xticks([1 2 3]);
xticklabels({'Correct','Early','Late/Error'});
ylabel('Occurence ratio');
set(gca,'box','off');


%% compare the reaction time in different conditions for each session
close all;
folder='Z:\RujiaChen\Results\';
dateUsed{1}={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
dateUsed{2}={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
ExoRT=cell(1,2);
EndoRT=cell(1,2);
for imonkey=1:2
    figure;
    ExoRT{imonkey}=[];
    EndoRT{imonkey}=[];
    for idate=1:numel(dateUsed{imonkey})
        %     figure;
        subplot(3,5,idate);
        RTime=cell(1,2);
        date=dateUsed{imonkey}{idate};
        load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
        load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat']);
        idx1=find(CorrectTrialParam.is_mouse_trial==0&CorrectTrialParam.isexocue==1);
        Time_offset=transpose(Correct.ArrayOnsetT-Correct.ArrayT);
        if imonkey==1||(imonkey==2&&idate<=7)
            RTime{1}=CorrectTrialParam.RT(idx1)-Time_offset(idx1);
        else
            RTime{1}=CorrectTrialParam.RT(idx1);
        end
        ExoRT{imonkey}=[ExoRT{imonkey}; RTime{1}];
        idx2=find(CorrectTrialParam.is_mouse_trial==0&CorrectTrialParam.isexocue==0);
        if imonkey==1||(imonkey==2&&idate<=7)
            RTime{2}=CorrectTrialParam.RT(idx2)-Time_offset(idx2);
        else
            RTime{2}=CorrectTrialParam.RT(idx2);
        end
        EndoRT{imonkey}=[EndoRT{imonkey}; RTime{2}];
        plot(idx1, RTime{1}, 'r*'); hold on;
        plot(idx2, RTime{2}, 'b*'); hold on;
        xlim([0 400]);
        ylim([0 600]);
    end
end
%% plot the distribution of reaction times for different conditions
close all;
tt=[0 0];
for imonkey=1:2
    figure;
    subplot(1,2,1)
    cc=-25:50:500;
    hh=histogram(ExoRT{imonkey},cc);
    subplot(1,2,2)
    hh2=histogram(EndoRT{imonkey},cc);
    tt(imonkey)=ttest2(ExoRT{imonkey},EndoRT{imonkey});
    
    figure;
    bb=[hh.Values' hh2.Values'];
    xx=hh.BinEdges(1:end-1)+hh.BinWidth/2;
    bar(xx,bb);
    xlabel('Reaction time (ms)');
    ylabel('Occurence');
    legend('Exo','Endo');
    set(gca,'box','off');
end

%% plot the session averaged results
close all;
figure;
clr2=[1 1 1; 0.7 0.7 0.7];
tt=[0 0];
pvalue=[0 0];
for imonkey=1:2    
    idx=ExoRT{imonkey}>100;    
    mmExo=mean(ExoRT{imonkey}(idx,:));
    eeExo=std(ExoRT{imonkey}(idx,:));  %/sqrt(numel(ExoRT{imonkey}));
%     bar(1,mmExo); hold on;
    bb(itype)=bar(1+(imonkey-1)*3, mmExo,'facecolor',clr2(1,:),'edgecolor',clr2(2,:),'linewidth',1.5); hold on;
    errorbar(1+(imonkey-1)*3, mmExo, eeExo,'color',[0.6 0.6 0.6],'linewidth',1.5); hold on;
    
    idx2=EndoRT{imonkey}>100;
    [tt(imonkey),pvalue(imonkey)]=ttest2(ExoRT{imonkey}(idx),EndoRT{imonkey}(idx2));
    mmEndo=mean(EndoRT{imonkey}(idx2,:));
    eeEndo=std(EndoRT{imonkey}(idx2,:));  %/sqrt(numel(EndoRT{imonkey}));
%     bar(2,mmEndo); hold on;
    bb(itype)=bar(2+(imonkey-1)*3, mmEndo,'facecolor',clr2(2,:),'edgecolor',clr2(2,:),'linewidth',1.5); hold on;
    errorbar(2+(imonkey-1)*3, mmEndo, eeEndo,'color',[0.6 0.6 0.6],'linewidth',1.5); hold on;
    xticks([1.5 4.5]);
    xticklabels({'Mikey','Vasco'});
    ylabel('Reaction time (ms)')
    set(gca,'box','off');
end
