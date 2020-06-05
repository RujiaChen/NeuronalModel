%% select the array responsive units which can significantly distinguish the cue-on and cue-away conditions 
close all;
clear all;
TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape
monkeyset={'Mikey','Vasco'};
TimeCue=-600:1200;
TimeArray=-505:1305;  % for mikey 
folder='Z:\Rujia\Results\';
TimeZero='ToArray';
CueCondition='Exo';

for imonkey=1:2
    monkey=monkeyset{imonkey};
    if strcmp(monkey, 'Mikey')
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
    elseif strcmp(monkey, 'Vasco')
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
    end
    
    for idate=1:numel(dateUsed)
        date=dateUsed{idate};
        load([folder 'CorrectTrialParam_' date '.mat']);
        % load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
        load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);
%         load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
        load(['Z:\RujiaChen\Results\SaccadeTime_'  TimeZero '_' CueCondition '_' date '.mat']);
        Ncell=length(SpikeTrain.channel);
        RFOn=zeros(Ncell,1);
        SaccadeOn=zeros(Ncell,1);
        CueOn=zeros(Ncell,1);
        
        cueMeanResp = zeros(Ncell,4);
        motorMeanResp = zeros(Ncell,4);
        bVisResp=zeros(Ncell,1);
        bVisMotor=zeros(Ncell,1);
        bMotorResp=zeros(Ncell,1);
        for ich=1:Ncell    %1:96
            num=0;
            motorResp=[];BeforeMotor=[];
            bMSResp=zeros(1,4);
            bResp=zeros(1,4);
            bCueResp=zeros(1,4);
            for ipos=1:2
                for icueP=1:2
                    num=num+1;
                    ntrl=0;
                    idxTrl=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icueP) &CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0);
                    ArrayResp=squeeze(SpikeTrain.arrayAlign(idxTrl,ich,:));% spike response
                    CueResp=squeeze(SpikeTrain.cueAlign(idxTrl,ich,:));
                    %%%% test whether the current channel is responsive to array
                    if numel(idxTrl)>=15
                        Baseline=squeeze(mean(CueResp(:,TimeCue>-200&TimeCue<0),2));
                        cResp=squeeze(mean(CueResp(:,TimeCue>50&TimeCue<250),2));
                        mResp=squeeze(mean(ArrayResp(:,TimeArray>50&TimeArray<250),2));  % 50-300 for LFP ; 50-600 for spike to see whether array responsive
                        preArray=squeeze(mean(ArrayResp(:,TimeArray>-200&TimeArray<0),2));
                        for itrl=1:numel(idxTrl)
                            if SaccadeTime{ipos, icueP}(itrl)>400
                                ntrl=ntrl+1;
                                motorT=min([SaccadeTime{ipos, icueP}(itrl),1200]);
                                motorResp(ntrl)=squeeze(mean(ArrayResp(itrl,TimeArray>motorT-100&TimeArray<=motorT+100),2));  %[-250 0]   %Time>300&Time<600
                                BeforeMotor(ntrl)=squeeze(mean(ArrayResp(itrl,TimeArray>motorT-350&TimeArray<=motorT-200),2));
                            end
                        end
                        
                        [~,bMSResp(num)]=signrank(motorResp(1:ntrl),BeforeMotor(1:ntrl),'tail','right','alpha',0.05);
                        [~,bResp(num)]=signrank(mResp,preArray,'tail','right','alpha',0.05);
                        [~,bCueResp(num)]=signrank(cResp,Baseline,'tail','right','alpha',0.05);
                    end                    
                    cueMeanResp(ich,num) = mean(cResp-Baseline);   
                    motorMeanResp(ich,num)= mean(motorResp - BeforeMotor);
                end
            end

            if sum(bMSResp)>0
                [~,SaccadeOn(ich)]= max(motorMeanResp(ich,:));
            end
            if sum(bResp)>0||sum(bCueResp)>0
                visResp=1;
                [~,RFOn(ich)]= max(cueMeanResp(ich,:));
            else
                visResp=0;
            end
            if sum(bCueResp)>0
                [~,CueOn(ich)]= max(cueMeanResp(ich,:));
            end
            bMotorResp(ich)=sum(bMSResp)>0&visResp==0;
            bVisResp(ich)=sum(bMSResp)==0&visResp==1;
            bVisMotor(ich)=sum(bMSResp)>0&visResp==1;
        end
        
%         save([folder 'CueOn_' date '_SingleUnit.mat'],'CueOn', '-v7.3');
%         save([folder 'RFOn_' date '_SingleUnit.mat'],'RFOn', '-v7.3');        
%         save([folder 'SaccadeOn_' date '_SingleUnit.mat'],'SaccadeOn', '-v7.3');
%         save([folder 'bVisMotor_' date '_SingleUnit.mat'],'bVisMotor', '-v7.3');
%         save([folder 'bVisResp_' date '_SingleUnit.mat'],'bVisResp','-v7.3');
%         save([folder 'bMotorResp_' date '_SingleUnit.mat'],'bMotorResp','-v7.3');
    
        save([folder 'CueOn_' date '_5sigma.mat'],'CueOn', '-v7.3');
        save([folder 'SaccadeOn_' date '_5sigma.mat'],'SaccadeOn', '-v7.3');
        save([folder 'RFOn_' date '_5sigma.mat'],'RFOn', '-v7.3');
        save([folder 'bVisMotor_' date '_5sigma.mat'],'bVisMotor', '-v7.3');
        save([folder 'bVisResp_' date '_5sigma.mat'],'bVisResp','-v7.3');
        save([folder 'bMotorResp_' date '_5sigma.mat'],'bMotorResp','-v7.3');
    end
end

