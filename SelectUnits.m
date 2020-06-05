%% select the array responsive units which can significantly distinguish the cue-on and cue-away conditions 
close all;
clear all;
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};  % for Vasco
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
% dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco
% dateUsed={'111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape
monkeyset={'Mikey','Vasco'};
TimeCue=-600:1200;
% TimeArray=-905:905;
TimeArray=-505:1305;  % for mikey 
folder='Z:\RujiaChen\Results\';
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

%% plot the distribution of different cell type
% close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
imonkey=1;
monkey=monkeyset{imonkey};
if imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};     
%     dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
    channel=1:96;  % start with pulvinar
elseif imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
    load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');    
end
figure;
for idate=1:numel(dateUsed)   %[1: 5  7: 9]   %
    date=dateUsed{idate};
    load([folder 'bVisResp_' date '_5sigma.mat']);  %SingleUnit
    load([folder 'bMotorResp_' date '_5sigma.mat']);
    load([folder 'bVisMotor_' date '_5sigma.mat']);
    load([folder 'CueOn_' date '_5sigma.mat']);
    load([folder 'RFOn_' date '_5sigma.mat']);
    load([folder 'SaccadeOn_' date '_5sigma.mat']);
%     load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);
    for iarea=1:3        
        if imonkey==1
            isession=find(strcmp(RawInfo(1,:), date)==1);
            channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
            ContactID =unique(SpikeTrain.channel(channelID));
        else
            channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
            ContactID =unique(SpikeTrain.channel(channelID));
        end
       
        subplot(1,3,iarea);
        %%%%% for MUA activity, unsorted data
        idxVis=bVisResp(channelID)>0;
        idxVisMotor=bVisMotor(channelID)>0;   %sum(bVisResp(channelID,:),2)>0|
        idxMotor=bMotorResp(channelID)>0;
        idx0=idxVis==0&idxVisMotor==0&idxMotor==0;
        yy=1:numel(channelID);
        plot(idate*ones(1,sum(idxMotor)),yy(idxMotor), 'go','markerfacecolor','g'); hold on;
        plot(idate*ones(1,sum(idxVis)),yy(idxVis), 'ro','markerfacecolor','r'); hold on;
        plot(idate*ones(1,sum(idxVisMotor)),yy(idxVisMotor), 'bo','markerfacecolor','b'); hold on;
        plot(idate*ones(1,sum(idx0)),yy(idx0), 'ko','markerfacecolor','w'); hold on;
        set(gca,'Ydir','reverse');                        
    end
end


%% plot the distribution of RF and saccade field
% close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
imonkey=2;
monkey=monkeyset{imonkey};
if imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};     
    channel=1:96;  % start with pulvinar
elseif imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
    load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');    
end

Selection={'_SingleUnit','_5sigma'};
iselect=2;

for igroup=2:3
    figure;
    clr5=[0 0 0; eye(3);1 1 0];
    for idate=1:numel(dateUsed)   %[1: 5  7: 9]   %
        date=dateUsed{idate};
        load([folder 'bVisResp_' date Selection{iselect} '.mat']);  %SingleUnit
        load([folder 'bMotorResp_' date Selection{iselect} '.mat']);
        load([folder 'bVisMotor_' date Selection{iselect} '.mat']);
        load([folder 'CueOn_' date Selection{iselect} '.mat']);
        load([folder 'RFOn_' date Selection{iselect} '.mat']);
        load([folder 'SaccadeOn_' date Selection{iselect} '.mat']);
        if iselect==1
            load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']); % Single units
        elseif iselect==2
            load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']); % 5sigma MUA
        end
       
        for iarea=1:3
            if imonkey==1
                isession=find(strcmp(RawInfo(1,:), date)==1);
                channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
                ContactID =unique(SpikeTrain.channel(channelID));
            else
                channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
                ContactID =unique(SpikeTrain.channel(channelID));
            end
            if igroup==1
                Vresp=bVisResp(channelID).*CueOn(channelID);
                Mresp=bVisResp(channelID).*SaccadeOn(channelID);
            elseif igroup==2
                Vresp=bVisMotor(channelID).*CueOn(channelID);
                Mresp=bVisMotor(channelID).*SaccadeOn(channelID);
            elseif igroup==3
                Vresp=bMotorResp(channelID).*CueOn(channelID);
                Mresp=bMotorResp(channelID).*SaccadeOn(channelID);
%             elseif igroup==4
%                 celltype=(bVisResp(channelID)+bVisMotor(channelID)).*CueOn(channelID);
            end
            subplot(2,3,iarea);
            %%%%% for MUA activity, unsorted data
            for ii=0:4
                idx1=Vresp==ii;
                yy=1:numel(channelID);
                plot(idate*ones(1,sum(idx1)),yy(idx1), 'o','markeredgecolor',clr5(ii+1,:),'markerfacecolor',clr5(ii+1,:)); hold on;
            end
            
            set(gca,'Ydir','reverse');
            
            subplot(2,3,3+iarea);
            %%%%% for MUA activity, unsorted data
            for ii=0:4
                idx1=Mresp==ii;
                yy=1:numel(channelID);
                plot(idate*ones(1,sum(idx1)),yy(idx1), 'o','markeredgecolor',clr5(ii+1,:),'markerfacecolor',clr5(ii+1,:)); hold on;
            end
            
            set(gca,'Ydir','reverse');
        end
    end
end
