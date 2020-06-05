%%%      Data preanalysis steps: 
%%%%   loadSpike,loadLFP, loadParams, SelectUnits, PopulationAnalysis
% clear all; 
close all;
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};  % for Vasco
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
% dateUsed={'112018', '112118', '010719', '010919', '011119', '011419', '011519', '011719', '011819', '012119'};  % for Vasco  '111918', 
% dateUsed={'111717', '111817','111917','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey

monkeyset={'Mikey','Vasco'};
imonkey=1;

if imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey 
    load('Z:\Rujia\Results\Mikey_RecordingInfo.mat');
    folder='Y:\anne\monkeys\pigz_plexondata\mikey_left_hemisphere\';  % for Mikey
elseif imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018','010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco
    folder='Y:\anne\monkeys\plexondata\';      
end
%%%% 011119 signal is not very good, spikes are too sparse
areaName='allArea';
for idate=1 %:numel(dateUsed)    
clear SpikeTrain; clear Correct;
date=dateUsed{idate};
% if ~exist(['Z:\RujiaChen\Results\Unsorted_SpikeTrain_' areaName  '_' date '_TriggerCorrected.mat'],'file')
fileName=[folder date '\merged_spl_noWB_001-5sigma-01.pl2'];
muaDataDirRoot=folder;  
sessionName='merged_spl_noWB_001-5sigma-01';
if imonkey==1
    isession=find(strcmp(RawInfo(1,:), date)==1);
    channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  %for Mikey 
elseif imonkey==2
    channel=1:96;   %for Vasco
end

isLoadSpikes=1;
isLoadMua=0;
isLoadLfp=0;
isLoadSpkc=0;
isLoadDirect=0;
spikeChannelPrefix = 'SPK_SPKC';
% spikeChannelPrefix = 'SPK';
spikeChannelsToLoad=channel;
muaChannelsToLoad=channel;
lfpChannelsToLoad = channel; 
spkcChannelsToLoad = channel;

directChannelsToLoad =3;
D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);

duration=D.blockStopTimes-D.blockStartTimes;
[~,gratingsTask3DIndices]=max(duration);

D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);
% blockName = strjoin(blockNames(gratingsTask3DIndices), '-');
% fprintf('Analyzing block names: %s.\n', blockName);

% %% read the events times in plexon
allEventTimes = cell2mat(D.events(1:11)');
allEventTimes = unique(allEventTimes);
BinaryCode=cell(size(allEventTimes));
EventCode=zeros(size(allEventTimes));
for i = 1:numel(allEventTimes)
    for j = 1:11
        if any(abs(D.events{j} - allEventTimes(i)) < 0.001)
            BinaryCode{i}=[BinaryCode{i};j];     
        end 
    end
    idigit=1;
    port=0;
    while idigit<=numel(BinaryCode{i})
        if BinaryCode{i}(idigit)<9
            EventCode(i)=2^(BinaryCode{i}(idigit)-1)+ EventCode(i);
        else
            EventCode(i)=2^(BinaryCode{i}(idigit)-1-8)+ EventCode(i);
            port=1;
        end
        idigit=idigit+1;
    end   
    if port==1
        EventCode(i)=EventCode(i)+128;
        port=0;
    end
end

% % %% Correct the trigger labels for Mikey
if imonkey==1
    idx=find(EventCode(1:end-1)==3);
    idx1=find(EventCode(idx+1)==132);
    EventCode(idx(idx1)+1)=4;

    EventCode(EventCode==134)=6;  % RESPONSE_SACCADE_STARTED = 6;
    EventCode(EventCode==133)=5; 

    idx=find(EventCode(1:end-1)==130);
    idx2=find(EventCode(idx-1)==129);
    EventCode(idx(idx2)-1)=1;
    EventCode(idx(idx2))=2;

    % idx2=find(ismember(EventCode(idx-1), [4 5]));
    % EventCode(idx(idx2))=130;

    idx=find(EventCode(1:end-1)==129);
    eventC=EventCode(idx-1);
    cc=unique(eventC);
    idx3=find(ismember(EventCode(idx-1), [1 5 6 128 129 130] ));
    EventCode(idx(idx3))=1;
end

%%%%% get the events time of interest
juiceEvent = allEventTimes(EventCode==128);
diffJuiceEvent = diff(juiceEvent);
idx=find(diffJuiceEvent<1);
FirstJuiceEvent=juiceEvent(idx)*1000;
RewardT=FirstJuiceEvent;

DelayT=allEventTimes(EventCode==3)*1000;
CueT = allEventTimes(ismember(EventCode,[8 16 32 64]))*1000;
allCue=EventCode(ismember(EventCode,[8 16 32 64]));
ArrayT=allEventTimes(EventCode==4)*1000;  % array onset
StartT=allEventTimes(EventCode==2)*1000;   % trial start
InitT=allEventTimes(EventCode==1)*1000;   % trial start

%%%% get the events time for only correct trials
MaxTrlT=6000;
num=0;
numArray=0;
TrlId=zeros(1,numel(RewardT));
load('Z:\RujiaChen\Results\p_arrayOn_trigger_correct.mat');

if exist('RewardId', 'var')>0
    clear RewardId;    
    clear Correct;
    clear SpikeTrain;
end
RewardId=[];

for itrl=1:numel(RewardT)
   idxDelay=find(DelayT < RewardT(itrl) & RewardT(itrl)-DelayT<MaxTrlT,1,'last');
   idxCue=find(CueT<RewardT(itrl)& RewardT(itrl) - CueT<MaxTrlT,1,'last');
   idxArray=find(ArrayT<RewardT(itrl)& RewardT(itrl) - ArrayT < MaxTrlT,1,'last');
   idxTrl=find(StartT<RewardT(itrl)& RewardT(itrl) - StartT < MaxTrlT,1,'last');
   idxInit=find(InitT<RewardT(itrl),1,'last');
   
   if ~isempty(idxDelay)&&~isempty(idxCue) &&~isempty(idxTrl)&&~isempty(idxInit)
%        num=num+1;
%        RewardId(num)=itrl;
       Correct.Cue(itrl)=allCue(idxCue);
       Correct.CueT(itrl)=CueT(idxCue);
       Correct.DelayT(itrl)=DelayT(idxDelay);
       Correct.ArrayT(itrl)=ArrayT(idxArray);  % get the array on events for Vasco
       Correct.ArrayOnsetT(itrl)=ArrayT(idxArray)+floor(idxInit*p(1)*60/1000)*1000/60;         
   end    
end
% Correct.RewardT=RewardT(RewardId);
Correct.RewardT=RewardT;
save(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat'],'Correct','-v7.3');

% %% check the duration of correct trials
% delayDu=Correct.ArrayT-Correct.DelayT;
% cueDu=Correct.DelayT-Correct.CueT;
% figure;
% subplot(1,2,1);
% hist(delayDu);
% xlabel('Delay duration');
% ylabel('Number of correct trials');
% subplot(1,2,2);
% [hh,pp]=hist(Correct.Cue,[8 16 32 64]);
% bar(1:4,hh);
% xlabel('Cue Type');
% ylabel('Number of correct trials');

% %% get the spike time of all correct trials 
NCorrect=numel(Correct.Cue);
% load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);

SpikeTrain.cueAlign = zeros(NCorrect, numel(D.allSpikeStructs),1801);
SpikeTrain.arrayAlign = zeros(NCorrect, numel(D.allSpikeStructs),1811);
SpikeTrain.rewardAlign = zeros(NCorrect, numel(D.allSpikeStructs),1501);

for ich=1:numel(D.allSpikeStructs)
    for itrl=1:NCorrect
        Duration=[Correct.CueT(itrl)-700  Correct.RewardT(itrl)+600];        
        idx00=D.allSpikeStructs{ich}.ts*1000>Duration(1)&D.allSpikeStructs{ich}.ts*1000<Duration(2);
        if sum(idx00)>0            
            spkTRaw=D.allSpikeStructs{ich}.ts(idx00)*1000;
            spkTcue=spkTRaw-Correct.CueT(itrl);
            BinT=-600:1200;   % cue delay=[600/900 1200]
            SpikeTrain.cueAlign(itrl,ich,:)=hist(spkTcue,BinT);

            if imonkey==1|| (imonkey==2 &&idate<=7)
                spkTarray=spkTRaw-Correct.ArrayOnsetT(itrl);
            elseif imonkey==2 && idate>7
                spkTarray=spkTRaw-Correct.ArrayT(itrl);   % after the trigger is corrected in the experimental code
            end
            BinT=-505:1305;  %hold delay=[500/600 900]
            SpikeTrain.arrayAlign(itrl,ich,:)=hist(spkTarray,BinT);
            
            spkTreward=spkTRaw-Correct.RewardT(itrl);
            BinT=-1000:500;
            SpikeTrain.rewardAlign(itrl, ich,:)=hist(spkTreward,BinT);
            
        end
    end 
    
    SpikeTrain.channel(ich)=D.allSpikeStructs{ich}.channelID;
    SpikeTrain.waveform(ich,:)=mean(D.allSpikeStructs{ich}.wf,1);            
end

%%%%  save the spike trains 
% save(['Z:\RujiaChen\Results\Unsorted_SpikeTrain_' areaName  '_' date '_5sigma_TriggerCorrected.mat'], 'SpikeTrain', '-v7.3' );
save(['Z:\RujiaChen\Results\Sorted_SpikeTrain_' areaName  '_' date '_TriggerCorrected.mat'],'SpikeTrain','-v7.3');

end

%% compare the waveform and select the units with best signals
folder='Z:\RujiaChen\Results\';
dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey 
clrN=[eye(3); 0 0 0; 0.5 0.5 0.5];
for idate=1:numel(dateUsed)
%     figure;
    date=dateUsed{idate};
    load([folder 'SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
    ElecID=unique(SpikeTrain.channel);
    Ncell=numel(ElecID);
    SpikeTrain.SU=zeros(1, numel(SpikeTrain.channel));
    for ielec=1:Ncell
%         subplot(ceil(sqrt(Ncell)), ceil(sqrt(Ncell)), ielec);
        idxCh=find(SpikeTrain.channel==ElecID(ielec));
        Wf=zeros(numel(idxCh), size(SpikeTrain.waveform,3));
        WfScaleRight=zeros(1, numel(idxCh));        
        WfScale=zeros(1, numel(idxCh));        
        for jj=1:numel(idxCh)
            Wf(jj,:)=squeeze(mean(SpikeTrain.waveform(:, idxCh(jj),:), 1));
%             plot(Wf(jj,:), 'color', clrN(jj,:)); hold on;
            [minV , minT]=min(Wf(jj,:));      
            WfScaleRight(jj)=max(Wf(jj,minT:end))-minV;
            WfScale(jj)=max(Wf(jj, :))-minV;
        end
        [maxR, idxR]=max(WfScaleRight);
        [maxAll, idxAll]=max(WfScale);
        if idxR==idxAll
            SpikeTrain.SU(idxCh(idxR))=1;
        elseif numel(idxCh)<=2
            SpikeTrain.SU(idxCh(idxAll))=2;
        else
            SpikeTrain.SU(idxCh(idxAll))=2;
            SpikeTrain.SU(idxCh(idxR))=1;
        end
    end
    save([folder 'SpikeTrain_allArea_' date '_TriggerCorrected.mat'], 'SpikeTrain', '-v7.3');
    
end


%% scatter plot of spikes in one channel
close all;
ElecID=unique(SpikeTrain.channel);
Time=-905:905;  %-300:1200;
plotRange=[-100 300];
for ielec=90:96 %:numel(ElecID)
    figure;
    idxCh=find(SpikeTrain.channel==ElecID(ielec));
    for jj=1:numel(idxCh)
        subplot(1, numel(idxCh), jj);
        data=squeeze(SpikeTrain.arrayAlign(:,idxCh(jj),:));
        % iTrl=find(CorrectTrialParam.rfyi==TargetPosY(2)&CorrectTrialParam.cueloc==cueLoc(2)&CorrectTrialParam.isexocue==1);
        % iTrl=find(Correct.Cue==16);
        for ii=1:size(data,1)  %numel(iTrl)  %size(data) %%% scatter plot for all trials in one condition
            mResp=squeeze(data(ii,:));
            idx1=find(Time>plotRange(1)&Time<plotRange(2)&mResp==1);
            plot(Time(idx1), ones(1,numel(idx1))*ii,'k.','linewidth', 1 ); hold on;
        end
        line([0 0],[0 ii+1],'linestyle','--','color','k');
    end
end

%% combine all the units in one channel and plot for each channel
clear all;
area='allArea';
date='122217';
load(['Z:\RujiaChen\Results\Unsorted_SpikeTrain_' area  '_' date '_TriggerCorrected.mat']);
% % load(['Z:\RujiaChen\Results\SpikeTrain_' area  '_' date '.mat']);
% load(['Z:\RujiaChen\Results\RFtrl_' date '.mat']);
% load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '.mat']);
load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat'])
%% plot the PSTH for an example channel
ChannelId=unique(SpikeTrain.channel);
cueTy=[8 16 32 64];
TargetPosition=[200 -200];
CuePsX=[200 -200; 200 -200];
CuePsY=[200 -200; -200 200];
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
isession=find(strcmp(RawInfo(1,:), date)==1);
% close all;
for iarea=3
    figure;
    clr=[1 0 0; 0 0 1; 1 0 0 ; 0 0 1];
    ifig=0;
    for ich=RawInfo{iarea+1, isession}   %32*(iarea-1)+1:32*iarea  %numel(ChannelId)
        idxCh=find(SpikeTrain.channel==ChannelId(ich));
        MUA.cue(ich,:,:)=squeeze(sum(SpikeTrain.cueAlign(:,idxCh,:),2));
        MUA.array(ich,:,:)=squeeze(sum(SpikeTrain.arrayAlign(:,idxCh,:),2));
        
        subplot(6,6, ich-32*(iarea-1));
        for ipos=1:2
            for icueTy=1:2
                %             ifig=ifig+1;
                %             subplot(2,2,ifig);
                idx0=find(Correct.Cue==cueTy(icueTy));  %&Correct.RF(:,2)'==TargetPosition(ipos));
                plotRange=[-250 550];
                Time=-300:1200;
                
                %             for ii=1:numel(idx0)  %%% scatter plot for all trials in one condition
                %                 mResp=squeeze(MUA.cue(ich,idx0(ii),:));
                %                 idx1=find(Time>plotRange(1)&Time<plotRange(2)&mResp'==1);
                %                 plot(Time(idx1), ones(1,numel(idx1))*ii,'.','color',clr(icueTy,:),'linewidth', 1 ,'DisplayName',num2str(cueTy)); hold on;
                %             end
                %             line([0 0],[0 ii+1],'linestyle','--','color','k');
                
%                 mResp=squeeze(mean(MUA.cue(ich,idx0,:),2));
                mResp=squeeze(mean(MUA.array(ich,idx0,:),2));
                %             mResp=mResp-mean(mResp(Time>-250&Time<0));
                mResp=conv(mResp,ones(1,50)/50,'same');
                idx1=Time>plotRange(1)&Time<plotRange(2);
                plot(Time(idx1), mResp(idx1),'color',clr(icueTy,:),'linewidth', 1.5 ,'DisplayName',num2str(cueTy)); hold on;
                
                title(['(' num2str(CuePsX(ipos,icueTy))  ', ' num2str(CuePsY(ipos,icueTy)) ')'])
                box off;
                xlabel('Time to cue onset');
                ylabel('Trial number');
                xlim([-250 550]);
%                 ylim([-0.02 0.08]);
            end
        end
    end
end








    
    
    







    