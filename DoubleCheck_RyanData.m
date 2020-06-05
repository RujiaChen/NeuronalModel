clear all;
folder='Y:\ryanly\McCartney\merged\M20170311\';
eventArray=32;
eventCue=8;

areaName='LIP';

fileName=[folder '\M20170311-all_merged_noWB-d1_35.00mm_d2_35.00mm.pl2'];
muaDataDirRoot=folder;   %042118
sessionName='M20170311-all_merged_noWB-d1_35.00mm_d2_35.00mm';
areaName='PUL';
%%%%% blocks=6,7,9,10,11,12;

isLoadSpikes=1;
isLoadMua=0;
isLoadLfp=0;
isLoadSpkc=0;
isLoadDirect=0;
% spikeChannelPrefix = 'SPK_SPKC';
spikeChannelPrefix = 'SPK';
channel=1:32;
spikeChannelsToLoad=channel;
muaChannelsToLoad=channel;
lfpChannelsToLoad = channel;
spkcChannelsToLoad = channel;
directChannelsToLoad = channel;
D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);
    
    
%%
gratingsTask3DIndices=6;
D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);

%% read the events times in plexon
allEventTimes = cell2mat(D.events([1:8])');
allEventTimes = unique(allEventTimes);
BinaryCode=cell(size(allEventTimes));
EventCode=zeros(size(allEventTimes));
for i = 1:numel(allEventTimes)
    for j = 1:8
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

%%%%% get the events time of interest
juiceEvent = allEventTimes(EventCode==128);
diffJuiceEvent = diff(juiceEvent);
idx=find(diffJuiceEvent>1);
num=0;
mm=cat(2,diffJuiceEvent(idx-3),diffJuiceEvent(idx-2),diffJuiceEvent(idx-1));
nn=sum(mm<1,2);
FirstJuiceEvent=juiceEvent(idx(nn==3)-3)*1000;
% 
RewardT=FirstJuiceEvent;
% DelayT=allEventTimes(EventCode==3)*1000;
CueT = allEventTimes(ismember(EventCode,[8 10]))*1000;
allCue=EventCode(ismember(EventCode,[8 10]));
ArrayT=allEventTimes(EventCode==32)*1000;  % array onset

%%%% get the event code of correct trials
MaxTrlT=5000;
num=0;
numArray=0;
TrlId=zeros(1,numel(RewardT));
for itrl=1:numel(RewardT)
%    idxDelay=find(DelayT < RewardT(itrl) & RewardT(itrl)-DelayT<MaxTrlT,1,'last');
   idxCue=find(CueT<RewardT(itrl)& RewardT(itrl) - CueT<MaxTrlT,1,'last' );
   idxArray=find(ArrayT<RewardT(itrl)& RewardT(itrl) - ArrayT < MaxTrlT,1,'last' );
%    idxTrl=find(StartT<RewardT(itrl)& RewardT(itrl) - StartT < MaxTrlT,1,'last' );
   if ~isempty(idxCue) &&~isempty(idxArray)
       num=num+1;
       RewardId(num)=itrl;
       Correct.Cue(num)=allCue(idxCue);
       Correct.CueT(num)=CueT(idxCue);
       Correct.ArrayT(num)=ArrayT(idxArray);  % get the array on events for Vasco
   end
end
Correct.RewardT=RewardT(RewardId);

%%
NCorrect=numel(Correct.Cue);
clear SpikeTrain;
% load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);

for itrl=1:NCorrect
    nCh=0;
    Duration=[Correct.CueT(itrl)-300  Correct.RewardT(itrl)];
    for ich=1:numel(D.allSpikeStructs)
        idx00=D.allSpikeStructs{ich}.ts*1000>Duration(1)&D.allSpikeStructs{ich}.ts*1000<Duration(2);
        if ~isempty(idx00)
            nCh=nCh+1;
            spkTRaw=D.allSpikeStructs{ich}.ts(idx00)*1000;
            spkTcue=spkTRaw-Correct.CueT(itrl);
            BinT=-300:602;
            SpikeTrain.cueAlign(itrl,nCh,:)=hist(spkTcue,BinT);

%             spkTarray=spkTRaw-Correct.DelayT(itrl)-CorrectTrialParam.delay_duration(itrl);
            spkTarray=spkTRaw-Correct.ArrayT(itrl);
            BinT=-300:602;
            SpikeTrain.arrayAlign(itrl,nCh,:)=hist(spkTarray,BinT);
            
            spkTreward=spkTRaw-Correct.RewardT(itrl);
            BinT=-1000:2;
            SpikeTrain.rewardAlign(itrl,nCh,:)=hist(spkTreward,BinT);
            SpikeTrain.channel(nCh)=D.allSpikeStructs{ich}.channelID;
        end
    end 
end

%% scatter plot of spikes in one channel
figure;
Time=-300:602;
plotRange=[-200 500];
for idxCh=16  %1:36 %numel(SpikeTrain.channel)
%     subplot(6,6,idxCh);
    % idxCh=find(SpikeTrain.channel==77);
    data=squeeze(SpikeTrain.arrayAlign(:,idxCh,:));
    % iTrl=find(CorrectTrialParam.rfyi==TargetPosY(2)&CorrectTrialParam.cueloc==cueLoc(2)&CorrectTrialParam.isexocue==1);
    for ii=1:size(data) %%% scatter plot for all trials in one condition
        mResp=squeeze(data(ii,:));
        idx1=find(Time>plotRange(1)&Time<plotRange(2)&mResp==1);
        plot(Time(idx1), ones(1,numel(idx1))*ii,'.','linewidth', 1 ); hold on;
    end
    line([0 0],[0 ii+1],'linestyle','--','color','k');

end













