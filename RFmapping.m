clear all; 
close all;
% folder='Z:\RujiaChen\Vasco\';
monkeyset={'Mikey','Vasco'};
imonkey=2;
monkey=monkeyset{imonkey};
if imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'}; 
    folder='Y:\anne\monkeys\plexondata\';  %for vasco
    channel=1:96;  % start with pulvinar
elseif imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
    folder='Y:\anne\monkeys\pigz_plexondata\mikey_left_hemisphere\';  % for Mikey
    load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
    isession=find(strcmp(RawInfo(1,:), date)==1);
    channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
end

for idate=1:numel(dateUsed)
date=dateUsed{idate};
fileName=[folder date '\merged_spl_noWB_001-01.pl2'];
muaDataDirRoot=['Z:\RujiaChen\Vasco\' date '\'];   %042118
sessionName='merged_spl_noWB_001-01';

areaName='allArea';   %'PUL';  %
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
directChannelsToLoad = channel;
D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);
    
% % %% trim only the rf mapping data
duration=D.blockStopTimes-D.blockStartTimes;
[~,ii]=max(duration);

%%%% get the params in he pre mapping session
RFtask=ii-1;
% while duration(RFtask)>1500
%     RFtask=RFtask-1;
% end
sessionOrder='pre';
D = trimSpikeTimesAndEvents(D, RFtask);

%%%% get the params in the post mapping session
% RFtask2=ii+1;
% if RFtask2<=numel(duration)&&duration(RFtask2)<0.5
%     sessionOrder='post';
%     D = trimSpikeTimesAndEvents(D, RFtask); 
% end

% %% read the events times in plexon
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

%%%%% get all the flash time 
CueOnT=allEventTimes(EventCode==16)*1000;
FlashOnT=allEventTimes(EventCode==32)*1000;

save(['Z:\RujiaChen\Results\RFmapping_' sessionOrder '_' date '_CueOnT.mat'],'CueOnT','-v7.3');
save(['Z:\RujiaChen\Results\RFmapping_' sessionOrder '_' date '_FlashOnT.mat'],'FlashOnT','-v7.3');


% % %% get the spike trains of all trials (16 flash/trl)
clear SpikeTrain;
Ntrl=length(CueOnT);   %length(FlashOnT);   %
% LFPresp=zeros(Ntrl, size(D.adjLfps,1), 2801);
% Resp=zeros(Ntrl,numel(D.allSpikeStructs));
for itrl=1:Ntrl
    Duration=[CueOnT(itrl)-300 CueOnT(itrl)+2600];  %% the param is recorded by trials
%     Duration=[FlashOnT(itrl)-100 FlashOnT(itrl)+200];
    
%     lfpIndicesToKeep1 = false(1, size(D.adjLfps, 2));
%     blockLfpIndices=max(1,floor(Duration(1))):floor(Duration(2));
%     lfpIndicesToKeep1(blockLfpIndices) = true;
%     LFPresp(itrl,:,:)=D.adjLfps(:,lfpIndicesToKeep1);    
    
    for ich=1:numel(D.allSpikeStructs)
        idx00=find(D.allSpikeStructs{ich}.ts*1000>Duration(1)& D.allSpikeStructs{ich}.ts*1000<Duration(2));
        if ~isempty(idx00)
            spkTime = D.allSpikeStructs{ich}.ts(idx00)*1000-CueOnT(itrl);
            BinT=-300:5:2600; % 16*150 for first session OR 8*300 for second session TO VASCO
            SpikeTrain(itrl,ich,:)=hist(spkTime,BinT);
%             Resp(itrl,ich)=numel(D.allSpikeStructs{ich}.ts(idx00)*1000);
        end
    end 
end
 
save(['Z:\RujiaChen\Results\SpikeTrain_RFmapping_' sessionOrder '_' areaName  '_' date '.mat'],'SpikeTrain','-v7.3');
% save(['Z:\RujiaChen\Results\LFPresp_RFmapping_' sessionOrder '_' areaName  '_' date '.mat'],'LFPresp','-v7.3');

end
%% plot the time course of visual response across all trials 
% close all;
TT=-100:200;
for iarea=2
    figure;
    ifig=0;
    for ich=(iarea-1)*32+1:32*iarea
        mTrain=squeeze(mean(SpikeTrain(:,ich,:),1)) ;
%         mTrain=squeeze(mean(LFPresp(:,ich,:),1)) ;
        mTrain=mTrain-mean(mTrain(1:100));
        ifig=ifig+1;
        subplot(6,6,ifig);
        mTrain=conv(mTrain,ones(1,30)/30,'same');
        plot(TT,mTrain); hold on;
        axis tight;
    end
end

%%
clear all; 
close all;
% folder='Z:\RujiaChen\Vasco\';
monkeyset={'Mikey','Vasco'};
imonkey=2;
monkey=monkeyset{imonkey};
if imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'}; 
    channel=1:96; 
    logFolder='Y:\anne\monkeys\presentation files\Vasco\logfiles\';  %for vasco
elseif imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
    load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
    isession=find(strcmp(RawInfo(1,:), date)==1);
    channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % for Mikey
    logFolder='Y:\anne\monkeys\Mikey\logfiles\' ;
end
areaName='allArea';
sessionOrder='pre';
Time=-300:5:2600;  % 100(200) ms onset duration, 50(100) ms interval, spikes binned in 5-ms steps
folder='Z:\RujiaChen\Results\';
for idate=8 %8:numel(dateUsed)
    date=dateUsed{idate};
    load([folder 'SpikeTrain_RFmapping_' sessionOrder '_' areaName  '_' date '.mat']);
    logDate=['20' date(5:6) date(1:4)];
    logPath=[logFolder logDate '\'];
    filename=dir([logPath 'multi_flash_rf_mapping_results_*.json']);
    
    % load(['Z:\RujiaChen\Results\LFPresp_RFmapping_' sessionOrder '_' areaName  '_' date '.mat']);
    load([folder 'RFmapping_' sessionOrder '_' date '_CueOnT.mat']);
    load([folder 'RFmapping_' sessionOrder '_' date '_FlashOnT.mat']);
    
    if numel(FlashOnT)/numel(CueOnT)>8
        NflashPerTrl=16;        
    else
        NflashPerTrl=8;
%         Time=-300:1300;
    end
    
    if ~isempty(filename)
        if strcmp(sessionOrder,'pre')==1
            iid=max([numel(filename)-1, 1]);
        elseif strcmp(sessionOrder,'post')==1
            iid=numel(filename);
        end
        
        trialStructs = loadjson(sprintf('%s/%s', logPath, filename(iid).name)); % for now, just the first file
        clear flashParams;
        assert(numel(trialStructs)==size(SpikeTrain,1),'The number of trials is different between logfile and plexon!');   %LFPresp
        itrl=0;
        flashCount=0;
        for ii = 1:numel(trialStructs)
            if ii<numel(trialStructs)
                flashT=ceil(FlashOnT(FlashOnT<CueOnT(ii+1)&FlashOnT>CueOnT(ii))-CueOnT(ii));
            else
                flashT=ceil(FlashOnT(FlashOnT>CueOnT(ii))-CueOnT(ii));
            end
            if numel(flashT)==NflashPerTrl %trialStructs{ii}.numFlashesShownTotal     %(numel(trialStructs{ii}.events)-1)
                itrl=itrl+1;
                for jj = 2:numel(trialStructs{ii}.events)
                    flashCount = flashCount + 1;
                    flashParams(flashCount) = trialStructs{ii}.events{jj}.flashParams;
                    idxT=find(Time>flashT(jj-1)+50&Time<flashT(jj-1)+200);
                    idx0=find(Time>-200&Time<0);
                    flashResp(flashCount,:)=squeeze(mean(SpikeTrain(ii,:,idxT),3))-squeeze(mean(SpikeTrain(ii,:,idx0),3));
                    %                 flashResp(flashCount,:)=squeeze(mean(LFPresp(ii,:,idxT),3))-squeeze(mean(LFPresp(ii,:,idx0),3));
                end
            end
        end
    end
    
    mm=struct2cell(flashParams);
    FlashInfo=squeeze(cell2mat(mm));
    stimID=unique(FlashInfo(1,:));
    Nstim=numel(stimID);
    mResp=zeros(Nstim,size(SpikeTrain,2));
    for istim=1:numel(stimID)
        idx=find(FlashInfo(1,:)==stimID(istim));
        mResp(istim,:)= mean(flashResp(idx,:),1);
        distStim(istim)=FlashInfo(2,idx(1));
        angleStim(istim)=FlashInfo(3,idx(1));
        diamStim(istim)=FlashInfo(4,idx(1));
    end
     
    center=[200,-200,200,-200; 200, -200,-200, 200];
    Rcircle=100;
     % %% get the distance of each stimulus in RF mapping to cue/target location
     stimLoc=zeros(2,Nstim);
    for istim=1:Nstim
        stimLoc(1, istim)=distStim(istim)*cos(angleStim(istim));
        stimLoc(2, istim)=distStim(istim)*sin(angleStim(istim));
    end
    
    Dist2Target=zeros(4,Nstim);
    cueLocResp=zeros(4,96);
    for ii=1:4
        Dist2Target(ii, :)=sqrt((stimLoc(1,:)-center(1, ii)).^2+(stimLoc(2,:)-center(2, ii)).^2);
        idx=Dist2Target(ii, :)<Rcircle;
        if ~isempty(idx)
            cueLocResp(ii,:)=mean(mResp(idx,:),1);
        end
    end
    %%%% plot the mean response to stimulus near to the cue/target location
%     load([folder 'CueOn_' date]);
%     for iarea=3
%         figure;
%         subplot(1,2,1);
%         NormResp=cueLocResp(:, 32*(iarea-1)+1:32*iarea)'./transpose(repmat(max(cueLocResp(:, 32*(iarea-1)+1:32*iarea)),4,1));
% %         imagesc(NormResp); hold on;
%         imagesc(cueLocResp(:, 32*(iarea-1)+1:32*iarea)'); hold on;
% %         caxis([0 1]);
%         colorbar;
%         axis xy;
%         iidx=find(CueOn(32*(iarea-1)+1:32*iarea,1)>0);
%         plot(CueOn(iidx+32*(iarea-1),1), iidx, '*b'); hold on;
%         iidx=find(CueOn(32*(iarea-1)+1:32*iarea,2)>0);
%         plot(CueOn(iidx+32*(iarea-1),2)+2, iidx, '*r'); hold on;
%         subplot(1,2,2);
% %         [rr, pp]=corr(cueLocResp(:, 32*(iarea-1)+1:32*iarea),'type','spearman');
% %         imagesc(rr);
%         [maxVal, pp]=max(cueLocResp(:, 32*(iarea-1)+1:32*iarea),[],1);
%         plot(pp, 1:32) ; hold on;
%     end    
     % %%  plot the response mapping for each channel
    clr=colormap(jet(Nstim));
    for iarea=1:3
        figure;
        nfig=0;
        for ich=32*(iarea-1)+1:32*iarea
            mm=mResp(:,ich);
%             mm=max(mm,0);
            mm=mm-min(mm);
            mm=mm./max(mm);
            nfig=nfig+1;
            
            if sum(mm)>0
                subplot(6,6,nfig);
                for ist=1:Nstim
                    if mm(ist)>0
%                         polarplot(angleStim(ist),distStim(ist),'o','markerfacecolor','b','markeredgecolor','b','markersize', 8*mm(ist)); hold on;
                        plot( stimLoc(1, ist),  stimLoc(2, ist), 'o','markerfacecolor','b','markeredgecolor','b','markersize', 8*mm(ist)); hold on;
                    end
                end
                
                viscircles(center',ones(4,1)*Rcircle,'color', 'r'); hold on;
                title(['ch' num2str(ich)]);
                axis equal;
            end
        end
    end
    save([folder 'meanResp_cueLoc_RFmap_' date '.mat'], 'cueLocResp','-v7.3');
end

%%
close all;
folder='Z:\RujiaChen\Results\';
dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};
for idate=1:numel(dateUsed)
    figure;
    date=dateUsed{idate};
    load([folder 'meanResp_cueLoc_RFmap_' date '.mat']);  % cueLocResp: response to gratings in RF mapping task
    load([folder 'Resp2cue_' date '.mat']);  % cueMeanResp:  response to cue stimuli at the flanker task
    load([folder 'CueOn_' date '.mat']);
    load([folder 'bArrayResp_' date '.mat']);
    subplot(1,2,1);
    idx=sum(bArrayResp,2)>0;
    data1=zeros(size(cueMeanResp));
    data1(idx,:)=cueMeanResp(idx,:);
    imagesc(data1);
    [~,maxIdx1]=max(cueMeanResp,[],2); hold on;
    yy=1:numel(maxIdx1);
    plot(maxIdx1(idx),yy(idx),'r*'); hold on;
    title('flanker task');
    colorbar;
    subplot(1,2,2);
    data2=zeros(size(cueLocResp));
    data2(:,idx)=cueLocResp(:, idx);
    imagesc(data2');
    colorbar;
    [~,maxIdx2]=max(cueLocResp,[],1); hold on;
    yy=1:numel(maxIdx2);
    plot(maxIdx2(idx),yy(idx),'r*'); hold on;
    title('RF mapping task'); 
    
    idx00=maxIdx1==maxIdx2'&sum(bArrayResp,2)>0;
    idx11=zeros(96,1);
    for ii=1:96
        idx11(ii)=bArrayResp(ii, ceil(maxIdx1(ii)/2))>0; 
    end
    plot(maxIdx2(idx00&idx11),yy(idx00&idx11),'y*'); hold on;
    MaxLoc=zeros(96,1);
    MaxLoc(idx00&idx11)=maxIdx1(idx00&idx11);
    
    save([folder 'Estimated_RFLoc_' date '.mat'],'MaxLoc','-v7.3'); 
end






