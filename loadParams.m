
%%%%% load the trial params for flanker sessions 
% clear all;
folder='Z:\RujiaChen\Results\';
% filename=dir([folder 'Unsorted_SpikeTrain_allArea_*_lowThreshold.mat']);
% dateUsed={'111717', '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey  ,'111917'
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};    % for Vasco
dateUsed={'112018', '112118', '010719', '010919', '011119', '011419', '011519', '011719', '011819', '012119'};  % for Vasco 
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118','110118','111618'};  % for Vasco
for ifile=1:numel(dateUsed)  %1:numel(filename)
%     dateID=strfind(filename(ifile).name,'_lowThreshold');
%     date=filename(ifile).name(dateID-6:dateID-1);
    date=dateUsed{ifile};
    logDate=['20' date(5:6) date(1:4)];
    logFolder=['Y:\anne\monkeys\presentation files\Vasco\logfiles\' logDate '\'];    % for Vasco
%     logFolder=['Y:\anne\monkeys\Mikey\logfiles\' logDate '\'];  % for mikey
    filepath=dir([logFolder 'flanker_trial_params*.txt']);
    filepath2=dir([logFolder 'flanker_trial_timing*.txt']);
    sizeFile=zeros(numel(filepath2),1);
    for ipath=1:numel(filepath2)
        sizeFile(ipath)=filepath2(ipath).bytes;
    end
    if numel(filepath2)>1
        [~, idxFile]=max(sizeFile);
    else
        idxFile=1;
    end
    if ~isempty(filepath)
        logPresentation=[logFolder filepath2(idxFile).name];
        fileID=fopen(logPresentation,'r');
        formatSpec='%d %d %d';
        dataRaw=textscan(fileID,formatSpec);
        ValidTrl=dataRaw{1}(dataRaw{2}==2);
        ArrayTPrest=dataRaw{3}(dataRaw{2}==4);
        save(['Z:\RujiaChen\Results\TrialTiming_' date '.mat'],'dataRaw','-v7.3');
        
        
        logPresentation=[logFolder filepath(idxFile).name];
        fileID=fopen(logPresentation,'r');
        tline=fgets(fileID);
        headInfo=textscan(tline,'%s','delimiter',' ');
        
        formatSpec='%f %f %f %f %f %f %f %f %f %f %f %s %s %s %f %f %s %f %f %f %f %f %f %f %f';
        startRow=1;
        
%         dataRaw=textscan(fileID,formatSpec);
        TrlInfo=textscan(fileID, formatSpec,  'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for ii=[12:14 17]
            TrlInfo{ii}=cellfun(@(x) strrep(x,'false','0'), TrlInfo{ii}, 'UniformOutput', false);
            TrlInfo{ii}=cellfun(@(x) strrep(x,'true','1'), TrlInfo{ii}, 'UniformOutput', false);
            TrlInfo{ii}=cellfun(@str2double, TrlInfo{ii});
        end
        ntrl=cellfun(@numel,TrlInfo);
        if unique(ntrl)>1
            validTrl=min(ntrl);
            TrlInfo=cellfun(@(x) x(1:validTrl),TrlInfo, 'UniformOutput', false);
        end  
        allInfo=cell2mat(TrlInfo);
%         TrlNum=unique(allInfo(:,1));
        TrlNum=ValidTrl;
        ValidInfo=zeros(numel(TrlNum),numel(TrlInfo));
        num=0;
        for ii=1:numel(TrlNum)
            num=num+1;
            idx00=find(allInfo(:,1)==ii,1,'last');
            ValidInfo(ii,:)=allInfo(idx00,:);    
        end
        for icell=1:numel(TrlInfo)
            TrlParam.(headInfo{1}{icell})=ValidInfo(:,icell);
        end
    end 
    save(['Z:\RujiaChen\Results\flanker_TrlParam_' date '.mat'],'TrlParam','-v7.3');
    
    idx=TrlParam.trial_response==-105;
    for icell=1:numel(headInfo{1})
        CorrectTrialParam.(headInfo{1}{icell})=TrlParam.(headInfo{1}{icell})(idx);   %CorrectTrlID(1:num)
    end
%     save(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat'],'CorrectTrialParam','-v7.3');
end

%% compare the reaction time in different conditions
close all;
folder='Z:\RujiaChen\Results\';
% dateUsed={'111717', '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey  ,'111917'
dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718','092618','092718','092818'};  % for Vasco
figure;
for idate=8 :numel(dateUsed)
%     figure;
    subplot(2,3,idate-7);
    date=dateUsed{idate};
    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
    idx1=find(CorrectTrialParam.is_mouse_trial==0&CorrectTrialParam.isexocue==1);
    RTime{1}=CorrectTrialParam.RT(idx1);
    idx2=find(CorrectTrialParam.is_mouse_trial==0&CorrectTrialParam.isexocue==0);
    RTime{2}=CorrectTrialParam.RT(idx2);
%     for ii=1:2
%         MeanRT(idate-10, ii)=mean(RTime{ii});        
%         subplot(1,2,ii);
%         hist(RTime{ii});
%     end
    plot(idx1, RTime{1}, 'r*'); hold on;
    plot(idx2, RTime{2}, 'b*'); hold on;
    xlim([0 400]);
    ylim([0 1000])
    
end

%%
% clear all;
% folder='Z:\RujiaChen\Results\';
% % filename=dir([folder 'Unsorted_SpikeTrain_allArea_*_lowThreshold.mat']);
% dateUsed={'111717', '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey  ,'111917'
% for ifile=1:numel(dateUsed)
% %     load([folder filename(ifile).name]) ;
% %     dateID=strfind(filename(ifile).name,'_lowThreshold');
% %     date=filename(ifile).name(dateID-6:dateID-1);
%     date=dateUsed{ifile};
%     load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat']);
%     
% %     folderData='Y:\anne\monkeys\plexondata\';   % for Vasco
%     folderData='Y:\anne\monkeys\pigz_plexondata\mikey_left_hemisphere\';  % for Mikey
%     fileName=[folderData date '\merged_spl_noWB_001-01.pl2'];
%     muaDataDirRoot=[folderData date '\'];   %042118
%     sessionName='merged_spl_noWB_001-01';
%     areaName='allArea';
%    
%     isLoadSpikes=1;
%     isLoadMua=0;
%     isLoadLfp=0;
%     isLoadSpkc=0;
%     isLoadDirect=0;
% %     spikeChannelPrefix = 'SPK_SPKC';
%     spikeChannelPrefix = 'SPK';
%     channel=33;
%     spikeChannelsToLoad=channel;
%     muaChannelsToLoad=channel;
%     lfpChannelsToLoad = channel;
%     spkcChannelsToLoad = channel;
%     directChannelsToLoad = channel;
%     D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
%         spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);
%     
%     duration = D.blockStopTimes-D.blockStartTimes;
%     [~,gratingsTask3DIndices]=max(duration);
%     D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);
%  
%     allEventTimes = cell2mat(D.events([1:11])');
%     allEventTimes = unique(allEventTimes);
%     BinaryCode=cell(size(allEventTimes));
%     EventCode=zeros(size(allEventTimes));
%     for i = 1:numel(allEventTimes)
%         for j = 1:11
%             if any(abs(D.events{j} - allEventTimes(i)) < 0.001)
%                 BinaryCode{i}=[BinaryCode{i};j];
%             end
%         end
%         idigit=1;
%         port=0;
%         while idigit<=numel(BinaryCode{i})
%             if BinaryCode{i}(idigit)<9
%                 EventCode(i)=2^(BinaryCode{i}(idigit)-1)+ EventCode(i);
%             else
%                 EventCode(i)=2^(BinaryCode{i}(idigit)-1-8)+ EventCode(i);
%                 port=1;
%             end
%             idigit=idigit+1;
%         end
%         if port==1
%             EventCode(i)=EventCode(i)+128;
%             port=0;
%         end
%     end
% 
%     %%%%%%% Correct the trigger labels for Mikey
%     idx=find(EventCode(1:end-1)==3);
%     idx1=find(EventCode(idx+1)==132);
%     EventCode(idx(idx1)+1)=4;
%     
%     EventCode(EventCode==134)=6;  % RESPONSE_SACCADE_STARTED = 6;
%     EventCode(EventCode==133)=5;
%     
%     idx=find(EventCode(1:end)==130);
%     idx2=find(EventCode(idx-1)==129);
%     EventCode(idx(idx2)-1)=1;
%     EventCode(idx(idx2))=2;
%     
% %     idx2=find(ismember(EventCode(idx-1), [4 5]));
% %     EventCode(idx(idx2))=130;
%     
%     idx=find(EventCode(1:end-1)==129);
%     eventC=EventCode(idx-1);
%     cc=unique(eventC);
%     idx3=find(ismember(EventCode(idx-1), [1 5  6 128 129 130] ));
%     EventCode(idx(idx3))=1;
%         
%     %%%%% get the events time of interest
%     juiceEvent = allEventTimes(EventCode==128);
%     diffJuiceEvent = diff(juiceEvent);
%     idx=find(diffJuiceEvent<1);
%     FirstJuiceEvent=juiceEvent(idx)*1000;
%     RewardT=FirstJuiceEvent;
%     
%     DelayT=allEventTimes(EventCode==3)*1000;
%     CueT = allEventTimes(ismember(EventCode,[8 16 32 64]))*1000;
%     allCue=EventCode(ismember(EventCode,[8 16 32 64]));
%     ArrayT=allEventTimes(EventCode==4)*1000;  % array onset
%     StartT=allEventTimes(EventCode==2)*1000;   % trial start
%     
%     load(['Z:\RujiaChen\Results\flanker_TrlParam_' date '.mat']);
% 
%     %%%% get the events time for only correct trials
%     MaxTrlT=5000;
%     num=0;
%     CorrectTrlID=zeros(1,numel(RewardT));
%     
%     if numel(TrlParam.ntr)~=numel(StartT)
% %         error('The number of valid trials in plexon and logfile is not same!');
%         StartT=StartT(1:numel(TrlParam.ntr));
%     end
%     
%     for itrl=1:numel(RewardT)
%         idxDelay=find(DelayT < RewardT(itrl) & RewardT(itrl)-DelayT<MaxTrlT,1,'last');
%         idxCue=find(CueT<RewardT(itrl)& RewardT(itrl) - CueT<MaxTrlT,1,'last' );
%         idxArray=find(ArrayT<RewardT(itrl)& RewardT(itrl) - ArrayT < MaxTrlT,1,'last' );
%         idxTrl=find(StartT<RewardT(itrl)& RewardT(itrl) - StartT < MaxTrlT,1,'last' );
%         if ~isempty(idxDelay)&&~isempty(idxCue) &&~isempty(idxTrl)
%             num=num+1;
%             CorrectTrlID(num)=idxTrl;
% %             Correct.Cue(num)=allCue(idxCue);
% %             Correct.CueT(num)=CueT(idxCue);
% %             Correct.DelayT(num)=DelayT(idxDelay);
% %             Correct.ArrayT(num)=ArrayT(idxArray);  % get the array on events for Vasco
% 
%         end
%     end
% end


