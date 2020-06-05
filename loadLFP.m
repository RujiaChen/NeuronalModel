% fileName='C:\Rujia\Mikey\Plexon\RawData\111717\merged_spl_noWB_001-02.pl2';
clear all; 
close all;
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
end
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118','110118','111618'};  % for Vasco
areaName='allArea';

%%
for idate=1:7   %numel(dateUsed)
    date=dateUsed{idate};
    if imonkey==1
        isession=find(strcmp(RawInfo(1,:), date)==1);
        channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}]; 
    end
    
%     if ~exist(['Z:\RujiaChen\Results\LFP_' areaName  '_' date '_TriggerCorrected.mat'], 'file')
    fileName=[folder date '\merged_spl_noWB_001-01.pl2'];
    muaDataDirRoot='Y:\anne\monkeys\plexondata\';  % for vasco   ['Z:\RujiaChen\Vasco\' date '\'];  
    sessionName='merged_spl_noWB_001-01';
    isLoadSpikes=0;
    isLoadMua=0;
    isLoadLfp=1;
    isLoadSpkc=0;
    isLoadDirect=0;
    spikeChannelPrefix = 'SPK_SPKC';  % 'SPK'
    spikeChannelsToLoad=channel;
    muaChannelsToLoad=channel;
    lfpChannelsToLoad = channel; 
    spkcChannelsToLoad = channel;
    directChannelsToLoad = 3;
    D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
            spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);

%         duration=D.blockStopTimes-D.blockStartTimes;
%         [~,gratingsTask3DIndices]=max(duration);
%         D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);

%     %% get the LFP of all correct trials 
    load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat']);
    NCorrect=numel(Correct.Cue);
    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);    
    clear LFPinfo;
    for itrl=1:NCorrect
        lfpIndicesToKeep1 = false(1, size(D.adjLfps, 2));
        lfpIndicesToKeep2 = false(1, size(D.adjLfps, 2));
        lfpIndicesToKeep3 = false(1, size(D.adjLfps, 2));
        
        if Correct.CueT(itrl)~=0
            blockLfpIndices=max(1,floor(Correct.CueT(itrl)-600)):floor(Correct.CueT(itrl)+1200);
            lfpIndicesToKeep1(blockLfpIndices) = true;
            LFPinfo.cueAlign(itrl,:,:)=D.adjLfps(:,lfpIndicesToKeep1);
            
            blockLfpIndices=max(1,floor(Correct.ArrayOnsetT(itrl)-900)):floor(Correct.ArrayOnsetT(itrl)+900); %for data before 092718
%             blockLfpIndices=max(1,floor(Correct.ArrayT(itrl)-900)):floor(Correct.ArrayT(itrl)+900);
            lfpIndicesToKeep2(blockLfpIndices) = true;
            LFPinfo.ArrayAlign(itrl,:,:)=D.adjLfps(:,lfpIndicesToKeep2);
            
            blockLfpIndices=max(1,floor(Correct.RewardT(itrl)-1000)):floor(Correct.RewardT(itrl)+500);
            lfpIndicesToKeep3(blockLfpIndices) = true;
            LFPinfo.RewardAlign(itrl,:,:)=D.adjLfps(:,lfpIndicesToKeep3);
        end
    end    
    save(['Z:\RujiaChen\Results\LFP_' areaName  '_' date '_TriggerCorrected.mat'],'LFPinfo','-v7.3');  
end

%% get the LFP signals of interest first, then calculate the mean LFP 
close all;
clear all;
monkeyset={'Mikey','Vasco'};
imonkey=1;
monkey=monkeyset{imonkey};
if imonkey==2
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'}; 
    channel=1:96;  
elseif imonkey==1
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
    load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');    
end
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118','110118','111618'};  % for Vasco
folder='Z:\RujiaChen\Results\';

TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
cuetype={'endo','exo'};

for igroup=1:3
    for isexo=0:1
        count=zeros(1,3);
        CueOnLFP=cell(1,3);
        CueOffLFP=cell(1,3);
        ArrayOnLFP=cell(1,3);
        ArrayOffLFP=cell(1,3);
        
        for idate=1:numel(dateUsed)
            date=dateUsed{idate};
            load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
%             load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);  %
            load([folder 'CorrectTrialParam_' date '.mat']);
            load([folder 'bArrayResp_' date '.mat']);
            load([folder 'bVisResp_' date '.mat']);
            load([folder 'bMotorResp_' date '.mat']);
            load([folder 'bVisMotor_' date '.mat']);
            load([folder 'CueOn_' date]);
            
            for iarea=1:3
                if imonkey==1
                    isession=find(strcmp(RawInfo(1,:), date)==1);                    
                    channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
                    channelID=find(ismember(channel, RawInfo{iarea+1,isession}));  % for Mikey
                else
                    channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
                end
                
                if ~isempty(channelID)
                    for ich=channelID
                        for ipos=1:2
                            if igroup==1
                                celltype=bVisResp(ich, ipos);
                            elseif igroup==2
                                celltype=bVisMotor(ich, ipos);
                            elseif igroup==3
                                celltype=bMotorResp(ich, ipos);
                            end
                            if celltype>0
                                count(iarea)=count(iarea)+1;
                                if CueOn(ich,ipos)>0
                                    icue=CueOn(ich,ipos);
                                else
                                    icue=imonkey;  % 1 for Mikey; 2 for Vasco
                                end
                                
                                idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                tt=squeeze(LFPinfo.cueAlign(idx1,ich,:));
                                CueOnLFP{iarea}(count(iarea),:)=mean(tt,1)-repmat(mean(tt,2),1,size(tt,2));
                                tt=squeeze(LFPinfo.ArrayAlign(idx1,ich,:));
                                ArrayOnLFP{iarea}(count(iarea),:)=mean(tt,1)-repmat(mean(tt,2),1,size(tt,2));
                                
                                idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                tt=squeeze(LFPinfo.cueAlign(idx2,ich,:));
                                CueOffLFP{iarea}(count(iarea),:)=mean(tt,1)-repmat(mean(tt,2),1,size(tt,2));
                                tt=squeeze(LFPinfo.ArrayAlign(idx2,ich,:));
                                ArrayOffLFP{iarea}(count(iarea),:)=mean(tt,1)-repmat(mean(tt,2),1,size(tt,2));
                                
                            end
                        end
                    end
                end
            end
        end
        
        save([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_CueOnResp_deMeaned.mat'],'CueOnLFP','-v7.3');
        save([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_CueOffResp_deMeaned.mat'],'CueOffLFP','-v7.3');
        save([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_ArrayOnResp_deMeaned.mat'],'ArrayOnLFP','-v7.3');
        save([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_ArrayOffResp_deMeaned.mat'],'ArrayOffLFP','-v7.3');
    end
end

%% plot the mean and std of population LFP 
% close all;
monkeyset={'Mikey','Vasco'};
igroup=2;
TimeCue=-600:1200;
TimeArray=-900:900;
for imonkey=1
    figure;
    monkey=monkeyset{imonkey};
    for iarea=1:3
        %     figure;
        for isexo=0:1
            load([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_CueOnResp_deMeaned.mat']);  %
            load([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_CueOffResp_deMeaned.mat']);
            load([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_ArrayOnResp_deMeaned.mat']);
            load([folder monkey '_LFP_' cuetype{isexo+1} '_' Group{igroup} '_ArrayOffResp_deMeaned.mat']);
            
            subplot(2,3,isexo*3+iarea);
            %         subplot(2,2,isexo*2+1);            
            mCue=mean(CueOnLFP{iarea},1);
            mm=mCue-mean(mCue);
            ee1=std(CueOnLFP{iarea},[],1)/sqrt(size(CueOnLFP{iarea},1));
            nn=mean(CueOffLFP{iarea},1);
            nn=nn-mean(nn);
            ee2=std(CueOffLFP{iarea},[],1)/sqrt(size(CueOffLFP{iarea},1));
            patchplot(TimeCue, mm, ee1,'r'); hold on;
            patchplot(TimeCue, nn, ee2,'b'); hold on;
            xlabel('Time to cue onset (ms)');
            ylabel('Evoked potential (LFP)');
            title(cuetype{isexo+1});
            axis tight;
            
%             subplot(2,2,isexo*2+2);           
%             mm=mean(ArrayOnLFP{iarea},1);  %-mean(mCue(TimeCue>-200&TimeCue<0));
%             ee1=std(ArrayOnLFP{iarea},[],1)/sqrt(size(ArrayOnLFP{iarea},1));
%             nn=mean(ArrayOffLFP{iarea},1);
%             ee2=std(ArrayOffLFP{iarea},[],1)/sqrt(size(ArrayOffLFP{iarea},1));
%             patchplot(TimeArray, mm, ee1, 'r'); hold on;
%             patchplot(TimeArray, nn, ee2, 'b'); hold on;
%             xlim([-300 600]);
%             xlabel('Time to array onset (ms)');
%             ylabel('Evoked potential (LFP)');
%             title(cuetype{isexo+1});
%             axis tight;
        end
    end
end
%%  get the power spectrum of LFP with hilbert transform
close all;
clear all;
monkeyset={'Mikey','Vasco'};
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118','110118','111618'};  % for Vasco
folder='Z:\RujiaChen\Results\';

TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
cuetype={'endo','exo'};
TimeCue=-600:1200;
Fs=1000;
for imonkey=1:2
    monkey=monkeyset{imonkey};
    if imonkey==2
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};
        channel=1:96;
    elseif imonkey==1
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
    end
    for igroup=1:2
        groupName=Group{igroup};
        for isexo=0:1
            count=zeros(1,3);
            for idate=1:numel(dateUsed)
                date=dateUsed{idate};
                load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
                %             load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);  %
                load([folder 'CorrectTrialParam_' date '.mat']);
                load([folder 'bArrayResp_' date '.mat']);
                load([folder 'bVisResp_' date '.mat']);
                load([folder 'bMotorResp_' date '.mat']);
                load([folder 'bVisMotor_' date '.mat']);
                load([folder 'CueOn_' date]);
                
                for iarea=1:3
                    if imonkey==1
                        isession=find(strcmp(RawInfo(1,:), date)==1);
                        channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
                        channelID=find(ismember(channel, RawInfo{iarea+1,isession}));  % for Mikey
                    else
                        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
                    end
                    disp(['date=' num2str(idate)]);
                    disp(['area=' num2str(iarea)]);
                    if ~isempty(channelID)
                        for ich=channelID
                            for ipos=1:2
                                if igroup==1
                                    celltype=bVisResp(ich, ipos);
                                elseif igroup==2
                                    celltype=bVisMotor(ich, ipos);
                                elseif igroup==3
                                    celltype=bMotorResp(ich, ipos);
                                end
                                if celltype>0
                                    if CueOn(ich,ipos)>0
                                        icue=CueOn(ich,ipos);
                                    else
                                        icue=imonkey;  % 1 for Mikey; 2 for Vasco
                                    end
                                    freq=3:2:55;
%                                     tCue=find(TimeCue>=-300&TimeCue<=600);
                                    
                                    idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                    data1=transpose(squeeze(LFPinfo.cueAlign(idx1,ich,:)));
%                                     refe1=transpose(squeeze(mean(LFPinfo.cueAlign(idx1,channelID,:),2)));
%                                     data1=data1-refe1;
                                          
                                    data2=transpose(squeeze(LFPinfo.ArrayAlign(idx1,ich,:)));
%                                     refe2=transpose(squeeze(mean(LFPinfo.ArrayAlign(idx1,channelID,:),2)));
%                                     data2=data2-refe2;
                                    

                                    idxBs=TimeCue>=-200&TimeCue<=0;
                                    baseResp=mean(data1(idxBs,:),1);
                                    bb=mean(baseResp);
                                    ee=std(baseResp);
                                    ValidTrl=abs(baseResp-bb)<5*ee;
                                   
                                    if sum(ValidTrl)>20
                                        num=0;
                                        count(iarea)=count(iarea)+1;                                        
                                        for ifreq=3:2:55
                                            num=num+1;                                                                     
                                            signal_filtered=bandpass(data2(:,ValidTrl), [ifreq-2 ifreq+2], Fs);
                                            signal_hilbert=hilbert(signal_filtered); 
                                            power_hb=mean(abs(signal_hilbert),2);
%                                             baseline_hb=mean(power_hb(1:300));
                                            Hilbert_power.CueOn{iarea}(count(iarea),num,:)=power_hb;  %(power_hb-baseline_hb)/baseline_hb;                                           
                                        end
                                        
                                        idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                        data1=transpose(squeeze(LFPinfo.cueAlign(idx2,ich,:)));
%                                         refe1=transpose(squeeze(mean(LFPinfo.cueAlign(idx2,channelID,:),2)));
%                                         data1=data1-refe1;
                                        
                                        data2=transpose(squeeze(LFPinfo.ArrayAlign(idx2,ich,:)));
%                                         refe2=transpose(squeeze(mean(LFPinfo.ArrayAlign(idx2,channelID,:),2)));
%                                         data2=data2-refe2;
                                        
                                        baseResp=mean(data1(idxBs,:),1);
                                        bb=mean(baseResp);
                                        ee=std(baseResp);
                                        ValidTrl=abs(baseResp-bb)<5*ee;
                                                                            
                                        num2=0;                                        
                                        for ifreq=3:2:55
                                            num2=num2+1;                                                                                          
                                            signal_filtered=bandpass(data2(:,ValidTrl), [ifreq-2 ifreq+2], Fs);
                                            signal_hilbert=hilbert(signal_filtered);  
                                            power_hb=mean(abs(signal_hilbert),2);
%                                             baseline_hb=mean(power_hb(1:300));
                                            Hilbert_power.CueOff{iarea}(count(iarea),num2,:)=power_hb;  %(power_hb-baseline_hb)/baseline_hb;                                                   
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
                        
            %         S.t=t-0.3;
            %         S.f=f_fft;
            Hilbert_power.f=3:2:55;
%             Hilbert_power.t=-600:1200;
            Hilbert_power.t=-900:900;
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_array_Hilbert_power_rawLFP.mat'],'Hilbert_power','-v7.3');
            clear Hilbert_power;
            %         save([folder monkey '_' cuetype{isexo+1} '_' groupName '_fft_power.mat'],'S','-v7.3');
        end
    end
end

%% plot the time course of mean power averaged across a frequency band
close all
AreaName={'LIP','PUL','FEF'};
freqBand=[5 5;8 20;21 35; 36 55];
monkeyset={'Mikey','Vasco'};
groupName={'VisResp','VisMotor','MotorResp'};
for imonkey=1:2
    monkey=monkeyset{imonkey};
    for isexo=0:1
        figure;  
        GroupData=cell(1,3);
        for igroup=1:2
            if exist([folder monkey '_' cuetype{isexo+1} '_' groupName{igroup} '_array_Hilbert_power_rawLFP.mat'],'file')   %
                load([folder monkey '_' cuetype{isexo+1} '_' groupName{igroup} '_array_Hilbert_power_rawLFP.mat']);   %  _referenceByAllChannels
                GroupData{igroup}=Hilbert_power; 
            end
        end
        data.t=Hilbert_power.t;
        data.f=Hilbert_power.f;
        for iarea=1:3
            data.CueOn{iarea}=cat(1,GroupData{1}.CueOn{iarea},GroupData{2}.CueOn{iarea});  %,GroupData{3}.CueOn{iarea}
            data.CueOff{iarea}=cat(1,GroupData{1}.CueOff{iarea},GroupData{2}.CueOff{iarea}); %,GroupData{3}.CueOff{iarea}
        end
        for ifreq=1:4
            for iarea=1:3
                subplot(4,3,(ifreq-1)*3+iarea);                
                t=data.t;
                f=data.f;               
                mmf=squeeze(mean(data.CueOn{iarea}(:,f>=freqBand(ifreq,1)&f<=freqBand(ifreq,2),:),2));
%                 mmf=squeeze(mean(data.CueOn{iarea}(:,:,t>300&t<500),3));
%                 mmf=log10(mmf);
                mm=mean(mmf,1);
                ee1=std(mmf,[],1)/sqrt(size(mmf,1));
                
                nnf=squeeze(mean(data.CueOff{iarea}(:,f>=freqBand(ifreq,1)&f<=freqBand(ifreq,2),:),2));
%                 nnf=squeeze(mean(data.CueOff{iarea}(:,:,t>300&t<500),3));
%                 nnf=log10(nnf);
                nn=mean(nnf,1);
                ee2=std(nnf,[],1)/sqrt(size(nnf,1));
                xx=t;
                patchplot(xx,mm,ee1,'r'); hold on;
                patchplot(xx,nn,ee2,'b'); hold on;
                title(AreaName{iarea});
                xlim([-500 700])
            end
        end
    end
end
%%  get the intra-areal spike_LFP coherence with fft
close all;
clear all;
monkeyset={'Mikey','Vasco'};
folder='Z:\RujiaChen\Results\';
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp'};    

params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];

TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape

TimeCue=-600:1200;
TimeArray=-900:900;
SpikeArray=-905:905;
for imonkey=1:2
    monkey=monkeyset{imonkey};
    if imonkey==2
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};
        channel=1:96;
    elseif imonkey==1
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
    end
    for igroup=3
        groupName=Group{igroup};
        for isexo=0:1
            count=zeros(1,3);
            clear S;
            clear C;
            for idate=1:numel(dateUsed)
                date=dateUsed{idate};
                load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
                load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);  %
                load([folder 'CorrectTrialParam_' date '.mat']);
                load([folder 'bArrayResp_' date '.mat']);
                load([folder 'bVisResp_' date '.mat']);
                load([folder 'bMotorResp_' date '.mat']);
                load([folder 'bVisMotor_' date '.mat']);
                load([folder 'CueOn_' date]);
                disp(['date=' num2str(idate)]);
                
                for iarea=1:3
                    disp(['area=' num2str(iarea)]);
                    if imonkey==1
                        isession=find(strcmp(RawInfo(1,:), date)==1);
                        channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
                        channelID=find(ismember(channel, RawInfo{iarea+1,isession}));  % for Mikey
                    else
                        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
                    end
       
                    if ~isempty(channelID)
                        for ich=channelID
                            for ipos=1:2
                                if igroup==1
                                    celltype=bVisResp(ich, ipos);
                                elseif igroup==2
                                    celltype=bVisMotor(ich, ipos);
                                elseif igroup==3
                                    celltype=bMotorResp(ich, ipos);
                                end
                                if celltype>0
                                    if CueOn(ich,ipos)>0
                                        icue=CueOn(ich,ipos);
                                    else
                                        icue=imonkey;  % 1 for Mikey; 2 for Vasco
                                    end                                    
%                                     tWin=find(TimeCue>=-300&TimeCue<=700);
                                    tWin=find(TimeArray>=-500&TimeArray<=700);
                                    tWin2=find(SpikeArray>=-500&SpikeArray<=700);
                                    
                                    DataSource1=LFPinfo.ArrayAlign;
                                    idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                    data1=transpose(squeeze(DataSource1(idx1,ich,tWin)));
                                    refe1=transpose(squeeze(mean(DataSource1(idx1,channelID,tWin),2)));
                                    data1=data1-refe1;
                                    baseResp=mean(data1(500:800,:),1);  % 500:800 on array condition, for array onset response; 1:300 on cue condition, for baseline
                                    bb=mean(baseResp);
                                    ee=std(baseResp);
                                    ValidTrl=abs(baseResp-bb)<3*ee;   % remove trials with outlier LFP
                                               
                                    DataSource2=SpikeTrain.arrayAlign;
                                    data0=transpose(squeeze(DataSource2(idx1,ich,tWin2)));
                                    ValidTrl=ValidTrl&mean(data0,1)~=0&~isnan(data0(1,:));   %
                                    if sum(ValidTrl)>20
                                        count(iarea)=count(iarea)+1;
                                        data=rmlinesc(data1(:,ValidTrl),params,[], 0, 60);
                                        [S.CueOn{iarea}(count(iarea),:,:),t, f_fft,~]=mtspecgramc(data, movingwin, params);
                                        [C.CueOn{iarea}(count(iarea),:,:),Phase.CueOn{iarea}(count(iarea),:,:),~,~,~,~,~,~]=cohgramcpb(data,data0(:, ValidTrl),movingwin,params);
                                        
                                        
                                        idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                        data1=transpose(squeeze(DataSource1(idx2,ich,tWin)));
                                        refe1=transpose(squeeze(mean(DataSource1(idx2,channelID,tWin),2)));
                                        data1=data1-refe1;
                                        baseResp=mean(data1(500:800,:),1);
                                        bb=mean(baseResp);
                                        ee=std(baseResp);
                                        ValidTrl=abs(baseResp-bb)<3*ee;   % remove trials with outlier LFP
                                        
                                        data0=transpose(squeeze(DataSource2(idx2,ich,tWin2)));                                        
                                        ValidTrl=ValidTrl&mean(data0,1)~=0&~isnan(data0(1,:));
                                        data=rmlinesc(data1(:,ValidTrl),params,[], 0, 60);
                                        [S.CueOff{iarea}(count(iarea),:,:),~, ~,~]=mtspecgramc(data, movingwin, params);
                                        [C.CueOff{iarea}(count(iarea),:,:),Phase.CueOff{iarea}(count(iarea),:,:),~,~,~,tc,fc,~]=cohgramcpb(data,data0(:, ValidTrl),movingwin,params);
                                        S.channel{iarea}(count(iarea),:)=[idate, ich];
                                        C.channel{iarea}(count(iarea),:)=[idate, ich];
                                    end  
                                  
                                end
                            end
                        end
                    end
                end
            end
            
            
            S.t=t-0.3;
            S.f=f_fft;
            C.phase=Phase;
            C.t=tc-0.3;
            C.f=fc;
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_coherence_referenceByAllChannels.mat'],'C','-v7.3');   % for S         
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_fft_power_referenceByAllChannels.mat'],'S','-v7.3');            
        end
    end
end

%%  plot the temporal-frequency power spectrum 
close all;
Group={'VisResp','VisMotor','MotorResp'}; 
AreaName={'LIP','PUL','FEF'};
igroup=2;
groupName=Group{igroup};
% zrange{1}=[0.05, 0.11;0.05, 0.12;0.05 0.12];
% zrange{2}=[0.08, 0.16;0.06, 0.14;0.05 0.12];
monkeyset={'Mikey','Vasco'};
folder='Z:\RujiaChen\Results\';
cuetype={'endo','exo'};

for isexo=0:1
    figure;
    for imonkey=1:2
        monkey=monkeyset{imonkey};        
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_coherence_referenceByAllChannels.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_fft_power_referenceByAllChannels.mat'])
        
        for iarea=1:3
            %             subplot(2,3,(imonkey-1)*3+iarea);
            subplot(3,2,(iarea-1)*2+imonkey);
            data=C;
            t=data.t-0.2;
            idx=~isnan(data.CueOn{iarea}(:,1,1));
            mm=transpose(squeeze(mean(data.CueOn{iarea}(idx,:,:),1)));
            nn=transpose(squeeze(mean(data.CueOff{iarea}(idx,:,:),1)));
%             base1=repmat(mean(mm(:,t<-0.1),2),1,size(mm,2));
%             base2=repmat(mean(nn(:,t<-0.1),2),1,size(nn,2));
%             mm=(mm-base1)./base1*100;
%             nn=(nn-base2)./base2*100;
            imagesc(t,data.f(data.f>10),mm(data.f>10,:));
            title(AreaName{iarea});
%             caxis([-15 15])
            %             caxis(zrange{imonkey}(iarea,:));
            axis xy;
            colorbar;
            
        end
    end
end

%%  get the inter-areal LFP_LFP coherence with fft
close all;
clear all;
monkeyset={'Mikey','Vasco'};
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118','110118','111618'};  % for Vasco
folder='Z:\RujiaChen\Results\';

params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];

TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
cuetype={'endo','exo'};
TimeCue=-600:1200;
TimeArray=-900:900;
tWin=find(TimeCue>=-300&TimeCue<=700); 
tArray=find(TimeArray>=-500&TimeArray<=600);
for imonkey=1:2
    monkey=monkeyset{imonkey};
    if imonkey==2
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};
        channel=1:96;
    elseif imonkey==1
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'}; % for Mikey
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
    end
    for igroup=1:3
        groupName=Group{igroup};
        for isexo=0:1
            count=zeros(2,3);
            clear Coh;
            for idate=1:numel(dateUsed)
                date=dateUsed{idate};
                load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
%                 load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);  %
                load([folder 'CorrectTrialParam_' date '.mat']);
                load([folder 'bArrayResp_' date '.mat']);
                load([folder 'bVisResp_' date '.mat']);
                load([folder 'bMotorResp_' date '.mat']);
                load([folder 'bVisMotor_' date '.mat']);
                load([folder 'CueOn_' date]);
                disp(['date=' num2str(idate)]);
                
                Data=LFPinfo.ArrayAlign;
                tAnalyse=tArray;
                for iarea=1:2
                    for jarea=(iarea+1):3
                        disp(['iarea=' num2str(iarea) ', jarea=' num2str(jarea)]);                        
                        if imonkey==1
                            isession=find(strcmp(RawInfo(1,:), date)==1);
                            channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
                            channelID1=find(ismember(channel, RawInfo{iarea+1,isession}));  % for Mikey
                            channelID2=find(ismember(channel, RawInfo{jarea+1,isession}));  % for Mikey
                        else
                            channelID1=(iarea-1)*32+1:iarea*32;  % for Vasco
                            channelID2=(jarea-1)*32+1:jarea*32;  % for Vasco
                        end
                        
%                         if ~isempty(channelID1) && ~isempty(channelID2)
                        if numel(channelID1)>=3&&numel(channelID2)>=3
                            for ich=channelID1(2:end-1)
                                for jch=channelID2(2:end-1)
                                    for ipos=1:2
                                        if igroup==1
                                            celltype=bVisResp(ich, ipos)&bVisResp(jch, ipos);
                                        elseif igroup==2
                                            celltype=bVisMotor(ich, ipos)&bVisMotor(jch, ipos);
                                        elseif igroup==3
                                            celltype=bMotorResp(ich, ipos)&bMotorResp(jch, ipos);
                                        end
                                        if celltype>0
                                            %                                         if CueOn(ich,ipos)>0
                                            %                                             icue=CueOn(ich,ipos);
                                            %                                         else
                                            icue=imonkey;  % 1 for Mikey; 2 for Vasco
                                            %                                         end                                                                                       
                                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;                                            
                                            data1=transpose(squeeze(Data(idx1,ich,tAnalyse)));
%                                             refe1=transpose(squeeze(mean(LFPinfo.cueAlign(idx1,channelID1,tCue),2)));
%                                             data1=data1-refe1; 
                                            data_i1=transpose(squeeze(Data(idx1,ich-1,tAnalyse)));
                                            data_i2=transpose(squeeze(Data(idx1,ich+1,tAnalyse)));
                                            data1=data_i1+data_i2-2*data1;  %% get the CSD based on the channel i-1 and i+1
                                                                                       
                                            baseResp=mean(data1(1:300,:),1);
                                            bb=mean(baseResp);
                                            ee=std(baseResp);
                                            ValidTrl1=abs(baseResp-bb)<3*ee;  % remove trials with outlier LFP
                                            
                                            data2=transpose(squeeze(Data(idx1,jch,tAnalyse)));
%                                             refe2=transpose(squeeze(mean(LFPinfo.cueAlign(idx1,channelID2,tCue),2)));
%                                             data2=data2-refe2;                                            
                                            data_j1=transpose(squeeze(Data(idx1,jch-1,tAnalyse)));
                                            data_j2=transpose(squeeze(Data(idx1,jch+1,tAnalyse)));
                                            data2=data_j1+data_j2-2*data2;  %% get the CSD based on the channel i-1 and i+1
                                            
                                            baseResp=mean(data2(1:300,:),1);
                                            bb=mean(baseResp);
                                            ee=std(baseResp);
                                            ValidTrl2=abs(baseResp-bb)<3*ee;  % remove trials with outlier LFP
                                            %                                         data0=transpose(squeeze(SpikeTrain.cueAlign(idx1,ich,tCue)));
                                            %                                         ValidTrl=ValidTrl&mean(data0,1)~=0&~isnan(data0(1,:));   %
                                            if sum(ValidTrl1&ValidTrl2)>20 
                                                ValidTrl=ValidTrl1&ValidTrl2;
                                                count(iarea,jarea)=count(iarea,jarea)+1;
                                                data_i=rmlinesc(data1(:,ValidTrl),params,[], 0, 60);
                                                data_j=rmlinesc(data2(:,ValidTrl),params,[], 0, 60);                                             
                                                [Coh.CueOn{iarea,jarea}(count(iarea,jarea),:,:),Phase.CueOn{iarea,jarea}(count(iarea,jarea),:,:),~,~,~,t,f]=cohgramc(data_i, data_j, movingwin, params);
                                                
                                                
                                                idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                                                data1=transpose(squeeze(Data(idx2,ich,tAnalyse)));
%                                                 refe1=transpose(squeeze(mean(LFPinfo.cueAlign(idx2,channelID1,tCue),2)));
%                                                 data1=data1-refe1;
                                                data_i1=transpose(squeeze(Data(idx2,ich-1,tAnalyse)));
                                                data_i2=transpose(squeeze(Data(idx2,ich+1,tAnalyse)));
                                                data1=data_i1+data_i2-2*data1;  %% get the CSD based on the channel i-1 and i+1
                                            
                                                baseResp=mean(data1(1:300,:),1);
                                                bb=mean(baseResp);
                                                ee=std(baseResp);
                                                ValidTrl1=abs(baseResp-bb)<3*ee;  % remove trials with outlier LFP
                                                
                                                data2=transpose(squeeze(Data(idx2,jch,tAnalyse)));
%                                                 refe2=transpose(squeeze(mean(LFPinfo.cueAlign(idx2,channelID2,tCue),2)));
%                                                 data2=data2-refe2;
                                                data_j1=transpose(squeeze(Data(idx2,jch-1,tAnalyse)));
                                                data_j2=transpose(squeeze(Data(idx2,jch+1,tAnalyse)));
                                                data2=data_j1+data_j2-2*data2;  %% get the CSD based on the channel i-1 and i+1
                                                
                                                baseResp=mean(data2(1:300,:),1);
                                                bb=mean(baseResp);
                                                ee=std(baseResp);
                                                ValidTrl2=abs(baseResp-bb)<3*ee;  % remove trials with outlier LFP
                                                ValidTrl=ValidTrl1&ValidTrl2;
                                                
                                                data_i=rmlinesc(data1(:,ValidTrl),params,[], 0, 60);
                                                data_j=rmlinesc(data2(:,ValidTrl),params,[], 0, 60);
                                                [Coh.CueOff{iarea,jarea}(count(iarea,jarea),:,:),Phase.CueOff{iarea,jarea}(count(iarea,jarea),:,:),~,~,~,~,~]=cohgramc(data_i, data_j, movingwin, params);
                                            end
                                        end
                                        
                                    end
                                end
                            end
                            disp(['count(' num2str(iarea) ', ' num2str(jarea) ')= ' num2str(count(iarea,jarea))]);
                        end
                    end
                end
            end
            
            Coh.phase=Phase;
            Coh.t=t-0.5;
            Coh.f=f;
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_interareal_coherence_CSD.mat'],'Coh','-v7.3');
    
        end
    end
end

%%  plot the inter-areal LFP-LFP coherence
close all;
AreaName={'LIP','PUL','FEF'};
Group={'VisResp','VisMotor','MotorResp'}; 
zrange=[-25 40; -25 50];
igroup=1;
groupName=Group{igroup};
for isexo=0:1
    figure;
    for imonkey=1:2
        monkey=monkeyset{imonkey};        
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ToArray_interareal_coherence_CSD.mat']);   %LFP_referenceByAllChannels
        num=0;
        for iarea=1:2
            for jarea=iarea+1:3
                num=num+1;
                subplot(3,2,(num-1)*2+imonkey);
                data=Coh;
                t=data.t;
                idx=~isnan(data.CueOn{iarea,jarea}(:,1,1));
                sCoh=squeeze(mean(mean(data.CueOn{iarea,jarea}(:,t>-0.05&t<0.05,:),2),3));
                ValidIdx=sCoh<0.2;   %abs(sCoh-mean(sCoh))<2.5*std(sCoh);
%                 plot(sCoh(ValidIdx),'*b'); hold on;
%                 plot(sCoh(ValidIdx==0),'*r'); hold on;
                
                mm=transpose(squeeze(mean(data.CueOn{iarea,jarea}(idx&ValidIdx,:,:),1)));
%                 ee=transpose(squeeze(std(mean(data.CueOn{iarea,jarea}(idx&ValidIdx,t<-0.07,:),2),[],1)));
%                 base1=repmat(mean(mm(:,t<-0.07),2),1,size(mm,2));
%                 mm=(mm-base1)./base1*100;
                nn=transpose(squeeze(mean(data.CueOff{iarea,jarea}(idx&ValidIdx,:,:),1)));
%                 base2=repmat(mean(nn(:,t<-0.07),2),1,size(nn,2));
%                 nn=(nn-base2)./base2*100;
                imagesc(t,data.f,mm);
%                 caxis([-8 8]);
                axis xy;
%                 xlim([0.2 0.7])
                colorbar;
                title([AreaName{iarea} '-' AreaName{jarea}]);
            end
        end        
    end
end

%%  plot the LFP traces for each session and find the contacts should be within the regions
close all;
AreaName={'LIP','PUL','FEF'};
iarea=3;
EstimateRange=zeros(10, 32);
siteN=[10 10 15];
for ii=1:size( mLFParray,1)
    figure;
    for jj=1:2
        idx=(iarea-1)*32+1:iarea*32;
        yy=31:-1:0;
        xx=-600:600;
        zz=mLFParray{ii,jj}(idx,:);
        [peakVal(ii,jj), peakCh(ii,jj)]=min(min(zz(:, xx>50 & xx<200), [], 2));
%         [peakMean(ii,jj), peakMeanCh(ii,jj)]=min(mean(zz(:, xx>50 & xx<200), 2));
%         baseline=mean(zz(:, xx>-250&xx<-50),2);
%         mm=zz-repmat(baseline, 1, size(zz,2));

        subplot(1,2,jj);
        imagesc(xx,yy,zz);  hold on;
        axis xy;
        line([-200 500], [yy(peakCh(ii,jj)) yy(peakCh(ii,jj))], 'color', 'k', 'linestyle', '--');
%         line([-200 500], [yy(peakMeanCh(ii,jj)) yy(peakMeanCh(ii,jj))], 'color', 'r', 'linestyle', '--');
        xlim([-100 300]);        
%         caxis([-0.5 0])
        colorbar;
    end
    if peakCh(ii,1)==peakCh(ii,2)
        tt=peakCh(ii,1);
    else
        tt=floor(min(peakCh(ii,:)));
    end
    rr=[floor(max([tt-siteN(iarea), 1]))  ceil(min([tt+siteN(iarea), 32]))];
    EstimateRange(ii, rr(1):rr(2))=1;

    title(dateUsed{ii});
end
EstimateRange([6, 10],:)=0;
save([folder 'EstimateRange_' AreaName{iarea} '.mat'], 'EstimateRange', '-v7.3');


%% remove the outliers && plot the coherence with sliding window 
close all;
folder='Z:\RujiaChen\Results\';
cueTy='exo';
load([folder 'LFPSpectrum_Mikey_SingleUnit_movWin_' cueTy '_VisualResp.mat']); 
load([folder 'Coherence_Mikey_SingleUnit_movWin_' cueTy '_VisualResp.mat']);  
load([folder 'TrlNumber_Mikey_SingleUnit_' cueTy '_VisualResp.mat'])
f=C.f;
t=C.t-0.9; % for array
% t=C.t-0.6;  % for cue
figure;
data=S;
conditions=fields(data);
area={'LIP', 'PUL', 'FEF'};
Ncell=zeros(1,3);
InvalidUnit=cell(1,3);
clear specPower;
for iarea=1:3   
    InvalidUnit{iarea}=zeros(4, size(data.(conditions{1}){iarea},1));
    for kk=1:4
        subplot(3,4,iarea*4-4+kk);
        idx1=find(transpose(~isnan(squeeze(mean(mean(data.(conditions{kk}){iarea},2),3)))));   %&TrlNumber{iarea}(1,:)<30
        Ncell(iarea)=numel(idx1);
        rr=squeeze(data.(conditions{kk}){iarea}(idx1,9,f>40));
%         plot(f(f>40), rr'); hold on;
        for icell=1:size(rr,1)
            specPower{iarea,kk}(:,icell)=abs(fft(rr(icell,:)));
        end
        mm=mean(specPower{iarea,kk},1);
        mVal=mean(mm);
        ee=std(mm);   
        [~, maxID]=max(specPower{iarea,kk}(7:12,:),[],1);
        idx00=mm-mVal>1.5*ee| maxID==4;
        InvalidUnit{iarea}(kk, idx1(idx00))=1;
        plot(specPower{iarea,kk}(:, idx00), 'k'); hold on;
        plot(specPower{iarea,kk}(:, idx00==0)); hold on;
        axis tight;
        title([area{iarea} '-' conditions{kk}(4:end)]);
    end
end
 
%%%%%%%% plot the population averaged coherence 
data=S;
figure;
for iarea=1:3
    subplot(3,3,(iarea-1)*3+1);
    kk=2;
    idx1=find(transpose(~isnan(squeeze(mean(mean(data.(conditions{kk}){iarea},2),3)))&~isnan(squeeze(mean(mean(data.(conditions{kk+2}){iarea},2),3))))&sum(InvalidUnit{iarea},1)==0);   %&C.attModulation{iarea}==1    %&TrlNumber{iarea}(1,:)<30   &data.SU{iarea}==1
    Ncell(iarea)=numel(idx1); 
    zz=squeeze(mean(data.(conditions{kk}){iarea}(idx1,:,:),1));
    xx=transpose(zz(:, f<100));
% %     xx=log(xx);  % for Spectrum
%     imagesc(t, f(f<100), xx); hold on;   
    
%     smoothKernal=ones(3)*0.5;
%     smoothKernal(2,2)=1;
%     smoothKernal=smoothKernal./5;
    smoothKernal=1;
    SmoothedData1=zeros(numel(idx1), size(xx,1), size(xx,2));
    SmoothedData2=zeros(numel(idx1), size(xx,1), size(xx,2));
    for icell=1:numel(idx1)
        SmoothedData1(icell,:,:)=transpose(conv2(squeeze(data.(conditions{kk}){iarea}(idx1(icell),:,f<100)), smoothKernal,'same'));
        SmoothedData2(icell,:,:)=transpose(conv2(squeeze(data.(conditions{kk+2}){iarea}(idx1(icell),:,f<100)), smoothKernal,'same'));
    end
    meanSmth=squeeze(mean(SmoothedData1,1));
    imagesc(t, f(f<100), meanSmth); hold on;    
    axis xy;
    cc=colorbar;
    title([area{iarea} '-' conditions{kk}(4:end)]);
    xlim([-0.7 0.7]);
    
    subplot(3,3,(iarea-1)*3+2);
    zz=squeeze(mean(data.(conditions{kk+2}){iarea}(idx1,:,:),1));
    yy=transpose(zz(:, f<100));
% %     yy=log(yy);  % for Spectrum
%     imagesc(t, f(f<100), yy); hold on;
    
    meanSmth=squeeze(mean(SmoothedData2,1));
    imagesc(t, f(f<100), meanSmth); hold on;
    
    axis xy;
    caxis(cc.Limits);
    colorbar;   
    title([area{iarea} '-' conditions{kk+2}(4:end)]);
   xlim([-0.7 0.7]);
    
    subplot(3,3,(iarea-1)*3+3);
    sigDiff=zeros(size(yy));
    for itime=1:size(yy,2)
        for ifreq=1:size(yy,1)
            data1=squeeze(SmoothedData1(:,ifreq, itime));
            data2=squeeze(SmoothedData2(:,ifreq, itime));
            sigDiff(ifreq, itime)=ttest(data1, data2);
        end
    end
    
%     diffCoh=log(xx)-log(yy);   % for spectrum
    diffCoh=xx-yy;
    diffCoh=conv2(diffCoh,smoothKernal,'same');
    imagesc(t, f(f<100), diffCoh); hold on;   %.*sigDiff
    axis xy;
%     caxis([-0.05 0.07]);
    colorbar;
xlim([-0.7 0.7]);
    title([area{iarea} '-(On-Away)']);
end

%%
%%%%%% plot the difference of coherence between cue-on and cue-away
%%%%%% conditions for each unit
close all;
for iarea=1:3
    figure;
    kk=2;
    idx1=find(transpose(~isnan(squeeze(mean(mean(data.(conditions{kk}){iarea},2),3)))&~isnan(squeeze(mean(mean(data.(conditions{kk+2}){iarea},2),3))))&data.SU{iarea}==1&sum(InvalidUnit{iarea},1)==0);   %    %&TrlNumber{iarea}(1,:)<30    &data.attModulation{iarea}==1
    Ncell(iarea)=numel(idx1);
    diffData=data.(conditions{kk}){iarea}(idx1,:,:)-data.(conditions{kk+2}){iarea}(idx1,:,:);
    for icell=1:min([Ncell(iarea), 25])  %size(diffData)
        subplot(5,5,icell);
        yy=transpose(squeeze(diffData(icell,:,f<100)));
        imagesc(t, f(f<100), yy); hold on;
        caxis([-0.2 0.2]);
    end
%     colorbar;
end




%% plot the coherence in delay or beforeonset period in exo or endo conditions
% clear all;
close all;
folder='Z:\RujiaChen\Results\';
% load([folder 'mLFP.mat']);
% load([folder 'LFPSpectrum_endo.mat']);
CueType='endo';
duration='delay';
load([folder 'LFPSpectrum_Mikey_SingleUnit_' duration '_' CueType '_ArrayResp.mat']);  % delay_endo_facilitateUnit   CueOn  delay
load([folder 'Coherence_Mikey_SingleUnit_' duration '_' CueType '_ArrayResp.mat']);  %beforeOnset     endo_AllTargetResponsive
% load([folder 'Coherence_endo_baseline.mat']);
% load([folder 'Coherence_endo_baseline_AllTargetResponsive.mat']);
load([folder 'f.mat']);
% load([folder 'Coherence_endo.mat']);
% load([folder 'Coherence_beforeOnset_exo.mat']);
% load([folder 'LFPSpecGram.mat']);
% load([folder 'Coherence_f.mat']);
% load([folder 'SpikeSpectrum_exo.mat']);  %_endo
% load([folder 'f_spike.mat']);
load([folder 'TrlNumber_Mikey_SingleUnit_delay_' CueType '_ArrayResp.mat']);
%%%%% plot the coherence in different areas and durations
figure;
clr=[1 0 0; 0 0 1];
data=C;  % plot the spectrum or coherence 
% Baseline=zeros(3,2,length(f));
for iarea=1:3
    conditions=fields(data);    
    sigDiff=zeros(length(f),2);
    for ifreq=1:length(f)
        idx1=TrlNumber{iarea}(1,:)>=20;
        idx2=TrlNumber{iarea}(2,:)>=20;
        sigDiff(ifreq,1)=ttest2(data.(conditions{1}){iarea}(idx1,ifreq), data.(conditions{3}){iarea}(idx2, ifreq)); 
        sigDiff(ifreq,2)=ttest2(data.(conditions{2}){iarea}(idx1,ifreq), data.(conditions{4}){iarea}(idx2, ifreq)); 
    end
    for ii=1:4  %numel(conditions)
%         subplot(2,2,ii);       
%         mm=mLFP.(conditions{ii}){iarea}-repmat(mean(mLFP.(conditions{ii}){iarea}(:,1:500),2),1,size(mLFP.(conditions{ii}){iarea},2));
%         mm=mm./repmat(max(abs(mm),[],2),1,size(mm,2));
% %         imagesc(mm);
%         meanLFP=mean(mm);
%         plot(meanLFP);
%         axis xy;
% %         colorbar;    
        ifig=(iarea-1)*2+2-mod(ii,2);        
        subplot(3,2,ifig);
        idx=find(TrlNumber{iarea}(1, :)>=20&TrlNumber{iarea}(2, :)>=20&data.SU{iarea}==1);  %&data.attModulation{iarea}==1);  %  
        yy=data.(conditions{ii}){iarea};
%         yy=yy-repmat(transpose(squeeze(Baseline(iarea,ceil(ii/2),:))), size(yy,1), 1);  % only for delay condition
        nn=mean(yy(idx,:),1);
        if mod(ii,2)==1
            Baseline(iarea, ceil(ii/2),:)=nn;
        end
        Ncell(iarea)=numel(idx);
        ee=std(yy(idx,:),0,1)/sqrt(numel(idx));
        idx1=f>0&f<100;
        xx=f(idx1);
        Vx=[xx fliplr(xx)  xx(1)];
        yy1=nn(idx1)-ee(idx1);
        yy2=nn(idx1)+ee(idx1);
        Vy=[yy1 fliplr(yy2) yy1(1)];
        patch(Vx,Vy, clr(ceil(ii/2),:),'edgecolor',clr(ceil(ii/2),:),'facealpha',0.5,'edgealpha',0.5); hold on;
        plot(xx, nn(idx1), 'color', clr(ceil(ii/2),:), 'linewidth', 1.2); hold on;   
        idxx=idx1'&sigDiff(:,2-mod(ii,2))==1;
        plot(f(idxx), sigDiff(idxx,2-mod(ii,2))*0, 'k*'); hold on;      
        axis tight;
        if ifig>4
            xlabel('Frequency (Hz)');
            ylabel('Power spectrum');
%             ylabel('Spike-field coherence');
        end
        if ifig==1
            title('Cue');
        elseif ifig==2
            title('Array');
        end
        set(gca,'box','off');
    end
end

% save([folder 'PowerSpec_baseline_Mikey_SingleUnit_' CueType '_ArrayResp.mat'],'Baseline', '-v7.3');
% save([folder 'Coherence_baseline_Mikey_SingleUnit_' CueType '_ArrayResp.mat'],'Baseline', '-v7.3');
%%  plot the coherence before array onset with the baseline (before cue onset) subtracted 

% load([folder 'Coherence_beforeOnset_exo_AllTargetResponsive.mat']);  %beforeOnset   
% load([folder 'Coherence_beforeOnset_exo_CueOn.mat']);  %beforeOnset     
load([folder 'LFPSpectrum_Mikey_SingleUnit_beforeOnset_exo_ArrayResp.mat']);  % delay_endo_facilitateUnit   CueOn  delay
load([folder 'Coherence_Mikey_SingleUnit_beforeOnset_exo_ArrayResp.mat']);  %beforeOnset     endo_AllTargetResponsive
% close all;
figure;
data=C;
for iarea=2
    conditions=fields(data);    
    sigDiff=zeros(length(f),1);
%     idx1=TrlNumber{iarea}(1,:)>=30;
%     idx2=TrlNumber{iarea}(2,:)>=30;
    yy=data.(conditions{2}){iarea}-data.(conditions{1}){iarea};  % subtract the baseline before cue onset   (idx1,:)-
    zz=data.(conditions{4}){iarea}-data.(conditions{3}){iarea};
    for ifreq=1:length(f)
        sigDiff(ifreq)=ttest2(yy(:,ifreq), zz(:,ifreq));         
    end
     SingleSession=zeros(2,size(yy,1), size(yy,2));
%     subplot(3,1,iarea);
    for ii=2:2:4
%         idx=TrlNumber{iarea}(ceil(ii/2), :)>=30;
        NormData=data.(conditions{ii}){iarea};  %(idx,:);  %-data.(conditions{ii-1}){iarea}(idx,:);  % subtract the baseline before cue onset
        SingleSession(round(ii/2), :,:) = NormData;
        nn=mean(NormData, 1);
        ee=std(NormData,0,1)/sqrt(size(NormData,1));
%         ee=std(data.(conditions{ii}){iarea}(idx,:),0,1)/sqrt(size(data.(conditions{ii}){iarea}(idx,:),1));
        idx1=f>0&f<100;
        xx=f(idx1);
        Vx=[xx fliplr(xx)  xx(1)];
        yy1=nn(idx1)-ee(idx1);
        yy2=nn(idx1)+ee(idx1);
        Vy=[yy1 fliplr(yy2) yy1(1)];
        patch(Vx,Vy, clr(ceil(ii/2),:),'edgecolor',clr(ceil(ii/2),:),'facealpha',0.5,'edgealpha',0.5); hold on;
        plot(xx, nn(idx1), 'color', clr(ceil(ii/2),:), 'linewidth', 1.2); hold on;   
        axis tight;
        if iarea==3
            xlabel('Frequency (Hz)');
        elseif iarea==2
            ylabel('Spike-field coherence');
        end
        set(gca,'box','off');
    end
    idxx=idx1'&sigDiff==1;
    plot(f(idxx), sigDiff(idxx)*0, 'k*','linewidth',1); hold on;
end


%%%%% plot the coherence for each session
clr2=[1 0 0; 0 0 1];
for isession=1:size(SingleSession,2)
    subplot(5,7, isession);
    for ii=1:2
        yy=squeeze(SingleSession(ii, isession,:));
        idx= f>0&f<100;
        plot(f(idx), yy(idx), 'color', clr2(ii,:));  hold on;
    end
    
end
%% plot the mean temporal-frequency specgram of different conditions
for iarea=1:3
    figure;
    subplot(2,2,1)
    mm=squeeze(mean(S.exoOnCue{iarea},1));
    idx=f<100;
    imagesc(t*1000-600,f(idx),log10(mm(:,idx)'));
    axis xy;
%     caxis([0 0.3]);
    colorbar;
    subplot(2,2,2)
    nn=squeeze(mean(S.exoOffCue{iarea},1));
    idx=f<100;
    imagesc(t*1000-600,f(idx),mm(:,idx)'-nn(:,idx)');
    axis xy;
    % caxis([0 0.5]);
    colorbar;
    
    subplot(2,2,3)
    mm=squeeze(mean(S.exoOnArray{iarea},1));
    idx=f<100;
    imagesc(t*1000-600,f(idx),log10(mm(:,idx)'));
    axis xy;
%     caxis([0 0.3]);
    colorbar;
    subplot(2,2,4)
    nn=squeeze(mean(S.exoOffArray{iarea},1));
    idx=f<100;
    imagesc(t*1000-600,f(idx),mm(:,idx)'-nn(:,idx)');
    axis xy;
    % caxis([0 0.5]);
    colorbar;
end