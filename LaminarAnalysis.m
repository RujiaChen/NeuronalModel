
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';

imonkey=2;
monkey=monkeyset{imonkey};
if strcmp(monkey, 'Mikey')
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
elseif strcmp(monkey, 'Vasco')
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
end

Time=-600:1200;   % for cue spike&LFP
Time1=-505:1305;  % for array spike
iarea=3;
TimeZero='ToArray';
CueCondition='Exo';

close all;
for ifile = 1:numel(dateUsed) %1 :numel(filename)
    date=dateUsed{ifile};
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);
    load([folder 'CorrectTrialParam_' date '.mat']);
    load([folder 'SaccadeTime_'  TimeZero '_' CueCondition '_' date '.mat']);
    load([folder 'RFOn_' date '.mat']);
    load([folder 'SaccadeOn_' date '.mat']);
                    
    if imonkey==1
        load([folder 'Mikey_RecordingInfo.mat']);
        isession=find(strcmp(RawInfo(1,:), date)==1);
        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
    else
        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
    end
    
    TimePlt=-300:600;
    Time_sac = -600:200;
    idxT = Time>=TimePlt(1)&Time<=TimePlt(end);
    idxT_array = Time1>=TimePlt(1)&Time1<=TimePlt(end);
    CueResp = zeros(numel(channelID),length(TimePlt));
    ArrayResp = zeros(numel(channelID),length(TimePlt));
    SaccadeResp = zeros(numel(channelID),length(Time_sac));
    nch=0;
    for ich = channelID
        if RFOn(ich)>0
            ipos=ceil(RFOn(ich)/2);
            icue=2-mod(RFOn(ich),2);
        elseif SaccadeOn(ich)>0
            ipos=ceil(SaccadeOn(ich)/2);
            icue=2-mod(SaccadeOn(ich),2);
        else
            ipos=1;
            icue=imonkey;
        end
        nch = nch +1;
        
        idxTrl=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);   
        CueResp(nch,:) = squeeze(mean(SpikeTrain.cueAlign(idxTrl,ich,idxT),1));         
        ArrayResp(nch,:) = squeeze(mean(SpikeTrain.arrayAlign(idxTrl,ich,idxT_array),1));
        
        ntrl=0;
        clear motorResp;
        motorResp=[];                        
        for itrl=1:numel(idxTrl)
            if SaccadeTime{ipos, icue}(itrl)>400&&SaccadeTime{ipos, icue}(itrl)<=1100
                ntrl=ntrl+1;
                idxT_motor=Time1-SaccadeTime{ipos, icue}(itrl)>=Time_sac(1)&Time1-SaccadeTime{ipos, icue}(itrl)<=Time_sac(end);
                motorResp(ntrl,:)=squeeze(SpikeTrain.arrayAlign(idxTrl(itrl),ich,idxT_motor));
            end
        end        
        SaccadeResp(nch,:)=mean(motorResp,1);        
    end
    mBase=mean(CueResp(:,Time>-200&Time<0),2);
    CueResp=conv2(CueResp-repmat(mBase,1,size(CueResp,2)), ones(1,50)/50,'same');
    ArrayResp=conv2(ArrayResp-repmat(mBase,1,size(ArrayResp,2)), ones(1,50)/50,'same');
    SaccadeResp=conv2(SaccadeResp-repmat(mBase,1,size(SaccadeResp,2)), ones(1,50)/50,'same');
%     ArrayResp=conv2(ArrayResp, ones(1,50)/50,'same');   % for attention modulation
            
    figure;
    subplot(1,2,1);
    chPlt= SpikeTrain.channel(channelID);    
%     imagesc(TimePlt,chPlt,CueResp(:,idxT));    
    imagesc(TimePlt,chPlt,ArrayResp);
    colorbar;
    xlabel('Time to array onset (ms)');       
    
    subplot(1,2,2);
    imagesc(Time_sac,chPlt,SaccadeResp);
    colorbar;
    xlabel('Time to saccade onset (ms)');    
    
end

%% get the CSD of individule session in different areas
close all;
TimeCue=-600:1200;
TimeArray=-900:900;
for ifile = [3 5 9]  %numel(dateUsed) %1 :numel(filename)
    date=dateUsed{ifile};
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
    load(['Z:\Rujia\Results\CorrectTrialParam_' date '.mat']);
    load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);   
    if imonkey==1
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
        isession=find(strcmp(RawInfo(1,:), date)==1);
        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
    else
        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
    end
        
    icueP=imonkey;    
    idxTrl=find(CorrectTrialParam.cueloc==cueLoc(icueP) &CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0);     
    chPlt= SpikeTrain.channel(channelID);    

    CueLFP=squeeze(mean(LFPinfo.cueAlign(idxTrl,channelID,TimeCue>=-100&TimeCue<=300),1));   % for LFP
    ArrayLFP=squeeze(mean(LFPinfo.ArrayAlign(idxTrl,channelID,TimeArray>=-100&TimeArray<=300),1));  
    CSD = CSD_stand(CueLFP,1,1);
%     CSD = CSD_stand(ArrayLFP,1,1);
    
    
    TimePlt=-100:300;    
%     idxT=TimePlt>=80&TimePlt<=100;
%     mCSD=mean(CSD(:,idxT),2);
%     [~, sinkCenter]=min(mCSD);
%     borderID=-1;
%     for ielec=sinkCenter:sinkCenter+5
%         if borderID<0 && mCSD(ielec)<0&&mCSD(ielec+1)>0 %&&mCSD(ielec+2)>0 &&mCSD(ielec-1)<0
%             borderID=ielec;
%             break;
%         end
%     end
    figure;
    imagesc(TimePlt, chPlt, CSD); hold on;
    line([100 100], [chPlt(1) chPlt(end)]); hold on;
%     plot(mCSD, chPlt,'r'); hold on;
%     line([TimePlt(1) TimePlt(end)],[chPlt(borderID) chPlt(borderID)],'color','r','linewidth',1.5); hold on;
    caxis([-2 2]);
    colorbar;
    xlabel('Time to cue onset (ms)');        
end

%% get the laminar profile of temporal power spectrum with hilbert transform
close all;
TimeCue=-600:1200;
TimeArray=-900:900;
Fs=1000;
FreqBand=[3, 7;8, 15; 16, 25; 30, 60];
iarea=1;
for ifile = [1 3 5 6]  %numel(dateUsed) %1 :numel(filename)
    date=dateUsed{ifile};
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
    load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);   
    if imonkey==1
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
        isession=find(strcmp(RawInfo(1,:), date)==1);
        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
    else
        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
    end
        
    icueP=imonkey;    
    idxTrl=find(CorrectTrialParam.cueloc==cueLoc(icueP) &CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0);     
    chPlt= SpikeTrain.channel(channelID);  
    
    Hilbert_power=zeros(numel(chPlt), numel(TimeCue));
    for ich=1:numel(chPlt)
        data1=squeeze(LFPinfo.cueAlign(idxTrl,chPlt(ich),:));           
        refe1=squeeze(mean(LFPinfo.cueAlign(idxTrl,chPlt,:),2));
        data1=data1-refe1;
        
        data2=squeeze(LFPinfo.ArrayAlign(idxTrl,chPlt(ich),:)); 
        refe2=squeeze(mean(LFPinfo.ArrayAlign(idxTrl,chPlt,:),2));
        data2=data2-refe2;
        
        idxBs=TimeCue>=-200&TimeCue<=0;
        baseResp=mean(data1(:,idxBs),2);
        bb=mean(baseResp);
        ee=std(baseResp);
        ValidTrl=abs(baseResp-bb)<5*ee;
        
        signal_filtered=bandpass(transpose(data1(ValidTrl,:)), FreqBand(2,:), Fs);
        signal_hilbert=hilbert(signal_filtered);
        power_hb=mean(abs(signal_hilbert),2);
%         baseline_hb=mean(power_hb(TimeCue>=-200&TimeCue<=0));
        Hilbert_power(ich, :)=power_hb;  %(power_hb-baseline_hb)/baseline_hb;               
    end
    
    figure;
%     TimePlt=TimeCue(TimeCue>=-400&TimeCue<=600);
    TimePlt=TimeArray(TimeArray>=-500&TimeArray<=600);
    imagesc(TimePlt, chPlt, Hilbert_power(:,idxT)); hold on;
    line([0 0], [chPlt(1) chPlt(end)]); hold on;
    colorbar;
    xlabel('Time to cue onset (ms)');        
end

%% get the laminar profile of temporal power spectrum with FFT
close all;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
FreqBand=[3, 7;8, 15; 16, 25; 30, 60];
num=0;
TimeCue=-600:1200;
TimeArray=-900:900;
clear Spec;
for ifile = [3 5 9] %numel(dateUsed) %1 :numel(filename)   
    num=num+1;
    date=dateUsed{ifile};
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
    load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);   
    if imonkey==1
        load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
        isession=find(strcmp(RawInfo(1,:), date)==1);
        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
    else
        channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
    end
        
    icueP=imonkey;    
    idxTrl=find(CorrectTrialParam.cueloc==cueLoc(icueP) &CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0);     
    idxTrl2=find(CorrectTrialParam.cueloc==cueLoc(3-icueP) &CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0);
    chPlt= SpikeTrain.channel(channelID);  
    
    Hilbert_power=zeros(numel(chPlt), numel(TimeCue));
    for ich=1:numel(chPlt)
        data1=squeeze(LFPinfo.cueAlign(idxTrl,chPlt(ich),:));  
%         data1=data1-repmat(mean(data1(:,TimeCue>=-300&TimeCue<=0),2),1,numel(TimeCue));
%         refe1=squeeze(mean(LFPinfo.cueAlign(idxTrl,chPlt,:),2));
%         data1=data1-refe1;
        
        data2=squeeze(LFPinfo.ArrayAlign(idxTrl,chPlt(ich),:)); 
        data2=data2-repmat(mean(data1(:,TimeCue>=-300&TimeCue<=0),2),1,numel(TimeCue));
%         refe2=squeeze(mean(LFPinfo.ArrayAlign(idxTrl,chPlt,:),2));
%         data2=data2-refe2;
        
        idxBs=TimeCue>=-200&TimeCue<=0;
        baseResp=mean(data1(:,idxBs),2);
        bb=mean(baseResp);
        ee=std(baseResp);
        ValidTrl=abs(baseResp-bb)<5*ee;        
%         data=rmlinesc(transpose(data1(ValidTrl,:)),params,[], 0, 57:63);
        data=rmlinesc(transpose(data2(ValidTrl,:)),params,[], 0, 57:63);
%         [S.CueOn(ich,:,:),~, ~,~]=mtspecgramc(data, movingwin, params);
%         [Spec{num}.CueOn(ich,:),~, ~]=mtspectrumc(data1(TimeCue>=250&TimeCue<=600,:),  params);
        [Spec{num}.CueOn(ich,:),~, ~]=mtspectrumc(data(TimeArray>=250&TimeArray<=600,:),  params);
              
        data1=squeeze(LFPinfo.cueAlign(idxTrl2,chPlt(ich),:));
%         data1=data1-repmat(mean(data1(:,TimeCue>=-300&TimeCue<=0),2),1,numel(TimeCue));
%         refe1=squeeze(mean(LFPinfo.cueAlign(idxTrl2,chPlt,:),2));
%         data1=data1-refe1;
        
        data2=squeeze(LFPinfo.ArrayAlign(idxTrl2,chPlt(ich),:)); 
        data2=data2-repmat(mean(data1(:,TimeCue>=-300&TimeCue<=0),2),1,numel(TimeCue));
%         refe2=squeeze(mean(LFPinfo.ArrayAlign(idxTrl2,chPlt,:),2));
%         data2=data2-refe2;
        
        idxBs=TimeCue>=-200&TimeCue<=0;
        baseResp=mean(data1(:,idxBs),2);
        bb=mean(baseResp);
        ee=std(baseResp);
        ValidTrl=abs(baseResp-bb)<5*ee;        
        data=rmlinesc(transpose(data2(ValidTrl,:)),params,[], 0, 57:63);
%         [S.CueOff(ich,:,:),t, f,~]=mtspecgramc(data, movingwin, params);
%         [Spec{num}.CueOff(ich,:),f, ~]=mtspectrumc(data(TimeCue>=250&TimeCue<=600,:),  params);
        [Spec{num}.CueOff(ich,:),~, ~]=mtspectrumc(data(TimeArray>=250&TimeArray<=600,:), params);
    end
    Spec{num}.f=f;
%     S.t=t-0.6;
%     S.f=f; 

%     figure;
%     for ifreq=1:4
%     subplot(1,4,ifreq);
%         mm=squeeze(mean(S.CueOn(:,:,f>=FreqBand(ifreq,1)&f<=FreqBand(ifreq,2)),3));
% %         mm=mm-repmat(mean(mm(:,S.t>=-0.3&S.t<=0),2),1,length(S.t));
%         imagesc(S.t, chPlt, mm);  %(S.t>-0.3)  (chPlt>=70&chPlt<=92) 
%         colorbar;
%     end
   
end

%%
close all;
FreqBand=[3, 7;8, 15; 16, 25; 30, 70];
for num=1
    figure;
%     subplot(1,2,1);
%     imagesc(Spec{num}.f, chPlt, log(Spec{num}.CueOn));    
%     colorbar;
%     subplot(1,2,2);
%     imagesc(Spec{num}.f, chPlt, Spec{num}.CueOn-Spec{num}.CueOff);
%     %     caxis([0 0.2])
%     colorbar;
    for ifreq=1:4
        subplot(1,4,ifreq);
        idx=chPlt>=70&chPlt<=93;
        mm=mean(Spec{num}.CueOn(idx,f>=FreqBand(ifreq,1)&f<=FreqBand(ifreq,2)),2);
        nn=mean(Spec{num}.CueOff(idx,f>=FreqBand(ifreq,1)&f<=FreqBand(ifreq,2)),2);
        yy=(mm-nn)./mm;
        
%         mm=mm-min(mm);
%         mm=mm./max(mm);
%         yy=mm;
      
        plot(flipud(yy),chPlt(idx),'linewidth',2.0); hold on;        
        line([min(yy) max(yy)],[76 76],'linestyle','--'); hold on;
%         line([min(yy) max(yy)],[86 86],'linestyle','--'); hold on;
        yticks([66 71 76 81 86 91 96]);
        yticklabels({'95','90','85','80','75','70', '65'});
        set(gca,'box','off');
    end
end

%%
figure;
for ich=1:size(S.CueOn,1)
    subplot(6,6,ich);
    mm=transpose(squeeze(S.CueOn(ich,:,:)));
    mm=mm-repmat(mean(mm(:,S.t>=-0.3&S.t<=0),2),1,length(S.t));
%     imagesc(S.t, S.f, mm );
%     axis xy;
%     caxis([0 0.3])
%     colorbar;

    nn=mean(mm(:, S.t>=0&S.t<=0.25),2);
    plot(S.f, nn); hold on;
end
%%



