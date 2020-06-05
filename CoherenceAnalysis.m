%% get the intra-areal spike-field coherence
clear all;
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
Time=-600:1200; 
Time1=-505:1305;
Time1_LFP=-900:900;
idxT_cue=Time>=-600&Time<=1000;
idx_array=Time1>=-500&Time1<=900;
idx_array_LFP=Time1_LFP>=-500&Time1_LFP<=900;
Selection={'','_5sigma','_SingleUnit'};
iselect=2;                            
for igroup=1:3
    groupName=Group{igroup};   
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end        
        for isexo=0:1                    
            SpikeLFPInfo.ElecInfo=cell(numel(dateUsed),3);
            SpikeLFPInfo.CueOnSpike=cell(numel(dateUsed),3);
            SpikeLFPInfo.ArrayOnSpike=cell(numel(dateUsed),3);   
            SpikeLFPInfo.CueOffSpike=cell(numel(dateUsed),3);
            SpikeLFPInfo.ArrayOffSpike=cell(numel(dateUsed),3);
            SpikeLFPInfo.CueOnLFP=cell(numel(dateUsed),3);
            SpikeLFPInfo.ArrayOnLFP=cell(numel(dateUsed),3);
            SpikeLFPInfo.CueOffLFP=cell(numel(dateUsed),3);
            SpikeLFPInfo.ArrayOffLFP=cell(numel(dateUsed),3);
            
            for iarea=1:3                
                for idate = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{idate};
                    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
                    load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
                    if iselect==2
                        load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);
                    elseif iselect==3
                        load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                    end
                    load([folder 'bVisResp_' date Selection{iselect} '.mat']);
                    load([folder 'bMotorResp_' date Selection{iselect} '.mat']);
                    load([folder 'bVisMotor_' date Selection{iselect} '.mat']);
                    load([folder 'CueOn_' date Selection{iselect} '.mat']);
                    load([folder 'RFOn_' date Selection{iselect} '.mat']);
                    load([folder 'SaccadeOn_' date Selection{iselect} '.mat']);
                 
                    if imonkey==1
                        isession=find(strcmp(RawInfo(1,:), date)==1);   
                        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
                        LFPChannel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}]; 
                        AreaCh=find(ismember(LFPChannel,SpikeTrain.channel(channelID)));
                    else
                        channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
                        LFPChannel=1:96;
                        AreaCh=find(ismember(LFPChannel,SpikeTrain.channel(channelID)));
                    end
                    num=0;
                    for ich=channelID
                        if igroup==1
                            celltype=bVisResp(ich)*CueOn(ich);  %, ipos
                        elseif igroup==2
                            celltype=bVisMotor(ich)*CueOn(ich);
                        elseif igroup==3
                            celltype=bMotorResp(ich);
                        elseif igroup==4
                            celltype=(bVisResp(ich)+bVisMotor(ich))*CueOn(ich); 
                        end
                        if celltype>0
                            num=num+1;                                
                            %%%%%% prime position
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
                            chID=find(LFPChannel==SpikeTrain.channel(ich));  % identify LFP elecID
                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;

                            SpikeLFPInfo.CueOnSpike{idate,iarea}{num}=squeeze(SpikeTrain.cueAlign(idx1,ich,idxT_cue));   
                            SpikeLFPInfo.ArrayOnSpike{idate, iarea}{num}=squeeze(SpikeTrain.arrayAlign(idx1,ich,idx_array)); 
                            
                            ReferenceLFP1=squeeze(mean(LFPinfo.cueAlign(idx1,AreaCh,idxT_cue),2));
                            SpikeLFPInfo.CueOnLFP{idate,iarea}{num}=squeeze(LFPinfo.cueAlign(idx1,chID,idxT_cue))-ReferenceLFP1;
                            ReferenceLFP2=squeeze(mean(LFPinfo.ArrayAlign(idx1,AreaCh,idx_array_LFP),2));
                            SpikeLFPInfo.ArrayOnLFP{idate,iarea}{num}=squeeze(LFPinfo.ArrayAlign(idx1,chID,idx_array_LFP))-ReferenceLFP2; 
                            
                            idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            SpikeLFPInfo.CueOffSpike{idate,iarea}{num}=squeeze(SpikeTrain.cueAlign(idx2,ich,idxT_cue));
                            SpikeLFPInfo.ArrayOffSpike{idate,iarea}{num}=squeeze(SpikeTrain.arrayAlign(idx2,ich,idx_array)); 
                            
                            ReferenceLFP3=squeeze(mean(LFPinfo.cueAlign(idx2,AreaCh,idxT_cue),2));
                            SpikeLFPInfo.CueOffLFP{idate,iarea}{num}=squeeze(LFPinfo.cueAlign(idx2,chID,idxT_cue))-ReferenceLFP3;
                            ReferenceLFP4=squeeze(mean(LFPinfo.ArrayAlign(idx2,AreaCh,idx_array_LFP),2));
                            SpikeLFPInfo.ArrayOffLFP{idate,iarea}{num}=squeeze(LFPinfo.ArrayAlign(idx2,chID,idx_array_LFP))-ReferenceLFP4;   

                            SpikeLFPInfo.ElecInfo{idate, iarea}(num,:)=[SpikeTrain.channel(ich), RFOn(ich), SaccadeOn(ich)];
                        end
                    end
                end
            end
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_5sigma.mat'],'SpikeLFPInfo','-v7.3');
        end
    end
end

%% get the power spectrum of LFP by hilbert transform in different areas
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
centerFreq=3:2:55;
for igroup=4
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=1
           load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']);
            count=zeros(1,3);
            for iarea=1:3
               for idate=1:size(SpikeLFPInfo.ElecInfo,1) 
                   LFPElec=[];
                    for icell=1:numel(SpikeLFPInfo.CueOnSpike{idate,iarea})
                        if ~ismember(SpikeLFPInfo.ElecInfo{idate,iarea}(icell,1),LFPElec)
                            LFPElec=[LFPElec, SpikeLFPInfo.ElecInfo{idate,iarea}(icell,1)]; % to make sure no repeated LFPs are counted                            
                            count(iarea)=count(iarea)+1;
                            Vaiables=fieldnames(SpikeLFPInfo);
                            
                            for iVar=6:9
                                data=transpose(SpikeLFPInfo.(Vaiables{iVar}){idate,iarea}{icell});
                                data=data-repmat(mean(data,1),size(data,1),1);
                                num=0;
                                for ifreq=3:2:55
                                    num=num+1;
                                    signal_filtered=bandpass(data, [ifreq-2 ifreq+2], Fs);
                                    signal_hilbert=hilbert(signal_filtered);
                                    power_hb=mean(abs(signal_hilbert),2);
                                    %                                             baseline_hb=mean(power_hb(1:300));
                                    Hilbert_power.(Vaiables{iVar}(1:end-3)){iarea}(count(iarea),num,:)=power_hb;  %(power_hb-baseline_hb)/baseline_hb;
                                end
                            end
                        end
                    end
               end
            end     
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_Hilbert_power_AllChannelReference.mat'],'Hilbert_power','-v7.3');
        end
    end
end

%% plot the hilbert results
close all;
monkeyset={'Mikey','Vasco'};
AreaName={'LIP','PUL','FEF'};
CueTime=-600:1000;
ArrayTime=-500:900;
BandFreq=[16 25];
for imonkey=1:2
    monkey=monkeyset{imonkey};
    load(['Z:\Rujia\Results\' monkey '_exo_AllVisual_Hilbert_power_AllChannelReference.mat']);
    figure;
    for iarea=1:3
        aa=squeeze(mean(Hilbert_power.CueOn{iarea},1));
        bb=squeeze(mean(Hilbert_power.CueOff{iarea},1));
        mm=squeeze(mean(Hilbert_power.ArrayOn{iarea},1));
        nn=squeeze(mean(Hilbert_power.ArrayOff{iarea},1));
        
        subplot(2,3,iarea);
%         plot(CueTime, mean(aa(centerFreq>=BandFreq(1)&centerFreq<=BandFreq(2),:)),'r'); hold on;
%         plot(CueTime, mean(bb(centerFreq>=BandFreq(1)&centerFreq<=BandFreq(2),:)),'b'); hold on;
%         xlim([-300 600]);
        plot(centerFreq, mean(aa(:,CueTime>=200&CueTime<=400),2),'r'); hold on;
        plot(centerFreq, mean(bb(:,CueTime>=200&CueTime<=400),2),'b'); hold on;
                
%         imagesc(CueTime,centerFreq(centerFreq>10),log(aa(centerFreq>10,:)));
%         axis xy;        
%         ylim([10 55]);
        title([ AreaName{iarea} '-Cue']);
        
        subplot(2,3,3+iarea);
%         plot(ArrayTime, mean(mm(centerFreq>=BandFreq(1)&centerFreq<=BandFreq(2),:)),'r'); hold on;
%         plot(ArrayTime, mean(nn(centerFreq>=BandFreq(1)&centerFreq<=BandFreq(2),:)),'b'); hold on;
%         xlim([-200 600]);
        plot(centerFreq, mean(mm(:,ArrayTime>=200&ArrayTime<=400),2),'r'); hold on;
        plot(centerFreq, mean(nn(:,ArrayTime>=200&ArrayTime<=400),2),'b'); hold on;
%         imagesc(CueTime,centerFreq(centerFreq>10),log(bb(centerFreq>10,:)));
%         axis xy;        
%         ylim([10 55]);
        title([AreaName{iarea} '-Array']);
    end
end

%% get the power spectrum of LFP by fast fourier transform (fft) in different areas
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
t=cell(1,4);
for igroup=1:2
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']); %_AllChannelReference
            count=zeros(1,3);
            for iarea=1:3
               for idate=1:size(SpikeLFPInfo.ElecInfo,1) 
                   LFPElec=[];
                   for icell=1:numel(SpikeLFPInfo.CueOnSpike{idate,iarea})
                       if ~ismember(SpikeLFPInfo.ElecInfo{idate,iarea}(icell,1),LFPElec)
                           LFPElec=[LFPElec, SpikeLFPInfo.ElecInfo{idate,iarea}(icell,1)]; % to make sure no repeated LFPs are counted
                           count(iarea)=count(iarea)+1;
                           Variables=fieldnames(SpikeLFPInfo);
                           imonkey, iarea, idate; icell,
                           for iVar=6:9
                               data=transpose(SpikeLFPInfo.(Variables{iVar}){idate,iarea}{icell});
                               data=rmlinesc(data,params,[], 0, 60);
                               data=data-repmat(mean(data,1),size(data,1),1);
                               [S.(Variables{iVar}(1:end-3)){iarea}(count(iarea),:,:),t{iVar-5},f,~]=mtspecgramc(data, movingwin, params);
                           end
                       end
                   end
               end
            end
            
            S.t_cue=t{1}-0.6;
            S.t_array=t{2}-0.5;
            S.f=f;
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_fft_power_demeaned_AllChannelReference_cueResponsive_SU.mat'],'S','-v7.3');
        end
    end
end

%% plot the power spectrum of LFP by fft in different areas
close all;
monkeyset={'Mikey','Vasco'};
AreaName={'LIP','PUL','FEF'};
CueTime=-600:1000;
ArrayTime=-500:900;
BandFreq=[16 25];
climit{1}=[-20 20;-15 15; -30 30]; % for cue on/off condition percent
climit{2}=[-20 20;-5 5; -20 20];

% climit{1}=[-30 30;-20 20; -60 60]; % for array on/off condition percent
% climit{2}=[-30 30;-15 15; -30 30];

% climit=[-5 5; -0.8 0.8; -3 3]/100000;
for igroup=4
    groupName=Group{igroup};
    figure;
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        %%%%%%%% {'_demeaned_AllChannelReference_cueResponsive_SU'}
        load(['Z:\Rujia\Results\' monkey '_exo_' groupName '_fft_power_demeaned_AllChannelReference_cueResponsive_SU.mat']);         
        for iarea=1:3
            idx=S.f>20;
            Variables=fieldnames(S);
            dd=zeros(4,size(S.CueOn{iarea},1));
            for ii=1:4
                pp=squeeze(mean(S.(Variables{ii}){iarea},2));           
                ppf=abs(fft(transpose(pp(:,idx))));  
                dd(ii,:)=ppf(17,:)-ppf(16,:)<=0;
            end
            ValidID=sum(dd,1)==4;
            ivar=1;
            aa=transpose(squeeze(mean(S.(Variables{ivar}){iarea}(ValidID,:,:),1)));
            cc=transpose(squeeze(mean(S.(Variables{ivar+2}){iarea}(ValidID,:,:),1)));
            
            if ivar==1||ivar==3
                Time=S.t_cue; 
                baseline=repmat(mean(aa(:, Time>=-0.3& Time<=-0.05),2),1,size(aa,2));
            else
                Time=S.t_array; 
                bb=transpose(squeeze(mean(S.(Variables{ivar-1}){iarea}(ValidID,:,:),1)));
                baseline=repmat(mean(bb(:, Time>=-0.3& Time<=-0.05),2),1,size(aa,2));
            end            
            aa=(aa-baseline)./baseline*100;
            
            subplot(2,3,iarea+3*imonkey-3);     
            idx=S.f>8;
            idxT=Time>-0.35;            
            imagesc(Time(idxT), S.f(idx), aa(idx,idxT)); hold on;
%             imagesc(Time(idxT), S.f(idx), aa(idx,idxT)-cc(idx,idxT)); hold on;
            axis xy;
            title([AreaName{iarea} '-Cue']);
            if ivar==1||ivar==3
                xlabel('Time to cue onset (s)');
            else
                xlabel('Time to array onset (s)');
            end
            ylabel('Frequency (Hz)');
%             caxis(climit(iarea,:));
            caxis(climit{imonkey}(iarea,:));
            colorbar;           
        end
    end
end

%% get the temporal structure of spike-field coherence within each areas
close all;
clear all;
folder='Z:\Rujia\Results\';
folder1='C:\Data\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
CellTypeName={'NarrowUnits','BroadUnits'};
                                                          
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];

TimeCue=-600:1000;
Nsmooth=50;
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
iCellType=2;
celltype=CellTypeName{iCellType};  
for igroup=4
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};        
        load([folder 'Visual_latency_fitting_' num2str(Nsmooth)  'ms_smooth_exo_cue_' monkey '_' groupName '.mat']);        
        for isexo=0:1
%            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']);  
           load([folder1 monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_DeSpike_linear_cueResponsive_SU.mat']);  
           load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);           
                                                   
           Variables=fieldnames(SpikeLFPInfo);           
           CellType=cell(size(SpikeLFPInfo.ElecInfo,1),3);           
           for iarea=1:3
               CountUnit=0;
               for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                   if ~isempty(SpikeLFPInfo.ElecInfo{idate, iarea})'
                       NN=size(SpikeLFPInfo.ElecInfo{idate, iarea},1);
                       CellType{idate,iarea}=SpikeShape.CellType{iarea}(CountUnit+1:CountUnit+NN);
                       CountUnit=CountUnit+NN;                       
                   end
               end
           end
           count=zeros(3,4);
           for iarea=1:3
               for ivar=6:9
                   Ncell=0;
                   for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                       if ~isempty(SpikeLFPInfo.CueOnSpike{idate,iarea})
                           for icell=1:numel(SpikeLFPInfo.CueOnSpike{idate,iarea})                                                                
                               if CellType{idate,iarea}(icell)==iCellType  
                                   Ncell=Ncell+1;
                                   data1=transpose(SpikeLFPInfo.(Variables{ivar}){idate,iarea}{icell}); % linear noise has been removed for despike LFP
%                                    data1=rmlinesc(data1,params,[], 0, 60);
%                                    data1=data1-repmat(mean(data1,1),size(data1,1),1);

                                   Spike1=transpose(SpikeLFPInfo.(Variables{ivar-4}){idate,iarea}{icell});
                                   ValidTrl1=find(sum(Spike1,1)>=1); %(>=5)
                                   if numel(ValidTrl1)>=15  
                                       count(iarea,ivar-5)=count(iarea,ivar-5)+1;
                                       LFP1=data1(:,ValidTrl1);
                                       [C.(Variables{ivar}(1:end-3)){iarea}(count(iarea,ivar-5),:,:),Phase.(Variables{ivar}(1:end-3)){iarea}(count(iarea,ivar-5),:,:),~,~,~,t{ivar-5},f,~]=cohgramcpb(LFP1,Spike1(:,ValidTrl1),movingwin,params);                                                                                            
                                       Vislatency.(Variables{ivar}(1:end-3)){iarea}(count(iarea,ivar-5),:)=RespLatency{idate,iarea}(:,icell);
                                   end
                               end
                           end
                       end
                   end
               end
           end
           C.t_cue=t{1}-0.6;
           C.t_array=t{2}-0.5;
           C.f=f;          
           save([folder1 monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_despike_linear_SU.mat'],'C','-v7.3');
           save([folder1 monkey '_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_despike_linear_SU.mat'],'Vislatency','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence-phase_demeaned.mat'],'Phase','-v7.3');
        end
    end
end

%% compare the local spike-field coherence between early and late cell groups 
close all;
clear all;
folder='Z:\Rujia\Results\';
CellTypeName={'NarrowUnits','BroadUnits'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
AreaName={'LIP','PUL','FEF'};
LatRange=[0, 80;80, 200];
for iCellType=1
    celltype=CellTypeName{iCellType};
    for isexo=0:1
        figure;
        Ncell=zeros(1,3);
        C1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_SU.mat']);
        C2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_SU.mat']);
        V1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_SU.mat']);
        V2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_SU.mat']);
        for iarea=1:3
            Var='CueOn';  %'ArrayOn';  %
            coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});
            latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
            param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];
            Ncell(iarea)=length(latency);
            t_cue=C1.C.t_cue*1000;
            subplot(2,3,iarea);
            if iarea<=2
                idx=latency<=LatRange(1,2)&latency>LatRange(1,1);
            else
                idx=latency>LatRange(1,1);
            end
            idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
            NN=sum(idx&idx0);
            mm=transpose(squeeze(mean(coh(idx&idx0,:,:),1)));
            if strcmp(Var,'ArrayOn')
                t=C1.C.t_array*1000;
                coh1=cat(1,C1.C.CueOn{iarea},C2.C.CueOn{iarea});
                latency1=[squeeze(V1.Vislatency.CueOn{iarea}(:,1)); squeeze(V2.Vislatency.CueOn{iarea}(:,1))];
                param1=[squeeze(V1.Vislatency.CueOn{iarea}(:,2)); squeeze(V2.Vislatency.CueOn{iarea}(:,2))];
                if iarea==2
                    idx2=latency1<=LatRange(1,2)&latency1>LatRange(1,1);
                else
                    idx2=latency1>LatRange(1,1);
                end
                idx3=~isnan(squeeze(mean(mean(coh1,2),3)))&param1>0.4;
                nn=transpose(squeeze(mean(coh1(idx2&idx3,:,:),1)));
                baseline = repmat(mean(nn(:,t_cue>-200&t_cue<-50),2),1,size(mm,2));
            elseif strcmp(Var,'ArrayOff')
                t=C1.C.t_array*1000;
                coh1=cat(1,C1.C.CueOff{iarea},C2.C.CueOff{iarea});
                latency1=[squeeze(V1.Vislatency.CueOff{iarea}(:,1)); squeeze(V2.Vislatency.CueOff{iarea}(:,1))];
                param1=[squeeze(V1.Vislatency.CueOff{iarea}(:,2)); squeeze(V2.Vislatency.CueOff{iarea}(:,2))];
                if iarea==2
                    idx2=latency1<=LatRange(1,2)&latency1>LatRange(1,1);
                else
                    idx2=latency1>LatRange(1,1);
                end
                idx3=~isnan(squeeze(mean(mean(coh1,2),3)))&param1>0.4;
                nn=transpose(squeeze(mean(coh1(idx2&idx3,:,:),1)));
                baseline = repmat(mean(nn(:,t_cue>-200&t_cue<-50),2),1,size(mm,2));                                                
            else
                t=t_cue;
                baseline=repmat(mean(mm(:,t>-200&t<-50),2),1,size(mm,2));
            end
            mm=(mm-baseline)./baseline*100;
            imagesc(t,C1.C.f,mm); hold on;
            line([0,0],[0 100],'color','k','linestyle','--'); hold on;
            line([-200, 500],[15, 15],'color','k','linestyle','--'); hold on;
            line([-200, 500],[30, 30],'color','k','linestyle','--'); hold on;
            axis xy;
            xlim([-200, 500]);
            title([AreaName{iarea} ':' num2str(NN)]);
            caxis([-30,30]);
            colorbar;
            
            if iarea==2
                subplot(2,3,iarea+3);
                idx=latency>LatRange(2,1)&latency<=LatRange(2,2);
                NN=sum(idx&idx0);
                mm=transpose(squeeze(mean(coh(idx&idx0,:,:),1)));
                if strcmp(Var,'ArrayOn')||strcmp(Var,'ArrayOff')
                    t=C1.C.t_array*1000;
                    idx2=latency1>100;
                    nn=transpose(squeeze(mean(coh1(idx2&idx3,:,:),1)));
                    baseline = repmat(mean(nn(:,t_cue>-200&t_cue<-50),2),1,size(mm,2));
                else
                    t=t_cue;
                    baseline=repmat(mean(mm(:,t>-200&t<-50),2),1,size(mm,2));
                end                
                mm=(mm-baseline)./baseline*100;
                imagesc(t,C1.C.f,mm); hold on;
                line([0,0],[0 100],'color','k','linewidth',1,'linestyle','--'); hold on;
                line([-200, 500],[15, 15],'color','k','linestyle','--'); hold on;
                line([-200, 500],[30, 30],'color','k','linestyle','--'); hold on;
                axis xy;
                xlim([-200, 500]);
                title([AreaName{iarea} ':' num2str(NN)]);
                caxis([-30, 30]);
                colorbar;
            end
            %     subplot(1,3,iarea);
            %     hist(latency); hold on;
            %     line([70, 70],[0 25]); hold on;
            
        end
    end
end
%% plot the temporal structure of results of spike-field coherence within the same area for each monkey
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
t_cue=-0.4990:0.02:0.9010;
t_array=-0.3990:0.02:0.8010;

for igroup=1:2
    groupName=Group{igroup};
    for isexo=0:1
        figure;
        for imonkey=1:2
            monkey=monkeyset{imonkey};            
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_AllChannelReference_cueResponsive_SU.mat']);
            Vars=fieldnames(C);
    
            for iarea=1:3
                for ivar=2
                    idx1=~isnan(squeeze(sum(sum(C.(Vars{ivar}){iarea},2),3)));
                    aa=transpose(squeeze(mean(C.(Vars{ivar}){iarea}(idx1,:,:),1)));
                    
                    idx2=~isnan(squeeze(sum(sum(C.(Vars{ivar+2}){iarea},2),3)));
                    cc=transpose(squeeze(mean(C.(Vars{ivar+2}){iarea}(idx2,:,:),1)));
                    if ivar==1||ivar==3
                        Time=C.t_cue;
                        baseline=repmat(mean(aa(:, Time>=-0.3& Time<=-0.05),2),1,size(aa,2));
                        colorlim=[-30 30; -20 20; -40 40]; % for cue on/off
                    else
                        Time=C.t_array;
                        idx2=~isnan(squeeze(sum(sum(C.(Vars{ivar-1}){iarea},2),3)));
                        bb=transpose(squeeze(mean(C.(Vars{ivar-1}){iarea}(idx2,:,:),1)));
                        baseline=repmat(mean(bb(:, Time>=-0.3& Time<=-0.05),2),1,size(aa,2));
                        colorlim=[-30 30; -20 20; -30 30]; % for array on/off
                    end
%                     aa=(aa-baseline)./baseline*100;

                    subplot(2,3,(imonkey-1)*3+iarea);
                    imagesc(Time(Time>=-0.3), C.f, aa(:,Time>=-0.3)-cc(:,Time>=-0.3)); hold on;
                    if ivar==1||ivar==3
                        xlabel('Time to cue onset (s)');
                    else
                        xlabel('Time to array onset (s)');
                    end
                    ylabel('Frequency (Hz)');
                    line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
%                     caxis(colorlim(iarea,:));
                    axis xy;
                    title([AreaName{iarea} '-' Vars{ivar}]);
                    colorbar;
                end
            end
        end
    end
end

%% converge all SFC pattern WITHIN THE SAME AREA from two monkeys together
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
CellTypeName={'NarrowUnits','BroadUnits'};
for igroup=4
    groupName=Group{igroup};
    for icelltype=2
    figure;
    for isexo=0
%         C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_AllChannelReference_cueResponsive_SU.mat']);  %demeaned
%         C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_AllChannelReference_cueResponsive_SU.mat']);
        C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']); 
        C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  
        Vars=fieldnames(C1.C);
        Time1=C1.C.t_cue;
        Time2=C1.C.t_array;
        for iarea=2
            ivar=1;
%             data1=cat(1, C1.C.(Vars{ivar}){iarea}, C2.C.(Vars{ivar}){iarea}); 
            data1=cat(1, C1.C.(Vars{ivar}){iarea, iarea}, C2.C.(Vars{ivar}){iarea,iarea});% for celltype-specific SFC only
            idx1=~isnan(squeeze(sum(sum(data1,2),3)));
            aa=transpose(squeeze(mean(data1(idx1,:,:),1)));
            %%%%%% to remove the baseline(coherence before cue onset)
            if ivar==1||ivar==3
                Time=C1.C.t_cue;
                baseline=repmat(mean(aa(:, Time1>=-0.25& Time1<= -0.1),2),1,size(aa,2));
                colorlim=[-15 15; -15 15; -20 20]; % for cue on/off
            else
                Time=C1.C.t_array;
                data2=cat(1, C1.C.(Vars{ivar-1}){iarea}, C2.C.(Vars{ivar-1}){iarea});
                idx2=~isnan(squeeze(sum(sum(data2,2),3)));
                bb=transpose(squeeze(mean(data2(idx2,:,:),1)));
                baseline=repmat(mean(bb(:, Time1>=-0.25& Time1<= -0.1),2),1,size(aa,2));
                colorlim=[-30 30; -20 20; -30 30]; % for array on/off
            end
            aa=(aa-baseline)./baseline*100;
            kernal=ones(4)/16;
            aa=conv2(aa,kernal,'same');
            
            %%%%%%%%%%%%%%% compare cue/array on vs off
            data2=cat(1, C1.C.(Vars{ivar+2}){iarea}, C2.C.(Vars{ivar+2}){iarea});
            idx2=~isnan(squeeze(sum(sum(data2,2),3)));
            cc=transpose(squeeze(mean(data2(idx2,:,:),1)));
            
            subplot(2,3,isexo*3+iarea);
            imagesc(Time(Time>=-0.35), C1.C.f, aa(:,Time>=-0.35)); hold on;  %-cc(:,Time>=-0.35)
            
            
            if ivar==1||ivar==3
                xlabel('Time to cue onset (s)');
            else
                xlabel('Time to array onset (s)');
            end
            ylabel('Frequency (Hz)');
            line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
            line([-0.4 0.7], [15 15],'color','k'); hold on;
            line([-0.4 0.7], [25 25],'color','k'); hold on;
            xlim([-0.35 0.6]);
            ylim([2 100]);
            caxis(colorlim(iarea,:));
            axis xy;
            title([AreaName{iarea} '-' Vars{ivar}]);
            colorbar;
        end
    end
    end
end

%%
%%%%%% plot the frequency averaged of SFC WITHIN THE SAME regions,
%%%%%% iarea==jarea
close all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
FigOrder=[3 1 5 4 2 6];
clr2=[1 0.4 0; 0 0 0];
clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
mmSFC=cell(1,2);
CellTypeName={'NarrowUnits','BroadUnits'};
% figure;
% num=0;
ll={'-','--'};

for igroup=1
    groupName=Group{igroup};
    for isexo=0:1
        figure;
        num=0;
        for icelltype=1
            C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %
            C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %AttentionModulated
            Vars=fieldnames(C1.C);
            Time1=C1.C.t_cue;
            Time2=C1.C.t_array;
            f=C1.C.f;
            for iarea=1:3
                clr2=[clr3(iarea,:); 0 0 0];
                data=cell(1,2); mSFC=cell(1,2); ee=cell(1,2); SFC_f=cell(1,2);
                data{1}=cat(1, C1.C.CueOn{iarea,iarea}, C2.C.CueOn{iarea,iarea});
                data{2}=cat(1, C1.C.ArrayOn{iarea,iarea}, C2.C.ArrayOn{iarea,iarea});
                ss2Baseline=cell(1,2); % significance test compared to baseline
                if ~isempty(data{1})
                    idxF=f>=15&f<=30;
                    for ii=1:2
                        idx1=~isnan(squeeze(sum(sum(data{ii},2),3)));
                        SFC_f{ii}=squeeze(mean(data{ii}(idx1,:,idxF),3));
                        if ii<=1
                            sBaseline=mean(SFC_f{1}(:,Time1<-0.1&Time1>-0.25),2);
                            baseline=mean(sBaseline,1);
                        end
                        for itime=1:size(SFC_f{ii},2)
                            ss2Baseline{ii}(itime)=ttest2(SFC_f{ii}(:,itime), sBaseline);
                        end
                        %                             SFC_f{ii}=SFC_f{ii}-repmat(baseline{2-mod(ii,2)},size(SFC_f{ii},1),size(SFC_f{ii},2));
                        SFC_f{ii}=(SFC_f{ii}-baseline)/baseline*100;
                        mSFC{ii}=mean(SFC_f{ii},1);
                        ee{ii}=std(SFC_f{ii},0,1)/sqrt(sum(idx1));
                    end
                    for ialign=1:2
                        if ialign==1
                            Time=Time1;
                        else
                            Time=Time2;
                        end                        
                        subplot(2,3,iarea+(ialign-1)*3);
                        pp(1)=patchplot(Time, mSFC{ialign}, ee{ialign}, clr2(icelltype,:)); hold on;
                        ss2=ss2Baseline{ialign};
                        plot(Time(ss2==1),ss2(ss2==1)*-5*icelltype,'*','color', clr2(icelltype,:)); hold on;
                        plot([0 0],[-15 40],'--k'); hold on;
                        plot([-0.3 0.6],[0 0],'k'); hold on;
                        xlim([-0.3 0.5]);
                        %                         ylim([-15 15]);
                        yticks([-15 0 15 30 45]);                        
                        xticks(-0.2:0.2:0.6);
%                         yticklabels({'','',''});
%                         xticklabels({'','','',''});
                        set(gca,'box','off','linewidth',2.0,'fontsize',15);
                    end
                end
            end
        end
    end
end

%% plot the temporal structure of spike-field coherence within the same area for shuffled data
AreaName={'LIP','PUL','FEF'};
figure;
colorlim=[-10 20; -10 30; -20 60];
clear mCohe_Shuffled;
for igroup=4
    groupName=Group{igroup};
    OnOffDiff=cell(1,3);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_demeaned_shuffle.mat']); 
            for iarea=1:3                
                idx1=~isnan(squeeze(sum(sum(sum(C.CueOn{iarea},2),3),4)));
                idx2=~isnan(squeeze(sum(sum(sum(C.CueOff{iarea},2),3),4)));
                idx3=~isnan(squeeze(sum(sum(sum(C.ArrayOn{iarea},2),3),4)));
                idx4=~isnan(squeeze(sum(sum(sum(C.ArrayOff{iarea},2),3),4)));                
                
                aa=squeeze(mean(C.CueOn{iarea}(idx1,:,:,:),2));
                ee1=squeeze(std(C.CueOn{iarea}(idx1,:,:,:),0,2));
                mCohe_Shuffled.CueOn{iarea} =aa;                
                mCohe_Shuffled.eeCueOn{iarea} =ee1;
                
                bb=squeeze(mean(C.CueOff{iarea}(idx2,:,:,:),2));
                ee2=squeeze(std(C.CueOff{iarea}(idx2,:,:,:),0,2));
                mCohe_Shuffled.CueOff{iarea} =bb;
                mCohe_Shuffled.eeCueOff{iarea} =ee2;
                
                mm=squeeze(mean(C.ArrayOn{iarea}(idx3,:,:,:),2));
                ee3=squeeze(std(C.ArrayOn{iarea}(idx3,:,:,:),0,2));
                mCohe_Shuffled.ArrayOn{iarea} =mm;
                mCohe_Shuffled.eeArrayOn{iarea} =ee3;
                
                nn=squeeze(mean(C.ArrayOff{iarea}(idx4,:,:,:),2));
                ee4=squeeze(std(C.ArrayOff{iarea}(idx4,:,:,:),0,2));
                mCohe_Shuffled.ArrayOff{iarea} =nn;
                mCohe_Shuffled.eeArrayOff{iarea} =ee4;
%                 OnOffDiff{iarea}(imonkey,:,:) =  mm-nn;
%                 
%                 subplot(2,3,(imonkey-1)*3+iarea);
% %                 imagesc(C.t_cue(C.t_cue>=-0.3), C.f, aa(:,C.t_cue>=-0.3)); hold on;
% %                 xlabel('Time to cue onset (s)');                
%                 imagesc(C.t_array, C.f, aa); hold on;
%                 xlabel('Time to array onset (s)');                  
%                 ylabel('Frequency (Hz)');
% %                 caxis(colorlim(iarea,:));
% %                 caxis([-0.015 0.015]);
%                 axis xy;
%                 title([AreaName{iarea} '-On']);
%                 colorbar;
            end  
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_mean_coherence_shuffled.mat'],'mCohe_Shuffled','-v7.3');            
        end
    end
end



%% get the spike-field coherence by fft between different area
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 200];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
T_cue=-600:1000;
T_array=-500:900;
Nsmooth=50;
SpkRange=[5 10; 12 25];
CellTypeName={'NarrowUnits','BroadUnits'};
for icelltype=1
celltype=CellTypeName{icelltype};
for igroup=4
    groupName=Group{igroup};
    for imonkey=2
        monkey=monkeyset{imonkey};
        for isexo=1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']); 
%             load([folder1 monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_DeSpike_linear_cueResponsive_SU.mat']);  
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_cueResponsive_SU.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_cueResponsive_SU.mat']);
            load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
            load([folder 'Visual_latency_fitting_' num2str(Nsmooth)  'ms_smooth_exo_cue_' monkey '_' groupName '.mat']);
        
            PeakLoc=cell(size(SpikeLFPInfo.ElecInfo,1),3);
            WaveAmp=cell(size(SpikeLFPInfo.ElecInfo,1),3);
            CellType = cell(size(SpikeLFPInfo.ElecInfo,1),3);
            pCueDelay=cell(size(SpikeLFPInfo.ElecInfo,1),3);
            pArrayDelay=cell(size(SpikeLFPInfo.ElecInfo,1),3);
            for iarea=1:3
                CountUnit=0;
                for isession=1:size(SpikeLFPInfo.ElecInfo,1)
                    if ~isempty(SpikeLFPInfo.ElecInfo{isession, iarea})
                        NN=size(SpikeLFPInfo.ElecInfo{isession, iarea},1);
%                         PeakLoc{isession,iarea}=SpikeShape.PeakRight{iarea}(CountUnit+1:CountUnit+NN);
%                         WaveAmp{isession,iarea}=SpikeShape.ValueRight{iarea}(CountUnit+1:CountUnit+NN);                        
                        CellType{isession,iarea}=SpikeShape.CellType{iarea}(CountUnit+1:CountUnit+NN);
                        pCueDelay{isession,iarea}=pValueCueDelay{iarea}(CountUnit+1:CountUnit+NN);
                        pArrayDelay{isession,iarea}=pValueArrayDelay{iarea}(CountUnit+1:CountUnit+NN);
                        CountUnit=CountUnit+NN;                        
                    end                    
                end
            end

%             clear C;
%             clear Phase;
%             clear PairDate;
            Vars=fieldnames(SpikeLFPInfo);           
            for ivar=6:9
                count=zeros(3);  
                if ivar==6||ivar==8
                    idxT=T_cue>=300&T_cue<=800;
                elseif ivar==7||ivar==9
                    idxT=T_array>=250&T_array<=500;
                end
                for iarea=1:3                       
                    for jarea=1:3                        
%                         if imonkey==1
%                             StartDate=4;
%                         else
%                             StartDate=1;
%                         end
                        for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                            ivar, idate,                            
                            if ~isempty(SpikeLFPInfo.(Vars{ivar}){idate,iarea})&&~isempty(SpikeLFPInfo.(Vars{ivar}){idate,jarea})
                                for icell=1:numel(SpikeLFPInfo.(Vars{ivar}){idate,iarea})                                    
%                                     Cellidx1= PeakLoc{idate,iarea}(icell)>= SpkRange(icelltype,1)&PeakLoc{idate,iarea}(icell)<= SpkRange(icelltype,2);  % & WaveAmp{idate,iarea}(icell)>=0.3 ;                                      
                                    if CellType{idate,iarea}(icell)==icelltype  %pCueDelay{idate,iarea}(icell)<0.05||pArrayDelay{idate,iarea}(icell)<0.05
%                                         LFP1=transpose(SpikeLFPInfo.(Vars{ivar}){idate,iarea}{icell});
%                                         LFP1=rmlinesc(LFP1,params,[], 0, 60);
%                                         LFP1=LFP1-repmat(mean(LFP1,1),size(LFP1,1),1);
                                        Spike1=transpose(SpikeLFPInfo.(Vars{ivar-4}){idate,iarea}{icell});                                                                               
                                        for jcell=1:numel(SpikeLFPInfo.(Vars{ivar}){idate,jarea})
%                                             Cellidx2= PeakLoc{idate,jarea}(jcell)>= SpkRange(icelltype,1)&PeakLoc{idate,jarea}(jcell)<= SpkRange(icelltype,2) & WaveAmp{idate,jarea}(jcell)>=0.3 ;   
                                            if SpikeLFPInfo.ElecInfo{idate,iarea}(icell,2)== SpikeLFPInfo.ElecInfo{idate,jarea}(jcell,2) %&& Cellidx2   %size(Spike1,2)==size(SpikeLFPInfo.(Vars{ivar-4}){idate,jarea}{jcell},1)
                                                LFP2=transpose(SpikeLFPInfo.(Vars{ivar}){idate,jarea}{jcell});
                                                LFP2=rmlinesc(LFP2,params,[], 0, 60);
                                                LFP2=LFP2-repmat(mean(LFP2,1),size(LFP2,1),1);
                                                Spike2=transpose(SpikeLFPInfo.(Vars{ivar-4}){idate,jarea}{jcell});
                                                
                                                ValidTrl=sum(Spike1,1)>=1&sum(Spike2,1)>=1;
                                                if sum(ValidTrl)>=15
                                                    count(iarea,jarea)=count(iarea,jarea)+1;         
%                                                     PairDate.(Vars{ivar}(1:end-3)){iarea,jarea}(count(iarea,jarea))=idate;                                                
                                                    Vislatency.(Vars{ivar}(1:end-3)){iarea,jarea}(count(iarea,jarea),:)=RespLatency{idate,iarea}(:,icell);
                                                    [C.(Vars{ivar}(1:end-3)){iarea,jarea}(count(iarea,jarea),:,:),~,~,~,~,t{ivar-5},f,~]=cohgramcpb(LFP2,Spike1,movingwin,params);
                                                   
                                                    %[C.(Vars{ivar}(1:end-3)){jarea,iarea}(count(jarea,iarea),:,:),~,~,~,~,t{ivar-5},f,~]=cohgramcpb(LFP1,Spike2,movingwin,params);
                                                    
%                                                     [C.(Vars{ivar}(1:end-3)){iarea,jarea}(count(iarea,jarea),:),~,~,~,~,~]=coherencycpb(LFP2(idxT,:),Spike1(idxT,:),params);
%                                                     [C.(Vars{ivar}(1:end-3)){jarea,iarea}(count(jarea,iarea),:),~,~,~,~,f]=coherencycpb(LFP1(idxT,:),Spike2(idxT,:),params);
                                                end
                                            end
                                        end
                                    end
                                end                                
                            end
                        end                                                                  
                    end
                end
            end            
            C.f=f;           
%             C.cue_delay=[300 800];
%             C.array_delay=[250 500];
            C.t_cue=t{1}-0.6;
            C.t_array=t{2}-0.5;  
%              save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_PairDate_' CellTypeName{icelltype} '_SU.mat'],'PairDate','-v7.3'); 
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_Baseline_AttentionModulated_SU.mat'],'C_base','-v7.3');
%              save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_Delay_cueResponsive_SU.mat'],'C','-v7.3');
%              save([folder 'f_spike-field-coherence_Delay.mat'],'f','-v7.3');
         
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat'],'Vislatency','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat'],'C','-v7.3');  %cueResponsive
        end
    end
end
end

%% plot the spike_field_coherence between different areas for each individual monkey
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
t_cue=-0.4990:0.02:0.9010;
t_array=-0.3990:0.02:0.8010;
for igroup=4
    groupName=Group{igroup};
%     figure;
    for imonkey=2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_SU.mat']); %cueResponsive_5sigma  AttentionModulated
            Vars=fieldnames(C);
            ivar=1;
            figure;
            for iarea=1:3
                for jarea=1:3  
                    if ~isempty(C.(Vars{ivar}){iarea,jarea})
                        idx1=~isnan(squeeze(sum(sum(C.(Vars{ivar}){iarea,jarea},2),3)));
                        aa=transpose(squeeze(mean(C.(Vars{ivar}){iarea,jarea}(idx1,:,:),1)));
                        if ivar<=2
                            idx1=~isnan(squeeze(sum(sum(C.(Vars{ivar+2}){iarea,jarea},2),3)));
                            cc=transpose(squeeze(mean(C.(Vars{ivar+2}){iarea,jarea}(idx1,:,:),1)));
                        end
                        
                        if ivar==1||ivar==3
                            Time=C.t_cue;
                            baseline_aa=repmat(mean(aa(:, Time>=-0.3& Time<=-0.1),2),1,size(aa,2));
                            if ivar==1
                                baseline_cc=repmat(mean(cc(:, Time>=-0.3& Time<=-0.1),2),1,size(cc,2));
                            end
                        else
                            Time=C.t_array;
                            idx2=~isnan(squeeze(sum(sum(C.(Vars{ivar-1}){iarea,jarea},2),3)));
                            bb=transpose(squeeze(mean(C.(Vars{ivar-1}){iarea,jarea}(idx2,:,:),1)));
                            baseline_aa=repmat(mean(bb(:, Time>=-0.3& Time<=-0.1),2),1,size(aa,2));
                            if ivar==2
                                idx2=~isnan(squeeze(sum(sum(C.(Vars{ivar+1}){iarea,jarea},2),3)));
                                bb=transpose(squeeze(mean(C.(Vars{ivar+1}){iarea,jarea}(idx2,:,:),1)));
                                baseline_cc=repmat(mean(bb(:, Time>=-0.3& Time<=-0.1),2),1,size(cc,2));
                            end
                            
                        end
                        aa=(aa-baseline_aa)./baseline_aa*100;
                        cc=(cc-baseline_cc)./baseline_cc*100;
                        aa=conv2(aa,ones(4)/16,'same');
                        cc=conv2(cc,ones(4)/16,'same');
                        
                        subplot(3,3,3*(iarea-1)+jarea);
                        imagesc(Time(Time<=0.7), C.f, aa(:,Time<=0.7)); hold on;  %-cc(:,Time>=-0.3)
                        if ivar==1||ivar==3
                            xlabel('Time to cue onset (s)');
                        else
                            xlabel('Time to array onset (s)');
                        end
                        ylabel('Frequency (Hz)');
                        line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
                        caxis([-12 12]);
                        ylim([5 80]);
                        xlim([-0.3 0.7]);
                        axis xy;
                        title([AreaName{iarea} '-' AreaName{jarea}]);
                        colorbar;
                    end
                end
            end
        end
    end
end

%% converge the spike_field_coherence between different areas of two monkeys together
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\Rujia\Results\';
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
LatencyRange=[20, 100; 100, 200];
CellTypeName={'NarrowUnits','BroadUnits'};
for icelltype=1
    celltype=CellTypeName{icelltype};
    Nunit=zeros(2,2,2);
%%%%% plot the temporal-frequency structure of SFC    
    for isexo=0:1
        figure;
        C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);  %cueResponsive
        C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);     %AttentionModulated  
        V1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);  
        V2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);    
        
  
        Vars=fieldnames(C1.C);
        ivar=2;
        Time1=C1.C.t_cue;
        Time2=C1.C.t_array;                
        for iarea=3
            for jarea=2  %(iarea+1):3
                data1=cat(1, C1.C.(Vars{ivar}){iarea,jarea}, C2.C.(Vars{ivar}){iarea,jarea}); 
                data2=cat(1, C1.C.(Vars{ivar}){jarea,iarea}, C2.C.(Vars{ivar}){jarea,iarea});
                param1 = [V1.Vislatency.(Vars{ivar}){iarea,jarea}(:,2);V2.Vislatency.(Vars{ivar}){iarea,jarea}(:,2)];
                Latency1 = [V1.Vislatency.(Vars{ivar}){iarea,jarea}(:,1);V2.Vislatency.(Vars{ivar}){iarea,jarea}(:,1)];
                
                param2 = [V1.Vislatency.(Vars{ivar}){jarea,iarea}(:,2);V2.Vislatency.(Vars{ivar}){jarea,iarea}(:,2)];
                Latency2 = [V1.Vislatency.(Vars{ivar}){jarea,iarea}(:,1);V2.Vislatency.(Vars{ivar}){jarea,iarea}(:,1)];
                
                if ~isempty(data1)
                    for irange=1:2
                        idx0=param1>0.4&Latency1>=20;%&Latency1<=LatencyRange(irange,2)&
                        idx1=~isnan(squeeze(sum(sum(data1,2),3)))&idx0;
                        Nunit(isexo+1,irange,1)=sum(idx1);
                        aa=transpose(squeeze(mean(data1(idx1,:,:),1)));
                        
                        idx0=param2>0.4&Latency2<=LatencyRange(irange,2)&Latency2>LatencyRange(irange,1);
                        idx1=~isnan(squeeze(sum(sum(data2,2),3)))&idx0;
                        Nunit(isexo+1,irange,2)=sum(idx1);
                        cc=transpose(squeeze(mean(data2(idx1,:,:),1)));
                        if ivar==1||ivar==3
                            Time=C1.C.t_cue;
                            baseline_aa=repmat(mean(aa(:, Time1>=-0.25& Time1<=-0.1),2),1,size(aa,2));
                            baseline_cc=repmat(mean(cc(:, Time1>=-0.25& Time1<=-0.1),2),1,size(cc,2));
                        else
                            Time=C1.C.t_array;
                            data3=cat(1, C1.C.(Vars{ivar-1}){iarea,jarea}, C2.C.(Vars{ivar-1}){iarea,jarea});                            
                            param3 = [V1.Vislatency.(Vars{ivar-1}){iarea,jarea}(:,2);V2.Vislatency.(Vars{ivar-1}){iarea,jarea}(:,2)];
                            Latency3 = [V1.Vislatency.(Vars{ivar-1}){iarea,jarea}(:,1);V2.Vislatency.(Vars{ivar-1}){iarea,jarea}(:,1)];
                            idx0=param3>0.4&Latency3>20;  %Latency3<=LatencyRange(irange,2)&
                            idx2=~isnan(squeeze(sum(sum(data3,2),3)))&idx0;
                            bb=transpose(squeeze(mean(data3(idx2,:,:),1)));
                            baseline_aa=repmat(mean(bb(:, Time1>=-0.25& Time1<=-0.1),2),1,size(aa,2));
                                                        
                            data4=cat(1, C1.C.(Vars{ivar-1}){jarea,iarea}, C2.C.(Vars{ivar-1}){jarea,iarea});
                            param4 = [V1.Vislatency.(Vars{ivar-1}){jarea,iarea}(:,2);V2.Vislatency.(Vars{ivar-1}){jarea,iarea}(:,2)];
                            Latency4 = [V1.Vislatency.(Vars{ivar-1}){jarea,iarea}(:,1);V2.Vislatency.(Vars{ivar-1}){jarea,iarea}(:,1)];
                            idx0=param4>0.4&Latency4<=LatencyRange(irange,2)&Latency4>LatencyRange(irange,1);
                            idx3=~isnan(squeeze(sum(sum(data4,2),3)))&idx0;
                            dd=transpose(squeeze(mean(data4(idx3,:,:),1)));
                            baseline_cc=repmat(mean(dd(:, Time1>=-0.25& Time1<=-0.1),2),1,size(cc,2));
                        end
                        aa=(aa-baseline_aa)./baseline_aa*100;
                        mm{1}=conv2(aa,ones(4)/16,'same');
                        cc=(cc-baseline_cc)./baseline_cc*100;
                        mm{2}=conv2(cc,ones(4)/16,'same');
                        Direction={[AreaName{iarea} '->' AreaName{jarea}],[AreaName{jarea} '->' AreaName{iarea}]};
                      
%                         subplot(2,2,1+isexo*2);
                        for idir=1:2
                            subplot(2,2,(irange-1)*2+idir);
                            idxT=Time>=-0.3&Time<=0.5;  %
                            imagesc(Time(idxT), C1.C.f, mm{idir}(:,idxT)); hold on;  % -cc(:,idxT)
                            line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
                            line([-0.3 0.5], [15 15],'color','k','linestyle','--' ); hold on;
                            line([-0.3 0.5], [30 30],'color','k','linestyle','--' ); hold on;
                            caxis([-20 20]);
                            ylim([5 100]);
                            axis xy;
                            title(Direction{idir});
                            set(gca,'box','off','linewidth',2.0,'fontsize',15);
                            colorbar;
                        end
                    end                   
                end
            end
        end
    end
end

%%
%%%%%% plot the time averaged of SFC
% close all;
% for igroup=1:2
%     groupName=Group{igroup};
%     num=0;
%     figure;
%     FigOrder=[3 1 5 4 2 6];
%     clr2=[1 0.4 0; 0 0 0];
%     for isexo=0:1
%         C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_SU.mat']);  %
%         C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_SU.mat']);     %AttentionModulated  
%         Vars=fieldnames(C1.C);
%         ivar=1;
%         if ivar==1||ivar==3
%             Time=C1.C.t_cue;            
%         else
%             Time=C1.C.t_array;         
%         end                
%         for iarea=1:2
%             for jarea=(iarea+1):3
%                 data1=cat(1, C1.C.(Vars{ivar}){iarea,jarea}, C2.C.(Vars{ivar}){iarea,jarea}); 
%                 data2=cat(1, C1.C.(Vars{ivar}){jarea,iarea}, C2.C.(Vars{ivar}){jarea,iarea}); 
%                 if ~isempty(data1)
%                     idxT=Time>=0.4 &Time<=0.6;
%                     idx1=~isnan(squeeze(sum(sum(data1,2),3)));
%                     aa=squeeze(mean(data1(idx1,idxT,:),2));                   
%                     
%                     idx2=~isnan(squeeze(sum(sum(data2,2),3)));
%                     dd=squeeze(mean(data2(idx2,idxT,:),2));                    
%                     
%                     BsTime=Time>=-0.3& Time<=-0.1;
%                     if ivar==1||ivar==3                                
%                         baseline_aa=squeeze(mean(mean(data1(idx1,BsTime,:),1),2)); 
%                         baseline_dd=squeeze(mean(mean(data2(idx2,BsTime,:),1),2));
%                     else                     
%                         data3=cat(1, C1.C.(Vars{ivar-1}){iarea,jarea}, C2.C.(Vars{ivar-1}){iarea,jarea});
%                         idx3=~isnan(squeeze(sum(sum(data3,2),3)));                      
%                         baseline_aa=squeeze(mean(mean(data3(idx3,BsTime,:),1),2)); 
%                         
%                         data4=cat(1, C1.C.(Vars{ivar-1}){jarea,iarea}, C2.C.(Vars{ivar-1}){jarea,iarea});
%                         idx4=~isnan(squeeze(sum(sum(data4,2),3)));                      
%                         baseline_dd=squeeze(mean(mean(data4(idx4,BsTime,:),1),2));    
%                     end
%                     
%                     aa=(aa-repmat(baseline_aa',size(aa,1),1))./repmat(baseline_aa',size(aa,1),1)*100;
%                     dd=(dd-repmat(baseline_dd',size(dd,1),1))./repmat(baseline_dd',size(dd,1),1)*100;
%                     mSFC1=mean(aa,1);
%                     ee1=std(aa,0,1)/sqrt(sum(idx1));
%                     mSFC2=mean(dd,1);
%                     ee2=std(dd,0,1)/sqrt(sum(idx2));
%                   
%                     ss=zeros(1,size(aa,2));
%                     pvalue=zeros(1,size(aa,2));
%                     for ifreq=1:size(aa,2)
%                         [ss(ifreq),pvalue(ifreq)]=ttest2(aa(:,ifreq), dd(:,ifreq));
%                     end
%                                           
%                     num=num+1;
%                     subplot(3,2,FigOrder(num));
%                     f=C1.C.f;
%                     pp(1)=patchplot(f, mSFC1, ee1, clr2(1,:)); hold on;
%                     pp(2)=patchplot(f, mSFC2, ee2, clr2(2,:)); hold on;
%                     plot(f(ss==1),ss(ss==1)*0.065,'k*'); hold on;
%                     plot([0 100],[0 0],'k--'); hold on;
%                     
%                     xlim([10 70]);
% %                     ylim([0.06 0.09]);
% %                     ylim([-15 20]);
% %                     title([AreaName{iarea} '-' AreaName{jarea}]);
%                     if num>3
%                         ll=legend(pp,{[AreaName{iarea} '->' AreaName{jarea}],[AreaName{jarea} '->' AreaName{iarea}]});
%                         set(ll,'box','off');
%                     end
%                     set(gca,'box','off','linewidth',2.0,'fontsize',20);
%                 end
%             end
%         end
%     end
% end

%%
%%%%%% plot the frequency averaged of SFC between each two regions
close all;
FigOrder=[3 1 5 4 2 6];
clr2=[1 0.4 0; 0 0 0];
clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
mmSFC=cell(1,2);
CellTypeName={'NarrowUnits','BroadUnits'};

for icelltype=1:2
    figure;
    num=0;
    for isexo=0:1
        %     figure;
        %     num=0;
        for igroup=2
            groupName=Group{igroup};
            C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %_HighAmp
            C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %AttentionModulated
%             C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_5sigma.mat']);  %_HighAmp
%             C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_5sigma.mat']);  %AttentionModulated
            Vars=fieldnames(C1.C);
            Time1=C1.C.t_cue;
            Time2=C1.C.t_array;
            f=C1.C.f;
            for iarea=1 %:2
                for jarea=3 %:3
                    data=cell(1,4); mSFC=cell(1,4); ee=cell(1,4); SFC_f=cell(1,4);
                    data{1}=cat(1, C1.C.CueOn{iarea,jarea}, C2.C.CueOn{iarea,jarea});
                    data{2}=cat(1, C1.C.CueOn{jarea,iarea}, C2.C.CueOn{jarea,iarea});
                    data{3}=cat(1, C1.C.ArrayOn{iarea,jarea}, C2.C.ArrayOn{iarea,jarea});
                    data{4}=cat(1, C1.C.ArrayOn{jarea,iarea}, C2.C.ArrayOn{jarea,iarea});
                    ss2Baseline=cell(1,4);
                    if ~isempty(data{1})
                        idxF=f>=15&f<=30;
                        sBaseline=cell(1,2);
                        baseline=cell(1,2);
                        for ii=1:4
                            idx1=~isnan(squeeze(sum(sum(data{ii},2),3)));
                            SFC_f{ii}=squeeze(mean(data{ii}(idx1,:,idxF),3));
                            if ii<=2
                                sBaseline{ii}=mean(SFC_f{ii}(:,Time1<-0.1&Time1>-0.25),2);
                                baseline{ii}=mean(sBaseline{ii},1);
                            end
                            for itime=1:size(SFC_f{ii},2)
                                ss2Baseline{ii}(itime)=ttest2(SFC_f{ii}(:,itime), sBaseline{2-mod(ii,2)});
                            end
                            %                             SFC_f{ii}=SFC_f{ii}-repmat(baseline{2-mod(ii,2)},size(SFC_f{ii},1),size(SFC_f{ii},2));
                            SFC_f{ii}=(SFC_f{ii}-baseline{2-mod(ii,2)})/baseline{2-mod(ii,2)};
                            mSFC{ii}=mean(SFC_f{ii},1)*100;
                            ee{ii}=std(SFC_f{ii},0,1)/sqrt(sum(idx1))*100;
                        end
                        for ialign=1:2
                            if ialign==1
                                Time=Time1;
                            else
                                Time=Time2;
                            end
                            num=num+1;
                            %                         subplot(1,2,ialign);
                            subplot(2,2,num);
                            ss1=zeros(1,size(SFC_f{1+(ialign-1)*2},2));
                            for itime=1:size(SFC_f{1+(ialign-1)*2},2)
                                ss1(itime)=ttest2(SFC_f{1+(ialign-1)*2}(:,itime), SFC_f{2+(ialign-1)*2}(:,itime));
                            end                            
                            pp(1)=patchplot(Time, mSFC{1+(ialign-1)*2}, ee{1+(ialign-1)*2}, clr3(iarea,:)); hold on;
                            pp(2)=patchplot(Time, mSFC{2+(ialign-1)*2}, ee{2+(ialign-1)*2}, clr3(jarea,:)); hold on;
                            %                         pp(icelltype)=patchplot(Time, mSFC{2+(ialign-1)*2}, ee{2+(ialign-1)*2}, clr2(icelltype,:)); hold on;
                            ss2=ss2Baseline{1+(ialign-1)*2};
                            ss3=ss2Baseline{2+(ialign-1)*2};
                            ss1=ss1.*(ss2+ss3)>0;
                            plot(Time(ss1==1),ss1(ss1==1)*-5,'k*'); hold on;
                            plot(Time(ss2==1),ss2(ss2==1)*-10,'*','color', clr3(iarea,:)); hold on;
                            plot(Time(ss3==1),ss3(ss3==1)*-15,'*','color', clr3(jarea,:)); hold on;   %
%                             plot(Time(ss3==1),ss3(ss3==1)*-5*icelltype,'*','color', clr2(icelltype,:)); hold on;   %clr3(jarea,:)
                            plot([0 0],[-15 45],'--k'); hold on;
                            plot([-0.3 0.6],[0 0],'k'); hold on;
                            xlim([-0.3 0.5]);
                            ylim([-20 40]);
                            %                         ylim([-0.005 0.005]);
                            %                         if num<=3
                            %                             ll=legend(pp,{[AreaName{iarea} '->' AreaName{jarea}],[AreaName{jarea} '->' AreaName{iarea}]});
                            %                             set(ll,'box','off');
                            %                         end
                            yticks([-15 0 15 30 45]);
                            %                         yticklabels({'','',''});
                            xticks(-0.2:0.2:0.6);
                            %                         xticklabels({'','','',''});
                            set(gca,'box','off','linewidth',2.0,'fontsize',15);
                        end
                    end
                end
            end
        end
    end
end


%%
%%%%%% compare SFC in different areas
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
clr2=[1 0.4 0; 0 0 0];
% clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
clr3=[0 1 0;1 0 0; 0 0 1];
mmSFC=cell(1,2);
CellTypeName={'NarrowUnits','BroadUnits'};
igroup=4;
for icelltype=1:2
    figure;
    num=0;
    for isexo=0:1
        groupName=Group{igroup};
        C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %_HighAmp
        C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' CellTypeName{icelltype} '_SU.mat']);  %AttentionModulated
        C3=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_PairDate_' CellTypeName{icelltype} '_SU.mat']);
        C4=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_PairDate_' CellTypeName{icelltype} '_SU.mat']);
        Vars=fieldnames(C1.C);
        Time{1}=C1.C.t_cue; 
        Time{2}=C1.C.t_array;
        f=C1.C.f;
        data=cell(1,4);  SFC_f=cell(1,4);
        Dates=cell(1,4); 
        %             mSFC=cell(1,4); ee=cell(1,4);
        iarea=3;
        jarea=1;
        data{1}=cat(1, C1.C.CueOn{iarea,jarea}, C2.C.CueOn{iarea,jarea});
        data{2}=cat(1, C1.C.ArrayOn{iarea,jarea}, C2.C.ArrayOn{iarea,jarea});
       
        Dates{1}=[C3.PairDate.CueOn{iarea,jarea}+100, C4.PairDate.CueOn{iarea,jarea}];
        Dates{2}=[C3.PairDate.ArrayOn{iarea,jarea}+100, C4.PairDate.ArrayOn{iarea,jarea}];
        
        iarea=2;
        jarea=3;
        data{3}=cat(1, C1.C.CueOn{iarea,jarea}, C2.C.CueOn{iarea,jarea});
        data{4}=cat(1, C1.C.ArrayOn{iarea,jarea}, C2.C.ArrayOn{iarea,jarea});       
        Dates{3}=[C3.PairDate.CueOn{iarea,jarea}+100, C4.PairDate.CueOn{iarea,jarea}];
        Dates{4}=[C3.PairDate.ArrayOn{iarea,jarea}+100, C4.PairDate.ArrayOn{iarea,jarea}];
        ss2Baseline=cell(1,4);
        if ~isempty(data{1})
            idxF=f>=15&f<=30;
            for ii=1:2
                num=num+1;
                subplot(2,2,num);                
                idx1=~isnan(squeeze(sum(sum(data{ii},2),3)));
                SFC_f{ii}=squeeze(mean(data{ii}(idx1,:,idxF),3));
                Date_record1=Dates{ii}(idx1);
                sBaseline=mean(SFC_f{1}(:,Time{1}<-0.1&Time{1}>-0.25),2);
                baseline=mean(sBaseline,1);
                for itime=1:size(SFC_f{ii},2)
                    ss2Baseline{ii}(itime)=ttest2(SFC_f{ii}(:,itime), sBaseline);
                end
                normSFC=(SFC_f{ii}-baseline)/baseline*100;
                mSFC_delay=mean(normSFC(:,Time{ii}>=0.3&Time{ii}<=0.5),2);
%                 mSFC=mean(normSFC,1);
%                 ee=std(normSFC,0,1)/sqrt(sum(idx1));
%                 pp(1)=patchplot(Time{ii}, mSFC, ee, clr3(3,:)); hold on;
%                 ss=ss2Baseline{ii};
%                 plot(Time{ii}(ss==1),ss(ss==1)*-5,'*','color', clr3(3,:)); hold on;
%                 line([-0.3 0.5],[0 0],'color','k');
                
                idx1=~isnan(squeeze(sum(sum(data{2+ii},2),3)));
                SFC_f{ii+2}=squeeze(mean(data{2+ii}(idx1,:,idxF),3));
                Date_record2=Dates{ii+2}(idx1);
                sBaseline2=mean(SFC_f{3}(:,Time{1}<-0.1&Time{1}>-0.25),2);
                baseline2=mean(sBaseline2,1);
                for itime=1:size(SFC_f{ii+2},2)
                    ss2Baseline{ii+2}(itime)=ttest2(SFC_f{ii+2}(:,itime), sBaseline2);
                end
                normSFC=(SFC_f{ii+2}-baseline2)/baseline2*100;
                mSFC_delay2=mean(normSFC(:,Time{ii}>=0.3&Time{ii}<=0.5),2);
%                 mSFC2=mean(normSFC,1);
%                 ee2=std(normSFC,0,1)/sqrt(sum(idx1));
%                 pp(2)=patchplot(Time{ii}, mSFC2, ee2, clr3(2,:)); hold on;
%                 ss=ss2Baseline{ii+2};
%                 plot(Time{ii}(ss==1),ss(ss==1)*-10,'*','color', clr3(2,:)); hold on;
%                 
%                 xlim([-0.3 0.5]);
%                 xticks(-0.2:0.2:0.6);
%                 yticks([-15 0 15 30 45]);
%                 set(gca,'box','off','linewidth',2.0,'fontsize',15);
%                 legend(pp,{'FEF->LIP','PUL->FEF'});
%                 ylim([-20 30]);

                ValidDate=unique(Date_record1);
                Nday=numel(ValidDate);
                aveSFC1=zeros(1,Nday);
                aveSFC2=zeros(1,Nday);
                for jj=1:Nday
                    idxx=Date_record1==ValidDate(jj);
                    aveSFC1(jj)=mean(mSFC_delay(idxx));
                    idxx=Date_record2==ValidDate(jj);
                    aveSFC2(jj)=mean(mSFC_delay2(idxx));

                end
              plot(aveSFC1,aveSFC2, '*'); hold on;               
            end            
        end
    end
end

%% plot the spike_field_coherence between different areas during delay period for each individual monkey
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};

for igroup=4
    groupName=Group{igroup};
%     figure;
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
%             load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_cueResponsive_5sigma.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_Delay_cueResponsive_SU.mat']);
            load([folder 'f_spike-field-coherence_Delay.mat']);
%             load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_Baseline_AttentionModulated_SU.mat']);
%             load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_Delay_AttentionModulated_SU.mat']);
            Vars=fieldnames(C);
            ivar=2;
            figure;            
            for iarea=1:3
                for jarea=1:3  
                    if ~isempty(C.(Vars{ivar}){iarea,jarea})
                        idx1=~isnan(sum(C.(Vars{ivar}){iarea,jarea},2));
                        aa=mean(C.(Vars{ivar}){iarea,jarea}(idx1,:),1);
                        ee1=std(C.(Vars{ivar}){iarea,jarea}(idx1,:),0,1)/sqrt(sum(idx1));  
%                         idx0=~isnan(sum(C_base.CueOn{iarea,jarea},2));
%                         base1=mean(C_base.CueOn{iarea,jarea}(idx0,:),1);
%                         aa=aa-base1;
                        
                        if ivar<=2                            
                            idx2=~isnan(sum(C.(Vars{ivar+2}){iarea,jarea},2));
                            cc=mean(C.(Vars{ivar+2}){iarea,jarea}(idx2,:),1);
                            ee2=std(C.(Vars{ivar+2}){iarea,jarea}(idx2,:),0,1)/sqrt(sum(idx2));
%                             idx0=~isnan(sum(C_base.CueOff{iarea,jarea},2));
%                             base2=mean(C_base.CueOff{iarea,jarea}(idx0,:),1);
%                             cc=cc-base2;
                        end
                        
                        subplot(3,3,3*(iarea-1)+jarea);
                        patchplot(f{ivar},aa,ee1,[1,0,0]); hold on;
                        patchplot(f{ivar},cc,ee2,[0,0,1]);
                        xlabel('Frequency (Hz)');                        
                        ylabel('SFC');     
                        xlim([5 80])
                        title([AreaName{iarea} '-' AreaName{jarea}]);                        
                    end
                end
            end
        end
    end
end

%% converge the spike_field_coherence between different areas in delay period of two monkeys together
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};

for igroup=4
    groupName=Group{igroup};
    for isexo=0:1
%         C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_delay_AttentionModulated_SU.mat']);  
%         C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_delay_AttentionModulated_SU.mat']);
        C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_delay_cueResponsive_SU.mat']);  %
        C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_delay_cueResponsive_SU.mat']);
        load([folder 'f_spike-field-coherence_Delay.mat']);
        Vars=fieldnames(C1.C);
        ivar=1;
        figure;
        num=0;
        for iarea=1:2
            for jarea=(iarea+1):3
                data1=cat(1, C1.C.(Vars{ivar}){iarea,jarea}, C2.C.(Vars{ivar}){iarea,jarea});   
                data2=cat(1, C1.C.(Vars{ivar}){jarea,iarea}, C2.C.(Vars{ivar}){jarea,iarea});   
                if ~isempty(data1)
                    idx1=~isnan(sum(data1,2));
                    aa=mean(data1(idx1,:),1);     
                    aa=conv(aa,ones(1,5)/5,'same');
                    ee1=std(data1(idx1,:),0,1)/sqrt(sum(idx1));
                    
                    idx2=~isnan(sum(data2,2));
                    bb=mean(data2(idx2,:),1);
                    bb=conv(bb,ones(1,5)/5,'same');
                    ee2=std(data2(idx2,:),0,1)/sqrt(sum(idx2));
                    
%                     if ivar <= 2   %%%%% off condition
%                         data3=cat(1, C1.C.(Vars{ivar+2}){iarea,jarea}, C2.C.(Vars{ivar+2}){iarea,jarea});
%                         idx3=~isnan(sum(data3,2));                        
%                         cc=mean(data3(idx3,:),1);
%                         ee2=std(data3(idx3,:),0,1)/sqrt(sum(idx3));
%                         
%                         idxFreq=find(C1.C.f<=100);
%                         for itime=idxFreq
%                             [ss(itime),pp(itime)]=ttest2(data1(idx1,itime),data3(idx2,itime));
%                         end
%                     end                    
                        
                    num=num+1;
                    subplot(1,3,num);
                    patchplot(f{ivar},aa,ee1,[1 0 0]); hold on;
                    patchplot(f{ivar},bb,ee2,[0 0 1]); hold on;
%                     plot(C1.C.f(idxFreq),ss*0.06,'.k'); hold on;
                    xlabel('Frequency (Hz)');
                    ylabel('SFC');                    
                    title([AreaName{iarea} '-' AreaName{jarea}]);
                    xlim([8 80]);
%                     ylim([0.05 0.1]);
                end
            end
        end
    end
end

%% get the Granger causality between Spikes/LFPs from different areas
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
CellTypeName={'NarrowUnits','BroadUnits'};
iCellType=1;
mtm.Fs=1000; % sampling frequency
mtm.fpass=[0 150]; % band of frequencies to be kept
mtm.tapers=[3 5]; % taper parameters
% mtm.movingwin=[0.25 0.025];
mtm.trialave = 1;
mtm.pad=0;
mtm.err=[1 0.05];

t=cell(1,4);
T_cue=-600:1000;
T_array=-500:900;
SpkRange=[5 10; 12 25];

for igroup=1
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
           load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']);  %  _5sigma
           load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);           
           Variables=fieldnames(SpikeLFPInfo);
                                 
           CellType=cell(size(SpikeLFPInfo.ElecInfo,1),3);
           for iarea=1:3
               CountUnit=0;
               for isession=1:size(SpikeLFPInfo.ElecInfo,1)
                   if ~isempty(SpikeLFPInfo.ElecInfo{isession, iarea})
                       NN=size(SpikeLFPInfo.ElecInfo{isession, iarea},1);
                       CellType{isession,iarea}=SpikeShape.CellType{iarea}(CountUnit+1:CountUnit+NN);
                       CountUnit=CountUnit+NN;
                   end
               end
           end
    
           for ivar=2:5 % 2-5 are spikes, 6-9 are LFPs
               count=zeros(3);
               for iarea=1:2
                   for jarea=(iarea+1):3                       
                       for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                           if ~isempty(SpikeLFPInfo.(Variables{ivar}){idate,iarea})&& ~isempty(SpikeLFPInfo.(Variables{ivar}){idate,jarea})
                               for icell=1:numel(SpikeLFPInfo.(Variables{ivar}){idate,iarea})
                                   for jcell=1:numel(SpikeLFPInfo.(Variables{ivar}){idate,jarea})
                                       if SpikeLFPInfo.ElecInfo{idate,iarea}(icell,2)== SpikeLFPInfo.ElecInfo{idate,jarea}(jcell,2)                                           
                                           data1=transpose(SpikeLFPInfo.(Variables{ivar}){idate,iarea}{icell});
                                           data2=transpose(SpikeLFPInfo.(Variables{ivar}){idate,jarea}{jcell});
                                           ValidTrl=sum(data1,1)>=1 & sum(data2,1)>=1; %(>=5)                     
                                           
                                           %%%%%% for LFP
%                                             data1=transpose(SpikeLFPInfo.(Variables{ivar}){idate,iarea}{icell});
%                                             data2=transpose(SpikeLFPInfo.(Variables{ivar}){idate,jarea}{jcell});
%                                             ValidTrl = ~isnan(mean(data1, 1))&~isnan(mean(data2, 1))& mean(data1, 1)~=0 & mean(data2, 1)~=0;                                                                                         
                                           
                                            if sum(ValidTrl)>20 
                                                count(iarea,jarea)=count(iarea,jarea)+1;   
                                                Spike1=data1(:,ValidTrl);   %LFP1(itime:itime+250,:)
                                                Spike2=data2(:,ValidTrl);    %LFP2(itime:itime+250,:)                                               
                                                for itrl=1:size(Spike1,2)
                                                    Spike1(:, itrl)= conv(Spike1(:, itrl),ones(10,1)/10,'same');     % transform discrete signal to continuous signal                                                
                                                    Spike2(:, itrl)= conv(Spike2(:, itrl),ones(10,1)/10,'same');                                                    
                                                end
                                                
%                                                 LFP1=rmlinesc(data1(:,ValidTrl),mtm,[], 0, 60);
%                                                 LFP1=LFP1-repmat(mean(LFP1,1),size(LFP1,1),1);
%                                                 LFP2=rmlinesc(data2(:,ValidTrl),mtm,[], 0, 60);
%                                                 LFP2=LFP2-repmat(mean(LFP2,1),size(LFP2,1),1);
                                                num=0;
                                                for itime=1:25:size(Spike1,1)-250   %size(Spike1,1)-250    %
                                                    num=num+1;                                                                                                        
                                                    if sum(sum(Spike1(itime:itime+250,:),1))~=0&& sum(sum(Spike2(itime:itime+250,:),1))~=0
                                                        [GCx2y,GCy2x,f] = LFP_PairwiseGranger(Spike1(itime:itime+250,:),Spike2(itime:itime+250,:),mtm);  %[GCx2y,GCy2x,f] =spikes_PairwiseGranger(Spike1(itime:itime+250,:),Spike2(itime:itime+250,:),mtm);
                                                        GC.(Variables{ivar}(1:end-5)){iarea,jarea}(count(iarea,jarea),num,:)=GCx2y;
                                                        GC.(Variables{ivar}(1:end-5)){jarea,iarea}(count(iarea,jarea),num,:)=GCy2x;
                                                    else
                                                        GC.(Variables{ivar}(1:end-5)){iarea,jarea}(count(iarea,jarea),num,:)=zeros(1,length(f));
                                                        GC.(Variables{ivar}(1:end-5)){jarea,iarea}(count(iarea,jarea),num,:)=zeros(1,length(f));
                                                    end
                                                    
%                                                     [GCx2y,GCy2x,f] = LFP_PairwiseGranger(LFP1(itime:itime+250,:),LFP2(itime:itime+250,:),mtm); 
%                                                     GC.(Variables{ivar}(1:end-3)){iarea,jarea}(count(iarea,jarea),num,:)=GCx2y;
%                                                     GC.(Variables{ivar}(1:end-3)){jarea,iarea}(count(iarea,jarea),num,:)=GCy2x;
                                                end
                                                SpikeType.(Variables{ivar}(1:end-5)){iarea,jarea}(count(iarea,jarea),:)=[CellType{idate,iarea}(icell),CellType{idate,jarea}(jcell)];
                                            end
                                       end
                                   end
                               end
                           end                           
                       end
                   end
               end                
           end
           GC.t_cue=T_cue(1:25:end-250)+125;
           GC.t_array=T_array(1:25:end-250)+125;
           GC.f=f;
           save([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_Spikes_between_areas_SU.mat'],'GC','-v7.3');
           save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeType_GC_Spikes_SU.mat'],'SpikeType','-v7.3');
        end
    end
end
%% plot the Granger causality results for cued and uncued conditions between a pair of regions 
clear all;
% close all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};

figure;
clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
for igroup=4
    groupName=Group{igroup};
    for imonkey=2
        monkey=monkeyset{imonkey};
        for isexo=0:1             
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_Spikes_between_areas_SU.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeType_GC_Spikes_SU.mat']);
            Variables=fieldnames(GC);
            ivar=2;
            
            if ivar==2|| ivar==4
                Var1=ivar-1;
            else
                Var1=ivar;
            end
            if ivar==1 || ivar==3
                time=GC.t_cue;
            else
                time=GC.t_array;
            end
            subplot(1,2,isexo+1);
            iarea=1;
            jarea=2 ; %iarea+1  % jarea>iarea
            clr2=[clr3(iarea,:); 0 0 0];  %clr3(jarea,:)
%             idx0=SpikeType.(Variables{ivar}){iarea,jarea}(:,1)>=1&SpikeType.(Variables{ivar}){iarea,jarea}(:,2)>=1;
            idx1=squeeze(mean(mean(GC.(Variables{ivar}){iarea,jarea},2),3))>0;  %&SpikeType.(Var){iarea,jarea}(:,1)==1;
            idx2=squeeze(mean(mean(GC.(Variables{ivar+2}){iarea,jarea},2),3))>0;  %&SpikeType.(Var){iarea,jarea}(:,2)==1;
            idx_f=GC.f>=10&GC.f<=35;
            GC_f{1}=squeeze(mean(GC.(Variables{ivar}){iarea,jarea}(idx1,:,idx_f),3));
            GC_f{2}=squeeze(mean(GC.(Variables{ivar+2}){iarea,jarea}(idx2,:,idx_f),3));
            
            GC_BS_f{1}=squeeze(mean(GC.(Variables{Var1}){iarea,jarea}(idx1,:,idx_f),3));
            GC_BS_f{2}=squeeze(mean(GC.(Variables{Var1+2}){iarea,jarea}(idx2,:,idx_f),3));
            sBase=cell(1,2);
            ss=zeros(2,length(time));
            for idir=1:2
                sBase{idir}=mean(GC_BS_f{idir}(:,GC.t_cue>=-250&GC.t_cue<=-50),2);
                for itime=1:length(time)
                    ss(idir,itime)=ttest(GC_f{idir}(:,itime), sBase{idir});
                end
                Base=mean(sBase{idir},1);
                GC_f{idir}=(GC_f{idir}-Base);  %/Base;
                mGC_f=mean(GC_f{idir},1);
                ee_f=std(GC_f{idir},0,1)/sqrt(sum(idx1));
                pp(idir)=patchplot(time, mGC_f, ee_f, clr2(idir,:)); hold on;
                idx00=ss(idir,:)==1;
                %                         plot(time(idx00),ss(idir,idx00)*-0.5*idir, '.','color', clr2(idir,:)); hold on;
                xlim([-200 600]);
                %                         ylim([-2 3]);
            end
            line([-200 600],[0 0],'color','k');
            sig=zeros(1, length(time));
            for itime=1:length(time)
                sig(itime)=ttest2(GC_f{1}(:,itime), GC_f{2}(:,itime));
            end
            idx11=sig==1;
            plot(time(idx11),sig(idx11)*0, '.','color', 'k'); hold on;
                    
%                     subplot(2,2,1+2*isexo);
%                     imagesc(time, GC.f, mGC1');
%                     axis xy;
%                     xlim([-200 600]);
%                     ylim([5 100]);
%                     title([num2str(iarea) '->' num2str(jarea)]);
%                     
%                     subplot(2,2,2+2*isexo);
%                     imagesc(time, GC.f, mGC2');
%                     axis xy;
%                     xlim([-200 600]);
%                     ylim([5 100]);
%                     title([num2str(jarea) '->' num2str(iarea)]);
              
        end
        legend(pp, {'Cued','Uncued'});
    end
end

%% plot the Granger causality for each directions between a pair of regions 
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
close all;
figure;
clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
for igroup=4
    groupName=Group{igroup};
    for imonkey=2
        monkey=monkeyset{imonkey};
        for isexo=0:1             
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_Spikes_between_areas_SU.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeType_GC_Spikes_SU.mat']);
            Variables=fieldnames(GC);
            ivar=2;            
            if ivar==2|| ivar==4
                Var1=ivar-1;
            else
                Var1=ivar;
            end
            if ivar==1 || ivar==3
                time=GC.t_cue;
            else
                time=GC.t_array;
            end
            subplot(1,2,isexo+1);
            iarea=1;
            jarea=2; %iarea+1  % jarea>iarea
            clr2=[clr3(iarea,:); clr3(jarea,:)];  %clr3(jarea,:)
            %                     idx0=SpikeType.(Variables{ivar}){iarea,jarea}(:,1)>=1&SpikeType.(Variables{ivar}){iarea,jarea}(:,2)>=1;
            idx0=squeeze(mean(mean(GC.(Variables{ivar}){iarea,jarea},2),3))>0 & squeeze(mean(mean(GC.(Variables{ivar}){jarea,iarea},2),3))>0 ; 
            idx1=idx0&SpikeType.(Variables{ivar}){iarea,jarea}(:,1)==1;  %&SpikeType.(Variables{ivar}){iarea,jarea}(:,2)==2;
            idx2=idx0&SpikeType.(Variables{ivar}){iarea,jarea}(:,2)==1;  %&SpikeType.(Variables{ivar}){iarea,jarea}(:,1)==2;
            idx_f=GC.f>=10&GC.f<=35;
            GC_f{1}=squeeze(mean(GC.(Variables{ivar}){iarea,jarea}(idx1,:,idx_f),3));
            GC_f{2}=squeeze(mean(GC.(Variables{ivar}){jarea,iarea}(idx2,:,idx_f),3));
            
            GC_BS_f{1}=squeeze(mean(GC.(Variables{Var1}){iarea,jarea}(idx1,:,idx_f),3));
            GC_BS_f{2}=squeeze(mean(GC.(Variables{Var1}){jarea,iarea}(idx2,:,idx_f),3));
            sBase=cell(1,2);
            ss=zeros(2,length(time));
            for idir=1:2
                sBase{idir}=mean(GC_BS_f{idir}(:,GC.t_cue>=-250&GC.t_cue<=-100),2);
                for itime=1:length(time)
                    ss(idir,itime)=ttest(GC_f{idir}(:,itime), sBase{idir});
                end
                Base=mean(sBase{idir},1);
                GC_f{idir}=(GC_f{idir}-Base);  %/Base;
                mGC_f=mean(GC_f{idir},1);
                ee_f=std(GC_f{idir},0,1)/sqrt(size(GC_f{idir},1));
                pp(idir)=patchplot(time, mGC_f, ee_f, clr2(idir,:)); hold on;
                idx00=ss(idir,:)==1;
                %                         plot(time(idx00),ss(idir,idx00)*-0.5*idir, '*','color', clr2(idir,:)); hold on;
                xlim([-200 600]);
                %                         ylim([-2 3]);
            end
            line([-200 600],[0 0],'color','k');
            sig=zeros(1, length(time));
            for itime=1:length(time)
                sig(itime)=ttest2(GC_f{1}(:,itime), GC_f{2}(:,itime));
            end
            idx11=sig==1;
            plot(time(idx11),sig(idx11)*0, '*','color', 'k'); hold on;
                    
%                     subplot(2,2,1+2*isexo);
%                     imagesc(time, GC.f, mGC1');
%                     axis xy;
%                     xlim([-200 600]);
%                     ylim([5 100]);
%                     title([num2str(iarea) '->' num2str(jarea)]);
%                     
%                     subplot(2,2,2+2*isexo);
%                     imagesc(time, GC.f, mGC2');
%                     axis xy;
%                     xlim([-200 600]);
%                     ylim([5 100]);
%                     title([num2str(jarea) '->' num2str(iarea)]);
              
        end
%         legend(pp, {'Cued','Uncued'});
    end
end

%% plot the Granger causality between LFPs for each directions between a pair of regions 
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
close all;
figure;
clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];
for igroup=4
    groupName=Group{igroup};
    for imonkey=2
        monkey=monkeyset{imonkey};
        for isexo=0:1             
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_LFP_between_areas_SU.mat']);           
            Variables=fieldnames(GC);
            ivar=1;            
            if ivar==2|| ivar==4
                Var1=ivar-1;
            else
                Var1=ivar;
            end
            if ivar==1 || ivar==3
                time=GC.t_cue;
            else
                time=GC.t_array;
            end
%             subplot(1,2,isexo+1);
            iarea=1;
            jarea=3; %iarea+1  % jarea>iarea
            clr2=[clr3(iarea,:); clr3(jarea,:)];  %clr3(jarea,:)           
            idx0=squeeze(mean(mean(GC.(Variables{ivar}){iarea,jarea},2),3))>0 & squeeze(mean(mean(GC.(Variables{ivar}){jarea,iarea},2),3))>0 ; 
            idx_f=GC.f>=10&GC.f<=35;
            
            mGC1=transpose(squeeze(mean(GC.(Variables{ivar}){iarea,jarea}(idx0,:,:),1)));
            mGC2=transpose(squeeze(mean(GC.(Variables{ivar}){jarea,iarea}(idx0,:,:),1)));
            subplot(2,2,isexo*2+1);
            imagesc(time, GC.f, mGC1); 
            axis xy;
            subplot(2,2,isexo*2+2);
            imagesc(time, GC.f, mGC2); 
            axis xy;
            
            
            
            
            GC_f{1}=squeeze(mean(GC.(Variables{ivar}){iarea,jarea}(idx0,:,idx_f),3));
            GC_f{2}=squeeze(mean(GC.(Variables{ivar}){jarea,iarea}(idx0,:,idx_f),3));
            
%             GC_BS_f{1}=squeeze(mean(GC.(Variables{Var1}){iarea,jarea}(idx0,:,idx_f),3));
%             GC_BS_f{2}=squeeze(mean(GC.(Variables{Var1}){jarea,iarea}(idx0,:,idx_f),3));
%             sBase=cell(1,2);
%             ss=zeros(2,length(time));
%             for idir=1:2
%                 sBase{idir}=mean(GC_BS_f{idir}(:,GC.t_cue>=-250&GC.t_cue<=-100),2);
%                 for itime=1:length(time)
%                     ss(idir,itime)=ttest(GC_f{idir}(:,itime), sBase{idir});
%                 end
%                 Base=mean(sBase{idir},1);
%                 GC_f{idir}=(GC_f{idir}-Base);  %/Base;
%                 mGC_f=mean(GC_f{idir},1);
%                 ee_f=std(GC_f{idir},0,1)/sqrt(size(GC_f{idir},1));
%                 pp(idir)=patchplot(time, mGC_f, ee_f, clr2(idir,:)); hold on;
%                 idx00=ss(idir,:)==1;
%                 %                         plot(time(idx00),ss(idir,idx00)*-0.5*idir, '*','color', clr2(idir,:)); hold on;
%                 xlim([-200 600]);
%                 %                         ylim([-2 3]);
%             end
%             line([-200 600],[0 0],'color','k');
%             sig=zeros(1, length(time));
%             for itime=1:length(time)
%                 sig(itime)=ttest2(GC_f{1}(:,itime), GC_f{2}(:,itime));
%             end
%             idx11=sig==1;
%             plot(time(idx11),sig(idx11)*0, '*','color', 'k'); hold on;                                 
        end
    end
end


%% get the spike-spike granger causality between different area
close all;
clear all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
T_cue=-600:1000;
T_array=-500:900;
for igroup=4
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_5sigma.mat']);            
            clear GCpair;
            Vars=fieldnames(SpikeLFPInfo);           
            for ivar=2:5
                count=zeros(3);
                for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                    ivar, idate,
                    for iarea=1:2   
                        for jarea=(iarea+1):3
                            if ~isempty(SpikeLFPInfo.(Vars{ivar}){idate,iarea})&&~isempty(SpikeLFPInfo.(Vars{ivar}){idate,jarea})
                                for icell=1:numel(SpikeLFPInfo.(Vars{ivar}){idate,iarea})
%                                     LFP1=transpose(SpikeLFPInfo.(Vars{ivar}){idate,iarea}{icell});
%                                     LFP1=rmlinesc(LFP1,params,[], 0, 60);
%                                     LFP1=LFP1-repmat(mean(LFP1,1),size(LFP1,1),1);
                                    Spike1=transpose(SpikeLFPInfo.(Vars{ivar}){idate,iarea}{icell});                                    
                                    for jcell=1:numel(SpikeLFPInfo.(Vars{ivar}){idate,jarea})
%                                         if SpikeLFPInfo.ElecInfo{idate,iarea}(icell,2)==SpikeLFPInfo.ElecInfo{idate,jarea}(jcell,2)    %size(Spike1,2)==size(SpikeLFPInfo.(Vars{ivar-4}){idate,jarea}{jcell},1)
%                                             LFP2=transpose(SpikeLFPInfo.(Vars{ivar}){idate,jarea}{jcell});
%                                             LFP2=rmlinesc(LFP2,params,[], 0, 60);
%                                             LFP2=LFP2-repmat(mean(LFP2,1),size(LFP2,1),1);
                                            Spike2=transpose(SpikeLFPInfo.(Vars{ivar}){idate,jarea}{jcell});
                                            ValidTrl=sum(Spike1,1)>=1&sum(Spike2,1)>=1;
                                            if sum(ValidTrl)>=15 
                                                count(iarea,jarea)=count(iarea,jarea)+1;
                                                count(jarea,iarea)=count(iarea,jarea); % pair numbers are the same in two directions
                                                
                                                itimeIdx=0;
                                                for itime=1:20:size(Spike1,1)-200
                                                    itimeIdx=itimeIdx+1;                                                    
                                                    if mean(mean(Spike1(itime:itime+200, :),2),1)>0 && mean(mean(Spike2(itime:itime+200, :),2),1)>0
                                                        [GCx2y,GCy2x,f] = spike_PairwiseGranger(Spike1(itime:itime+200,:),Spike2(itime:itime+200,:),params);
                                                    else
                                                        GCx2y=zeros(1,51);
                                                        GCy2x=zeros(1,51);
                                                    end
                                                    
                                                    if length(GCx2y)==51
                                                        GCpair.(Vars{ivar}(1:end-5)){iarea, jarea}(count(iarea,jarea),itimeIdx,:)=GCx2y;
                                                        GCpair.(Vars{ivar}(1:end-5)){jarea, iarea}(count(jarea,iarea),itimeIdx,:)=GCy2x;                                                        
                                                    end
                                                end                                               
                                            end
%                                         end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            GCpair.f=f;
            GCpair.t_cue=(0:20:1400)-500; % Spike signal length = 1601
            GCpair.t_array=(0:20:1200)-400;  
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_between_areas_cueResponsive_5sigma.mat'],'GCpair','-v7.3');            
        end
    end
end

%% plot the granger causality between different areas for each individual monkey
close all;
% clear all;
AreaName={'LIP','PUL','FEF'};
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
t_cue=-0.4990:0.02:0.9010;
t_array=-0.3990:0.02:0.8010;
for igroup=2
    groupName=Group{igroup};
%     figure;
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_GC_between_areas_cueResponsive_5sigma.mat']);
            Vars=fieldnames(GCpair);
            VarsNew=Vars([1 3 :5]);
            ivar=1;
            figure;
            for iarea=1:3
                for jarea=1:3  
                    if ~isempty(GCpair.(VarsNew{ivar}){iarea,jarea})
                        idx1=~isnan(squeeze(sum(sum(GCpair.(VarsNew{ivar}){iarea,jarea},2),3)));
                        aa=transpose(squeeze(mean(GCpair.(VarsNew{ivar}){iarea,jarea}(idx1,:,:),1)));
                        if ivar<=2
                            idx1=~isnan(squeeze(sum(sum(GCpair.(VarsNew{ivar+2}){iarea,jarea},2),3)));
                            cc=transpose(squeeze(mean(GCpair.(VarsNew{ivar+2}){iarea,jarea}(idx1,:,:),1)));
                        end
                        
                        if ivar==1||ivar==3
                            Time=GCpair.t_cue;
                            baseline_aa=repmat(mean(aa(:, Time>=-300& Time<=-50),2),1,size(aa,2));
                            if ivar==1
                                baseline_cc=repmat(mean(cc(:, Time>=-300& Time<=-50),2),1,size(cc,2));
                            end
                        elseif ivar==2||ivar==4
                            Time=GCpair.t_array;
                            idx2=~isnan(squeeze(sum(sum(GCpair.(VarsNew{ivar-1}){iarea,jarea},2),3)));
                            bb=transpose(squeeze(mean(GCpair.(VarsNew{ivar-1}){iarea,jarea}(idx2,:,:),1)));
                            baseline_aa=repmat(mean(bb(:, Time>=-300& Time<=-50),2),1,size(aa,2));
                            if ivar==2
                                idx2=~isnan(squeeze(sum(sum(GCpair.(VarsNew{ivar+1}){iarea,jarea},2),3)));
                                bb=transpose(squeeze(mean(GCpair.(VarsNew{ivar+1}){iarea,jarea}(idx2,:,:),1)));
                                baseline_cc=repmat(mean(bb(:, Time>=-300& Time<=-50),2),1,size(cc,2));
                            end
                            
                        end
                        aa=(aa-baseline_aa);  %./baseline_aa*100;
                        cc=(cc-baseline_cc);  %./baseline_cc*100;
                        
                        subplot(3,3,3*(iarea-1)+jarea);
                        imagesc(Time(Time>=-300&Time<500), GCpair.f(GCpair.f<90), cc(GCpair.f<90,Time>=-300&Time<500)); hold on;  %-cc(:,Time>=-0.3)
                        if ivar==1||ivar==3
                            xlabel('Time to cue onset (s)');
                        else
                            xlabel('Time to array onset (s)');
                        end
                        ylabel('Frequency (Hz)');
                        line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
%                         caxis([-15 15]);
                        
                        axis xy;
                        title([AreaName{iarea} '-' AreaName{jarea}]);
                        colorbar;
                    end
                end
            end
        end
    end
end


%% get the spike-field coherence within each areas during the delay period
close all;
clear all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
T_cue=-600:1000;
T_array=-500:900;
for igroup=1:3
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference.mat']);
            count=zeros(4,3);
            clear C;
            clear Phase;
            clear nSpike;
            Vars=fieldnames(SpikeLFPInfo);           
            for iarea=1:3
                for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                    for icell=1:numel(SpikeLFPInfo.CueOnSpike{idate,iarea})
                        for ivar=6:9
                            data1=transpose(SpikeLFPInfo.(Vars{ivar}){idate,iarea}{icell});
                            data1=rmlinesc(data1,params,[], 0, 60);
                            data1=data1-repmat(mean(data1,1),size(data1,1),1);
                            
                            idxT0=T_cue>=-400&T_cue<=0;
                            if ivar==6||ivar==8
                                idxT1=T_cue>=200&T_cue<=600;
                            else
                                idxT1=T_array>=200&T_array<=600;
                            end                            
                            Spike1=transpose(SpikeLFPInfo.(Vars{ivar-4}){idate,iarea}{icell});
                            ValidTrl1=sum(Spike1(idxT1,:),1)>=0; %(>=3 spikes, >=20 trials) (>=0 spikes, >=30 trials)
                            
                            if sum(ValidTrl1)>=30 %&& sum(ValidTrl2)>=15 && sum(ValidTrl3)>=15 && sum(ValidTrl4)>=15                                
                                count(ivar-5,iarea)=count(ivar-5,iarea)+1;
%                                 nSpike.(Vars{ivar}(1:end-3)){iarea}{count(ivar-5,iarea)}=sum(Spike1(idxT1,ValidTrl1),1);
                                LFP1=data1(idxT1,ValidTrl1);
                                [C.(Vars{ivar}(1:end-3)){iarea}(count(ivar-5,iarea),:),Phase.(Vars{ivar}(1:end-3)){iarea}(count(ivar-5,iarea),:),~,~,~,f]=coherencycpb(LFP1,Spike1(idxT1,ValidTrl1),params);
                                
                                %                             for irepeat=1:50
                                %                                 idxNew1=randperm(size(LFP1,2));
                                %                                 idxNew2=randperm(size(LFP2,2));
                                %                                 idxNew3=randperm(size(LFP3,2));
                                %                                 idxNew4=randperm(size(LFP4,2));
                                %
                                %                                 [C.CueOn{iarea}(count(iarea),irepeat,:,:),~,~,~,~,~,~,~]=cohgramcpb(LFP1(:,idxNew1),Spike1(:,ValidTrl1),movingwin,params);
                                %                                 [C.CueOff{iarea}(count(iarea),irepeat,:,:),~,~,~,~,t1,f1,~]=cohgramcpb(LFP2(:,idxNew2),Spike2(:,ValidTrl2),movingwin,params);
                                %                                 [C.ArrayOn{iarea}(count(iarea),irepeat,:,:),~,~,~,~,~,~,~]=cohgramcpb(LFP3(:,idxNew3),Spike3(:,ValidTrl3),movingwin,params);
                                %                                 [C.ArrayOff{iarea}(count(iarea),irepeat,:,:),~,~,~,~,t2,f2,~]=cohgramcpb(LFP4(:,idxNew4),Spike4(:,ValidTrl4),movingwin,params);
                                %                             end
                            end
                        end
                    end
                end
            end
            C.f=f;
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_delayed_spike-field-coherence_demeaned_AllChannelReference.mat'],'C','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_delayed_spike-field-phase_demeaned_AllChannelReference.mat'],'Phase','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_baseline_AllChannelReference.mat'],'C','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_delayed_nSpike.mat'],'nSpike','-v7.3');
        end
    end
end

%% plot the population results of power spectrum for each individual monkey
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
AreaName={'LIP','PUL','FEF'};
for igroup=4
    groupName=Group{igroup};
    for isexo=0:1
        figure;        
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            Coh=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_delayed_spike-field-coherence_demeaned_AllChannelReference.mat']);  %_AllChannelReference
            Base=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_baseline_AllChannelReference.mat']);  %_AllChannelReference

            for iarea=1:3
                data1=Coh.C.CueOn{iarea};   
                data2=Coh.C.CueOff{iarea};                
                idx00=~isnan(mean(Base.C.CueOn{iarea},2)); 
                baseline1=mean(Base.C.CueOn{iarea}(idx00,:),1);
                idx00=~isnan(mean(Base.C.CueOff{iarea},2));
                baseline2=mean(Base.C.CueOff{iarea}(idx00,:),1);
                
                idx1=~isnan(mean(data1,2));            
                aa=mean(data1(idx1,:),1);  %-baseline1;
                ee1=std(data1(idx1,:),0,1)/sqrt(sum(idx1));

                idx2=~isnan(mean(data2,2));
                bb=mean(data2(idx2,:),1);  %-baseline2;
                ee2=std(data2(idx2,:),0,1)/sqrt(sum(idx2));                       

                subplot(2,3,(imonkey-1)*3+iarea);
                patchplot(Coh.C.f, aa, ee1, [1 0 0]); hold on;
                patchplot(Coh.C.f, bb, ee2, [0 0 1]); hold on;
%                 plot(Coh.C.f, aa, 'r'); hold on;
%                 plot(Coh.C.f, bb, 'b'); hold on;
                xlim([0 60]);
            end
        end
    end
end

%% plot the population results of power spectrum with all data from two monkeys converged
% close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
AreaName={'LIP','PUL','FEF'};
figure;
for igroup=4
    groupName=Group{igroup};
    OnOffDiff=cell(1,3);
    for isexo=0:1
        C1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_delayed_spike-field-coherence_demeaned_AllChannelReference.mat']); %
        C2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_delayed_spike-field-coherence_demeaned_AllChannelReference.mat']);
                
        Vars=fieldnames(C1.C);
        ivar=2;
        for iarea=1:3
            data1=[C1.C.(Vars{ivar}){iarea}; C2.C.(Vars{ivar}){iarea}];
            data2=[C1.C.(Vars{ivar+2}){iarea}; C2.C.(Vars{ivar+2}){iarea}];
            idx1=~isnan(mean(data1,2));
            aa=mean(data1(idx1,:),1);
            ee1=std(data1(idx1,:),0,1)/sqrt(sum(idx1));
            
            idx2=~isnan(mean(data2,2));
            bb=mean(data2(idx2,:),1);
            ee2=std(data2(idx2,:),0,1)/sqrt(sum(idx2));
            
            subplot(2,3,isexo*3+iarea);
            patchplot(C1.C.f, aa, ee1, [1 0 0]); hold on;
            patchplot(C1.C.f, bb, ee2, [0 0 1]); hold on;
            xlim([0 60]);
        end
    end
end


%% get the LFP segments around individual spikes
close all;
clear all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Fs=1000;
params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];
T_cue=-600:1000;
T_array=-500:900;
for igroup=4
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference.mat']);
            count=zeros(1,3);
            LFPsegment.CueOn=cell(1,3);
            LFPsegment.CueOff=cell(1,3);
            LFPsegment.ArrayOn=cell(1,3);
            LFPsegment.ArrayOff=cell(1,3);
            clear SpikeTime;
            for iarea=1:3
                for idate=1:size(SpikeLFPInfo.ElecInfo,1)
                    for icell=1:numel(SpikeLFPInfo.CueOnSpike{idate,iarea})                        
                        data1=transpose(SpikeLFPInfo.CueOnLFP{idate,iarea}{icell});
                        data2=transpose(SpikeLFPInfo.CueOffLFP{idate,iarea}{icell});
                        data3=transpose(SpikeLFPInfo.ArrayOnLFP{idate,iarea}{icell});
                        data4=transpose(SpikeLFPInfo.ArrayOffLFP{idate,iarea}{icell});
                        
                        data1=rmlinesc(data1,params,[], 0, 60);
                        data2=rmlinesc(data2,params,[], 0, 60);
                        data3=rmlinesc(data3,params,[], 0, 60);
                        data4=rmlinesc(data4,params,[], 0, 60);

                        
                        Spike1=transpose(SpikeLFPInfo.CueOnSpike{idate,iarea}{icell});
                        Spike2=transpose(SpikeLFPInfo.CueOffSpike{idate,iarea}{icell});
                        Spike3=transpose(SpikeLFPInfo.ArrayOnSpike{idate,iarea}{icell});
                        Spike4=transpose(SpikeLFPInfo.ArrayOffSpike{idate,iarea}{icell});
                        
%                         idxT0=T_cue>=-400&T_cue<=0;
                        idxT1=T_cue>=200&T_cue<=600;
                        idxT2=T_array>=200&T_array<=600;
                        
                        ValidTrl1=find(sum(Spike1(idxT1,:),1)>0);
                        ValidTrl2=find(sum(Spike2(idxT1,:),1)>0);
                        ValidTrl3=find(sum(Spike3(idxT2,:),1)>0);
                        ValidTrl4=find(sum(Spike4(idxT2,:),1)>0);
                        
                        if numel(ValidTrl1)>=15&& numel(ValidTrl2)>=15&& numel(ValidTrl3)>=15&& numel(ValidTrl4)>=15
                            count(iarea)=count(iarea)+1;
                            nSpk=0;
                            for itrl=1:numel(ValidTrl1)
                                idxSpk=find(Spike1(idxT1,ValidTrl1(itrl))>0);
                                for ispk=1:numel(idxSpk)
                                    nSpk=nSpk+1;
                                    idxLFP=(idxSpk(ispk)+800-250):(idxSpk(ispk)+800+250);
                                    LFPsegment.CueOn{iarea}{count(iarea)}(:,nSpk)=data1(idxLFP,ValidTrl1(itrl));
%                                     SpikeTime.CueOn{iarea}{count(iarea)}(nSpk)=idxSpk(ispk)+199;
                                end
                            end
                            
                            nSpk=0;
                            for itrl=1:numel(ValidTrl2)
                                idxSpk=find(Spike2(idxT1,ValidTrl2(itrl))>0);
                                for ispk=1:numel(idxSpk)
                                    nSpk=nSpk+1;
                                    idxLFP=(idxSpk(ispk)+800-250):(idxSpk(ispk)+800+250);
                                    LFPsegment.CueOff{iarea}{count(iarea)}(:,nSpk)=data2(idxLFP,ValidTrl2(itrl));
%                                     SpikeTime.CueOff{iarea}{count(iarea)}(nSpk)=idxSpk(ispk)+199;
                                end
                            end
                            
                            nSpk=0;
                            for itrl=1:numel(ValidTrl3)
                                idxSpk=find(Spike3(idxT2,ValidTrl3(itrl))>0);
                                for ispk=1:numel(idxSpk)
                                    nSpk=nSpk+1;
                                    idxLFP=(idxSpk(ispk)+700-250):(idxSpk(ispk)+700+250);
                                    LFPsegment.ArrayOn{iarea}{count(iarea)}(:,nSpk)=data3(idxLFP,ValidTrl3(itrl));
%                                     SpikeTime.ArrayOn{iarea}{count(iarea)}(nSpk)=idxSpk(ispk)+199;
                                end
                            end
                            
                            nSpk=0;
                            for itrl=1:numel(ValidTrl4)
                                idxSpk=find(Spike4(idxT2,ValidTrl4(itrl))>0);
                                for ispk=1:numel(idxSpk)
                                    nSpk=nSpk+1;
                                    idxLFP=(idxSpk(ispk)+700-250):(idxSpk(ispk)+700+250);
                                    LFPsegment.ArrayOff{iarea}{count(iarea)}(:,nSpk)=data4(idxLFP,ValidTrl4(itrl));
%                                     SpikeTime.ArrayOff{iarea}{count(iarea)}(nSpk)=idxSpk(ispk)+199;
                                end
                            end
                        end
                        
                    end
                end
            end  
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_LFPsegment_around_spike_AllChannelReference.mat'],'LFPsegment','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeTime.mat'],'SpikeTime','-v7.3');
        end
    end
end

%% calculate the spike-LFP phase coupling 
Fs=1000;
Time=-250:250;
for igroup=4
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=0:1
            clear PPL;
            clear power_hb;
            clear SpkNum;
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_LFPsegment_around_spike_AllChannelReference.mat']);
%             load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeTime.mat']);
            Variables=fieldnames(SpikeTime);         
            Freq_center=4:2:57;
            
            for ii=1:numel(Variables)
                for iarea=1:3
                    imonkey, ii, iarea,
                    for icell=1:numel(LFPsegment.(Variables{ii}){iarea})  
                        num=0;
                        SpkNum.(Variables{ii}){iarea}(icell)=size(LFPsegment.(Variables{ii}){iarea}{icell},2);
%                         for ifreq=4:2:57
%                             num=num+1;
%                             signal_filtered=bandpass(LFPsegment.(Variables{ii}){iarea}{icell}, [ifreq-2 ifreq+2], Fs);
%                             signal_hilbert=hilbert(signal_filtered);
%                             power_hb.(Variables{ii}){iarea}(icell,num)=mean(abs(signal_hilbert(Time==0,:)));                           
%                             PPL.(Variables{ii}){iarea}(icell,num)=abs(mean(exp(1i.*angle(signal_hilbert(Time==0,:)))));
%                         end
                    end
                end
            end
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpkNum.mat'],'SpkNum','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_power_hb_AllChannelReference.mat'],'power_hb','-v7.3');  %_AllChannelReference
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pairwise_phase_locking_AllChannelReference.mat'],'PPL','-v7.3');
        end
    end
end

figure;
subplot(1,2,1);
hist(SpkNum.ArrayOn{3});
subplot(1,2,2);
hist(SpkNum.ArrayOff{3});

%% plot the pairwaise phase locking between spike and LFP within the same area, results are sensitive to spike counts (not reliable)
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
Freq_center=4:2:57;
for igroup=4
    groupName=Group{igroup};    
    for isexo=0:1
        PPL1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_pairwise_phase_locking_AllChannelReference.mat']);
        PPL2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_pairwise_phase_locking_AllChannelReference.mat']);
%         load([folder monkey '_' cuetype{isexo+1} '_' groupName '_power_hb_AllChannelReference.mat']);
        NSpk1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_SpkNum.mat']);
        NSpk2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_SpkNum.mat']);
        figure;
        Variables=fieldnames(PPL1.PPL);
        ii=1;
        for iarea=1:3
            data1=[PPL1.PPL.(Variables{ii}){iarea};  PPL2.PPL.(Variables{ii}){iarea}];
            Nspike=[NSpk1.SpkNum.(Variables{ii}){iarea} NSpk1.SpkNum.(Variables{ii}){iarea}];
            idx=Nspike>120;
            mm=mean(data1(idx,:),1);
            ee1=std(data1(idx,:),0,1)/sqrt(sum(idx));
            
            data2=[PPL1.PPL.(Variables{ii+1}){iarea};  PPL2.PPL.(Variables{ii+1}){iarea}];
            Nspike=[NSpk1.SpkNum.(Variables{ii+1}){iarea} NSpk1.SpkNum.(Variables{ii+1}){iarea}];
            idx=Nspike>120;
            nn=mean(data2(idx,:),1);
            ee2=std(data2(idx,:),0,1)/sqrt(sum(idx));
            
            %                 mm=mean(power_hb.(Variables{ii}){iarea},1);
            %                 ee1=std(power_hb.(Variables{ii}){iarea},0,1)/sqrt(size(power_hb.(Variables{ii}){iarea},1));
            %                 nn=mean(power_hb.(Variables{ii+1}){iarea},1);
            %                 ee2=std(power_hb.(Variables{ii+1}){iarea},0,1)/sqrt(size(power_hb.(Variables{ii+1}){iarea},1));
            subplot(1,3,iarea);
            patchplot(Freq_center, mm, ee1, 'r'); hold on;
            patchplot(Freq_center, nn, ee2, 'b'); hold on;
        end        
    end
end

