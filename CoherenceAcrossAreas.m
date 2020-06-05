%%%% coherence between different areas
close all;
clear all;
% dateUsed={'051018','042418','050718','051218','050118','050518','042118','051718','041718'};  %,'050418'
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};
% dateUsed={'111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey  '111717', '111917',
dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718', '092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
folder='Z:\RujiaChen\Results\';
AreaName={'LIP','PUL','FEF'};
params.tapers=[2,3];
params.pad=0;
params.Fs=1000;
params.fpass=[0 150];
params.err=[1 0.05];
params.trialave=1;
movingwin=[0.2, 0.02];

%% Get the Spikes and LFPs in different conditions for each area
TargetPosY=[200,-200];
cueLoc=[2,5];  % cue-on and off
RFtype=[1,2]; % different target shape

TimeUsed=[-250, 600];   %[-250 50];  %
monkey='Vasco';
% load([folder 'Facilitate_units.mat']);
% load([folder 'CuePosition.mat']);
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
load([folder 'VisualUnit.mat']);
load([folder 'CuePos.mat']);
load([folder 'RFLoc_' monkey '.mat']);
CuePosition=cell(numel(dateUsed), 3);
SelectChannel=cell(numel(dateUsed), 3);
for idate=1:numel(dateUsed)  %[1: 5  7: 9]   %
    date=dateUsed{idate};
%     load([folder 'CueOn_' date '.mat']);
    load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
    load([folder 'CorrectTrialParam_' date '.mat']);
    load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);   
    
    load([folder 'bArrayResp_' date '.mat']);
    load([folder 'CueOn_' date '.mat']);
    count=zeros(1,3);

    for iarea=1:3
%         load([folder 'EstimateRange_' AreaName{iarea} '.mat']);
%         isession=find(strcmp(RawInfo(1,:), date)==1);
%         channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96
            channel=1:96;
%         channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}(VisualUnit{idate, iarea}>0)));  % for Mikey  &SpikeTrain.SU==1
        
%         if ~isempty(channelID)
            for ich=(iarea-1)*32+2:iarea*32-1   % channelID  ##for Mikey
                nCond=0;
                if sum(bArrayResp(ich,:))>0    %sum(CueOn(ich,:))~=0  %&&EstimateRange(idate, ich-(iarea-1)*32)>0  %Facilitate(idate, ipos, ich)>0   %  
                    count(iarea)=count(iarea)+1;
                    SelectChannel{idate,iarea}(count(iarea))=SpikeTrain.channel(ich);
                    CuePosition{idate, iarea}(count(iarea),:)=zeros(1,4);
                    for ipos=1:2
                        for icue=1:2
                            nCond=nCond+1;
                            if CueOn(ich,ipos)==icue
                                CuePosition{idate, iarea}(count(iarea),nCond)=1;
                            elseif CueOn(ich,ipos)==0 && RFLoc{iarea}(iarea,ipos)==icue
                                CuePosition{idate, iarea}(count(iarea),nCond)=1;                                
                            end                    

                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==1 &CorrectTrialParam.is_mouse_trial==0;
                            if sum(idx1)>20
                                idxCh=find(channel==ich);  %SpikeTrain.channel(ich)
                                data1=transpose(squeeze(LFPinfo.cueAlign(idx1,idxCh,:)));   %%%raw LFP                                    
                                if ismember(channel(idxCh), [1 33 65])
                                    data2=transpose(squeeze(LFPinfo.cueAlign(idx1,idxCh,:))-squeeze(LFPinfo.cueAlign(idx1,idxCh+1,:)));
                                else
                                    data2=transpose(squeeze(LFPinfo.cueAlign(idx1,idxCh,:))-squeeze(LFPinfo.cueAlign(idx1,idxCh-1,:)));         % LFP corrected by neighbouring channel
                                end
%                                 data3=transpose(squeeze(LFPinfo.cueAlign(idx1,idxCh+1,:))+squeeze(LFPinfo.cueAlign(idx1,idxCh-1,:))-2*squeeze(LFPinfo.cueAlign(idx1,idxCh,:)));   
                                data0=transpose(squeeze(SpikeTrain.cueAlign(idx1,ich,:)));
                                                               
                                data=rmlinesc(data2,params,[], 0, 60);
%                                 CueCSD{idate,iarea}{nCond,count(iarea)}=data;  
                                CueLFP{idate,iarea}{nCond,count(iarea)}=data;  %(idxT1,:);
                                CueSpk{idate,iarea}{nCond,count(iarea)}=data0;  %(idxT2,:);


                                data1=transpose(squeeze(LFPinfo.ArrayAlign(idx1,idxCh,:)));  %%%raw LFP
                                if ismember(channel(idxCh), [1 33 65])
                                     data2=transpose(squeeze(LFPinfo.ArrayAlign(idx1,idxCh,:))-squeeze(LFPinfo.ArrayAlign(idx1,idxCh+1,:)));
                                else
                                    data2=transpose(squeeze(LFPinfo.ArrayAlign(idx1,idxCh,:))- squeeze(LFPinfo.ArrayAlign(idx1,idxCh-1,:))); % LFP corrected by neighbouring channal
                                end
%                                 data3=transpose(squeeze(LFPinfo.ArrayAlign(idx1,idxCh+1,:))+squeeze(LFPinfo.ArrayAlign(idx1,idxCh-1,:))-2*squeeze(LFPinfo.ArrayAlign(idx1,idxCh,:)));   
                                data0=transpose(squeeze(SpikeTrain.arrayAlign(idx1,ich,:)));

                                data=rmlinesc(data2,params,[], 0, 60);
%                                 ArrayCSD{idate,iarea}{nCond,count(iarea)}=data;                      
                                ArrayLFP{idate,iarea}{nCond,count(iarea)}=data;   %(idxT1,:);                               
                                ArraySpk{idate,iarea}{nCond,count(iarea)}=data0(6:end-5,:);  %idxT2
                            end
                        end
                    end
                end
            end
            fprintf('%d\n', idate);
%         end
    end
end

%%
cueTy='exo';
monkey='Vasco';  %'Mikey';
% Unit='VisualResp_unSorted';  % for spike
Unit='ArrayResp_unSorted';  % for spike
Unit1='ArrayResp_corrected_by_N-1';  % for LFP
save([folder  monkey  'CueLFP_' cueTy '_' Unit1 '.mat'], 'CueLFP', '-v7.3');
save([folder  monkey 'CueSpk_' cueTy '_' Unit '.mat'], 'CueSpk', '-v7.3');
save([folder  monkey 'ArrayLFP_' cueTy '_' Unit1 '.mat'], 'ArrayLFP', '-v7.3');
save([folder  monkey  'ArraySpk_' cueTy '_' Unit '.mat'], 'ArraySpk', '-v7.3');
% save([folder  monkey 'CuePostion_' cueTy '_CueOn.mat'], 'CuePostion', '-v7.3');
save([folder 'CuePosition_' monkey '.mat'], 'CuePosition', '-v7.3');  
%% compare the neural responses (spikes) in cue-on and cue-away conditions 
close all;
clear all;
folder='Z:\RujiaChen\Results\';
cueTy='exo';
monkey='Vasco';  %'Mikey';
Unit='ArrayResp_unSorted';  % for spike  VisualResp
load([folder  monkey 'CueSpk_' cueTy '_' Unit '.mat']);
load([folder  monkey 'ArraySpk_' cueTy '_' Unit '.mat']);
TimeCue=-600:1200;
TimeArray=-900:900;
TimeCoh=[-300 800];
%%
for idate=1:size(CueSpk,1)
    figure;
    for iarea=1:3
        for ipos=1:2
            subplot(3,4,(iarea-1)*4+ipos*2-1);
            RespDiff=zeros(size(CueSpk{idate,iarea},2), diff(TimeCoh)+1);
            for icell=1:size(CueSpk{idate,iarea},2)
                idx=TimeCue>=TimeCoh(1)&TimeCue<=TimeCoh(2);
                mm=mean(CueSpk{idate,iarea}{ipos*2-1,icell}(idx,:),2);
                nn=mean(CueSpk{idate,iarea}{ipos*2,icell}(idx,:),2);
                RespDiff(icell,:)=conv(nn-mm, ones(1,30)/30, 'same');
%                 plot(TimeCue(idx), RespDiff(icell,:)); hold on;
%                 xlim([-300 600])
            end
            plot(TimeCue(idx), mean(RespDiff,1)); hold on;
            line(TimeCoh, [0 0], 'color', [0 0 0]); hold on;
            xlim([-300 600]);
      
            subplot(3,4,(iarea-1)*4+ipos*2);
            RespDiff2=zeros(size(ArraySpk{idate,iarea},2), diff(TimeCoh)+1);
            for icell=1:size(ArraySpk{idate,iarea},2)
                idx=TimeArray>=TimeCoh(1)&TimeArray<=TimeCoh(2);
                mm=mean(ArraySpk{idate,iarea}{ipos*2-1,icell}(idx,:),2);
                nn=mean(ArraySpk{idate,iarea}{ipos*2,icell}(idx,:),2);
                RespDiff2(icell,:)=conv(nn-mm, ones(1,50)/50, 'same');
%                 plot(TimeArray(idx), RespDiff(icell,:)); hold on;
%                 xlim([-300 600]);
            end
            plot(TimeArray(idx), mean(RespDiff2,1)); hold on;
            line(TimeCoh, [0 0], 'color', [0 0 0]); hold on;
            xlim([-300 600])
        end
    end    
end

%% determine the RF location based on Cue response
RFLoc=cell(1,size(CueSpk,1));
for idate= 1:size(CueSpk,1)
    RFLoc{idate}=zeros(3,2);
    for iarea=1:3
        if ~isempty(CueSpk{idate,iarea})
            for ipos=1:2
                CueLoc=zeros(size(CueSpk{idate,iarea},2),1);
                for icell=1:size(CueSpk{idate,iarea},2)
                    idx=TimeCue>=50&TimeCue<=100;
                    mm=mean(CueSpk{idate,iarea}{ipos*2-1,icell}(idx,:),1);
                    nn=mean(CueSpk{idate,iarea}{ipos*2,icell}(idx,:),1);
                    if ttest2(mm,nn,'tail','right')==1
                        CueLoc(icell)=1;
                    elseif ttest2(mm,nn, 'tail','left')==1
                        CueLoc(icell)=2;
                    end
                end
                if sum(CueLoc)>0
                    if sum(CueLoc==1)>sum(CueLoc==2)
                        RFLoc{idate}(iarea,ipos)=1;
                    else
                        RFLoc{idate}(iarea,ipos)=2;
                    end
                end
            end
        end
    end
end
save([folder 'RFLoc_' monkey '.mat'], 'RFLoc', '-v7.3') ;

%% pooling array response based on the RF location
ArrayOn=cell(1,3);
ArrayOff=cell(1,3);
count=zeros(3,1);
for iarea=1:3
    for ipos=1:2
        for idate=1:numel(RFLoc)            
            if RFLoc{idate}(iarea,ipos)>0
                for icell=1:size(CueSpk{idate,iarea},2)                 %ArraySpk
                    mm=mean(CueSpk{idate,iarea}{ipos*2-2+RFLoc{idate}(iarea,ipos),icell},2);  % cue on condition
                    nn=mean(CueSpk{idate,iarea}{ipos*2+1-RFLoc{idate}(iarea,ipos),icell},2);
                    count(iarea)=count(iarea)+1;
                    ArrayOn{iarea}(count(iarea),:)=mm;
                    ArrayOff{iarea}(count(iarea),:)=nn;
                end
            end
        end
    end
end

%%
close all;
for iarea=1:3
    figure;
    mm=conv(mean(ArrayOn{iarea},1), ones(1,20)/20, 'same');
    nn= conv(mean(ArrayOff{iarea},1), ones(1,20)/20, 'same');
    plot(TimeArray, mm, 'r'); hold on;
    plot(TimeArray, nn, 'b'); hold on;
    xlim([-600 600])
end

%%  get the Coherence or Granger Causality between LFP-LFP/Spk_Spk of different areas
folder='Z:\RujiaChen\Results\';
cueTy='endo';   %'endo';
% monkey='Mikey';
monkey='Vasco';
% Unit='VisualResp_unSorted';  % for spike
Unit='ArrayResp_corrected_by_N-1';  % for LFP

% load([folder 'CuePostion_' cueTy '_CueOn.mat']);
% load([folder monkey 'ArrayLFP_' cueTy '_' Unit '.mat']);
% load([folder  monkey  'ArraySpk_' cueTy '_' Unit '.mat']);
load([folder  monkey 'CueLFP_' cueTy '_' Unit '.mat']);
% load([folder  monkey 'CueSpk_' cueTy '_' Unit '.mat']);

num=0; 
TimeCue=-600:1200;
TimeArray=-900:900;
TimeCoh=[-300 600];
source=CueLFP;  
% source=ArrayLFP;  
% source=ArraySpk;
% source=CueSpk; 
cohLFP=cell(1,4);
PhaseLFP=cell(1,4);
CellID=[];
for iarea=1:2
    for jarea=(iarea+1):3
        num=0;
        for idate=1:size(source,1)
            for icell=1:size(source{idate,iarea}, 2)
                for jcell=1:size(source{idate,jarea}, 2)
                    num=num+1;
                    CellID(num,:)=[icell, jcell];
                    for ii=1:size(source{idate,iarea}, 1)
                        if size(source{idate,iarea}{ii, icell}, 2)>=30 && size(source{idate,jarea}{ii, jcell}, 2)>=30
%                             idx=TimeArray>=TimeCoh(1)&TimeArray<=TimeCoh(2);
                            idx=TimeCue>=TimeCoh(1)&TimeCue<=TimeCoh(2);                                                     
                            
                            data1=source{idate,iarea}{ii, icell}(idx,:);
                            data2=source{idate,jarea}{ii, jcell}(idx,:);
                            %%%%% remove outlier trials
                            meanResp=mean(data1(300:600,:), 1)+mean(data2(300:600,:),1);
                            mu=mean(meanResp);
                            ee=std(meanResp);
                            idxValid=abs(meanResp-mu)<3*ee;
                            
                            [C,phi,S12,S1,S2,t,f]=cohgramc(data1(:, idxValid),data2(:, idxValid),movingwin,params);
%                             [C,phi,S12,S1,S2,t, f]=cohgrampb(data1,data2,movingwin,params);  %coherencyc(data1,data2,params);
                            cohLFP{ii}{iarea, jarea}(num,:,:)=C;
                            PhaseLFP{ii}{iarea, jarea}(num,:,:)=phi;

%                             itimeIdx=0;
%                             for itime=1:20:size(data1,1)-200
%                                 itimeIdx=itimeIdx+1;
%                                 if mean(mean(data1(itime:itime+200, :),2),1)>0&&mean(mean(data2(itime:itime+200, :),2),1)>0
%                                     [GCx2y,GCy2x,f] = spike_PairwiseGranger(data1(itime:itime+200, :),data2(itime:itime+200, :),params);
%                                 else
%                                     GCx2y=zeros(size(GCx2y));
%                                     GCy2x=zeros(size(GCy2x));
%                                 end
%                                 GCpair{ii}{iarea, jarea}(num,itimeIdx,:)=GCx2y;
%                                 GCpair{ii}{jarea, iarea}(num,itimeIdx,:)=GCy2x;
%                             end
                        end
                    end
                end
            end
        end     
    end
end

%%
folder='Z:\RujiaChen\Results\';
dataSource='cue_LFP';   %'array_LFP'; % 'SU'
% save([folder 'GC_cue_' dataSource '_' cueTy '_' monkey '_' Unit '.mat'] , 'GCpair' ,'-v7.3');
save([folder 'Coh_' dataSource '_' cueTy '_' monkey '_' Unit '.mat'] , 'cohLFP' ,'-v7.3');
%% plot the GC and GC difference between different areas
cueTy='exo';
load([folder 'GC_cue_SU_' cueTy '_' monkey '_' Unit '.mat']);
TimeCenter=[TimeCoh(1):20:TimeCoh(2)-200]+100;
close all;
Area={'LIP', 'PUL','FEF'};
CndName={'On', 'Away'};
for iarea=1:2
    for jarea=(iarea+1):3
        figure;
        for ii=1:2  %numel(GCpair)
            subplot(2,2,ii);
            idx1=squeeze(mean(mean(GCpair{ii}{iarea,jarea},2),3))~=0;
            mm=squeeze(mean(GCpair{ii}{iarea,jarea}(idx1,:,:),1));
            imagesc(TimeCenter, f, mm');
            axis xy;
            cc=colorbar;
            title([Area{iarea} '->' Area{jarea} ':' CndName{ii}]);
            
            subplot(2,2,2+ii);
            idx1=squeeze(mean(mean(GCpair{ii}{jarea,iarea},2),3))~=0;
            mm=squeeze(mean(GCpair{ii}{jarea,iarea}(idx1,:,:),1));
            imagesc(TimeCenter, f, mm');
            axis xy;
            caxis(cc.Limits);
            colorbar;    
            title([Area{jarea} '->' Area{iarea} ':'  CndName{ii}]);
        end
        
%         for ii=1:2   %numel(GCpair)
%             subplot(2,2,ii);
%             idx1=squeeze(mean(mean(GCpair{ii}{iarea,jarea},2),3))~=0;
%             mm=squeeze(mean(GCpair{ii*2-1}{iarea,jarea}(idx1,:,:),1));
%             nn=squeeze(mean(GCpair{ii*2}{iarea,jarea}(idx1,:,:),1));
%             hh=imagesc(TimeCenter, f, mm'-nn');
%             axis xy;
%             
%             caxis([-0.006 0.006]);
%             colorbar;
%             title([Area{iarea} '->' Area{jarea}]);
%             
%             subplot(2,2,2+ii);
%             idx1=squeeze(mean(mean(GCpair{ii*2-1}{jarea,iarea},2),3))~=0;
%             mm=squeeze(mean(GCpair{ii*2-1}{jarea,iarea}(idx1,:,:),1));
%             nn=squeeze(mean(GCpair{ii*2}{jarea,iarea}(idx1,:,:),1));
%             imagesc(TimeCenter, f, mm'-nn');
%             axis xy;
%             caxis([-0.006 0.006]);
%             colorbar;     
%             
%             title([Area{jarea} '->' Area{iarea}]);
%         end
    end
end

%% plot the coherence and coherence difference between different areas
close all;
cueTy='exo';  %'endo';   %
TimeWindow='cue';  % 'array';  %
Unit='ArrayResp_corrected_by_N-1';  % for LFP
load([folder 'Coh_' TimeWindow  '_LFP_' cueTy '_' monkey '_' Unit '.mat']);

figure;
ifig=0;
time=t-0.3;  
Source=cohLFP;
for iarea=1:2
    for jarea=(iarea+1) :3
%         figure;
%          for icnd=1:4
%              data=squeeze(mean(cohLFP{icnd}{iarea, jarea},2));
%              yy=zeros(size(data,1),1);
% %              subplot(2,2,icnd);
%              for ii=1:size(data,1)  
%                  specPower(ii,:)=abs(fft(data(ii,f>40)));
%                  yy(ii)=sum(specPower(ii,[7 23]),2);  % 13 17 
% %                  if yy(ii)<1.5
% %                  plot(specPower(ii,:)); hold on;
% %                  end               
%              end             
% %              hist(yy); hold on
%             subplot(2,2,icnd);
% %             data=squeeze(mean(cohLFP{icnd}{iarea, jarea}(:,3,:),3));
% %             mVal=mean(data);
% %             ee=std(data);
% %             idxx1=abs(data-mVal)<3*ee;   
%             idxx1=yy<1.5;
%             idx1=~isnan(squeeze(mean(mean(cohLFP{icnd}{iarea, jarea},2),3)));            
%             mm=squeeze(mean(cohLFP{icnd}{iarea, jarea}(idxx1&idx1,:,:),1));            
%             imagesc(time, f, log(mm)');
%             axis xy;
% %             caxis([-0.015 0.015])
%             colorbar;
%          end
                  
        for icnd=1:2
            ifig=ifig+1;
            subplot(3,2,ifig);
%             data=squeeze(mean(cohLFP{icnd*2-1}{iarea, jarea}(:,3,:),3));
%             mVal=mean(data);
%             ee=std(data);
%             idxx1=abs(data-mVal)<3*ee;
%             data=squeeze(mean(cohLFP{icnd*2}{iarea, jarea}(:,3,:),3));
%             mVal=mean(data);
%             ee=std(data);
%             idxx2=abs(data-mVal)<3*ee;
            
            idx1=~isnan(squeeze(mean(mean(Source{icnd*2-1}{iarea, jarea},2),3)));
            idx2=~isnan(squeeze(mean(mean(Source{icnd*2}{iarea, jarea},2),3)));
            mm=squeeze(mean(Source{icnd*2-1}{iarea, jarea}(idx1,:,:),1));
            nn=squeeze(mean(Source{icnd*2}{iarea, jarea}(idx2,:,:),1));
%             imagesc(time, f, log(nn)'-log(mm)');
            imagesc(time, f,nn'-mm');
            axis xy;
            ylim([0 100]);
%             caxis([-0.15 0.35]);   % for cue response
            caxis([-0.02 0.06]);            
            colorbar;
        end
    end
end
%% get the outlier sessions

%%%%% coherence within certain window 
% Time=-250:600;  
% TimeCoh=[-250 50];  %[300 600] ; 
% for icell=1:size(source{idate,iarea}, 2)
%     for jcell=1:size(source{idate,jarea}, 2)
%         idxCue=find(CuePostion{idate,iarea}(icell,:)==1&CuePostion{idate,iarea}(icell,:)==CuePostion{idate,jarea}(jcell,:));
%         if ~isempty(idxCue)
%             num=num+1;
%             CuePos{num}=idxCue;          
%             CellID(num,:)=[icell, jcell];
%             for ii=1:size(source{idate,iarea}, 1)
%                 idx=Time>=TimeCoh(1)&Time<=TimeCoh(2);
%                 data1=source{idate,iarea}{ii, icell}(idx,:);
%                 data2=source{idate,jarea}{ii, jcell}(idx,:);
%     %             [C,phi,S12,S1,S2,t,f]=cohgramc(data1,data2,movingwin,params);
%                  [C,phi,S12,S1,S2,f]=coherencyc(data1,data2,params);
%                 cohLFP{ii}(num,:)=C;
%                 PhaseLFP{ii}(num,:)=phi;
%             end
%         end
%     end
% end
% Coherence_baseline=cohLFP;

%%  get the coherence between Spk_LFP of different areas
cueTy='exo';
load([folder 'CuePostion_' cueTy '_CueOn.mat']);
load([folder 'ArrayLFP_' cueTy '_CueOn.mat']);
load([folder 'ArraySpk_' cueTy '_CueOn.mat']);

% load([folder 'CueLFP_exo_CueOn.mat']);
% load([folder 'CueSpk_exo_CueOn.mat']);

source1=ArrayLFP;  
source2=ArraySpk;
Time=-250:600;  
TimeCoh=[300 600] ; %  %

idate=1;
for iarea=1:3
    for jarea=1:3
        num=0;
        for icell=1:size(source1{idate,iarea}, 2)
            for jcell=1:size(source2{idate,jarea}, 2)
                idxCue=find(CuePostion{idate,iarea}(icell,:)==1&CuePostion{idate,iarea}(icell,:)==CuePostion{idate,jarea}(jcell,:));
                if ~isempty(idxCue)
                    num=num+1;
                    CuePos{iarea,jarea}{num}=idxCue;
                    CellID{iarea,jarea}(num,:)=[icell, jcell];
                    for ii=1:size(source1{idate,iarea}, 1)
                        idx=Time>=TimeCoh(1)&Time<=TimeCoh(2);
                        data1=source1{idate,iarea}{ii, icell}(idx,:);
                        data2=source2{idate,jarea}{ii, jcell}(idx,:);
                        %             [C,phi,S12,S1,S2,t,f]=cohgramc(data1,data2,movingwin,params);
                        [C,phi,S12,S1,S2,f]=coherencycpb(data1,data2,params);
                        cohLFP{iarea,jarea}{ii}(num,:)=C;
                        PhaseLFP{iarea,jarea}{ii}(num,:)=phi;
                    end
                end
            end
        end
    end
end
% Coherence_baseline=cohLFP;
%% plot the coherence for each pair
figure;
clr4=colormap(hsv(4));
for ipair=1:size(cohLFP{1},1)
    subplot(9, 9, ipair);
    for jj=1:numel(cohLFP)
        idxf=f<100;
        plot(f(idxf), cohLFP{jj}(ipair, idxf), 'color', clr4(jj,:)); hold on;
    end
%     title(['[' num2str(CellID(ipair,1)) ',' num2str(CellID(ipair,1)) ']']);
    title(num2str(CuePos{ipair}));
end
legend('1', '2', '3', '4');

%% plot the population coherence  averaged across all unit pairs with close RFs
close all;
icue=4;
idxx=cellfun(@(x) ~isempty(find(x==icue, 1))+1,  CuePos, 'UniformOutput', false);
figure;
clr2=[1 0 0 ; 0 0 1];
ValidPair=cell2mat(idxx)==2;  %&transpose(ismember(CellID(:,1), [10, 11, 12]) &ismember(CellID(:,2), [7, 12, 13]) );
sig=zeros(1, length(f));
for ifreq=1:numel(f)
    data1=cohLFP{icue-1};  %-Coherence_baseline{icue-1};
    data2=cohLFP{icue};  %-Coherence_baseline{icue};
    sig(ifreq)=ttest(data1(ValidPair,ifreq), data2(ValidPair,ifreq)); 
end

for ii=1:2
    data=cohLFP{icue-2+ii}(ValidPair,:);  %-Coherence_baseline{icue-2+ii}(ValidPair,:);
    mCoh=mean(data,1);
    ee=std(data,0, 1)/sqrt(sum(ValidPair));
%     subplot(1,2,1)
    idxf=f<100;
%     plot(f(idxf), mCoh(idxf), 'color', clr2(ii,:)); hold on;
    errorbar(f(idxf), mCoh(idxf),ee(idxf), 'color', clr2(ii,:)); hold on;
    plot(f(idxf),sig(idxf)*0.3, 'k*'); hold on;
    axis tight;
    set(gca, 'box' ,'off');
    ylabel('Coherence');
    xlabel('Frequency (Hz)');
    
%     subplot(1,2,2)
%     idxf=f>=10&f<100;
% %     plot(f(idxf), mCoh(idxf), 'color', clr2(ii,:)); hold on;
%     errorbar(f(idxf), mCoh(idxf),ee(idxf), 'color', clr2(ii,:)); hold on;
%     plot(f(idxf),sig(idxf), 'k*'); hold on;
%     axis tight;
%     set(gca, 'box' ,'off');
%     ylabel('Coherence');
%     xlabel('Frequency (Hz)');
end

%% plot the population averaged coherence no matter the RF distance
close all;
icue=2;
PairNum=zeros(3);
figure;
clr2=[1 0 0 ; 0 0 1];
for iarea=1:3
    for jarea=1:3
        subplot(3,3,(iarea-1)*3+jarea);
        sig=zeros(1, length(f));
        idxx=cellfun(@(x) ~isempty(find(x==icue, 1))+1,  CuePos{iarea, jarea}, 'UniformOutput', false);
        ValidPair=cell2mat(idxx)==2;
        PairNum(iarea, jarea)=sum(ValidPair);
        for ifreq=1:numel(f)
            data1=cohLFP{iarea, jarea}{icue-1};  %-Coherence_baseline{icue-1};
            data2=cohLFP{iarea, jarea}{icue};  %-Coherence_baseline{icue};
            sig(ifreq)=ttest(data1(ValidPair,ifreq), data2(ValidPair,ifreq));
        end
        for ii=1:2
            data=cohLFP{iarea, jarea}{icue-2+ii}(ValidPair,:);  %-Coherence_baseline{icue-2+ii}(ValidPair,:);   (ValidPair,:)
            mCoh=mean(data,1);
            ee=std(data,0, 1)/sqrt(size(data,1));
            %     subplot(1,2,1)
            idxf=f<100;
            %     plot(f(idxf), mCoh(idxf), 'color', clr2(ii,:)); hold on;
            errorbar(f(idxf), mCoh(idxf),ee(idxf), 'color', clr2(ii,:)); hold on;
            plot(f(idxf),sig(idxf)*0.15, 'k*'); hold on;
            axis tight;
            set(gca, 'box' ,'off');
            ylabel('Coherence');
            xlabel('Frequency (Hz)');
        end
    end
end
%% plot the PSTH for the cells used in coherence analysis
figure;
cellUsed=unique(CellID(ValidPair, 2));
for jj=1:numel(cellUsed)
    subplot(4,4,jj);
    for ii=1:2
        mSpk=mean(ArraySpk{jarea}{icue-2+ii, cellUsed(jj)},2);
        yy=conv(mSpk,ones(1,20)/20, 'same');
        plot(Time, yy, 'color', clr2(3-ii,:)); hold on;
    end
    title(num2str(cellUsed(jj)));
    axis tight;   
end






