%% get the mean response for array responsive unit
clear all;
close all;
folder='Z:\RujiaChen\Results\';
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
for igroup=1:2
    groupName=Group{igroup};    
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end        
        for isexo=0:1
            CueOnResp=cell(1,3);
            ArrayOnResp=cell(1,3);
            RewardOnResp=cell(1,3);
            CueOffResp=cell(1,3);
            ArrayOffResp=cell(1,3);
            RewardOffResp=cell(1,3);
            ElecInfo=cell(1,3);
            for iarea=1:3
                num=0;
                for ifile = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{ifile};
                    load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);                    
                    load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
                    load([folder 'bVisResp_' date '_5sigma.mat']);
                    load([folder 'bMotorResp_' date '_5sigma.mat']);
                    load([folder 'bVisMotor_' date '_5sigma.mat']);
                    load([folder 'CueOn_' date '_5sigma.mat']);
                    load([folder 'RFOn_' date '_5sigma.mat']);
                    load([folder 'SaccadeOn_' date '_5sigma.mat']);
                    if imonkey==1
                        isession=find(strcmp(RawInfo(1,:), date)==1);
                        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
                    else
                        channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
                    end
                    for ich=channelID
                        if igroup==1
                            celltype=bVisResp(ich)*CueOn(ich);
                        elseif igroup==2
                            celltype=bVisMotor(ich)*CueOn(ich);
                        elseif igroup==3
                            celltype=bMotorResp(ich);
                        elseif igroup==4
                            celltype=(bVisResp(ich)+bVisMotor(ich))*CueOn(ich);
                        end
                        if celltype>0
                            num=num+1;
                            if RFOn(ich)>0
                                ipos=ceil(RFOn(ich)/2);
                                icue=2-mod(RFOn(ich),2);
                            elseif SaccadeOn(ich)>0
                                ipos=ceil(SaccadeOn(ich)/2);
                                icue=2-mod(SaccadeOn(ich),2);
                            else
                                ipos=2;
                                icue=imonkey;
                            end
                            
                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            CueOnResp{iarea}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx1,ich,:),1));
                            ArrayOnResp{iarea}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx1,ich,:),1));
                            RewardOnResp{iarea}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx1,ich,:),1));
                            
                            idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            CueOffResp{iarea}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx2,ich,:),1));
                            ArrayOffResp{iarea}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx2,ich,:),1));
                            RewardOffResp{iarea}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx2,ich,:),1));
                            
                            ElecInfo{iarea}(num,:)=[SpikeTrain.channel(ich), ifile, RFOn(ich), SaccadeOn(ich)];
                        end
                    end
                end
            end
            
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp_cueResponsive_5sigma.mat'],'CueOnResp','-v7.3');  %
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp_cueResponsive_5sigma.mat'],'CueOffResp','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp_cueResponsive_5sigma.mat'],'ArrayOnResp','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp_cueResponsive_5sigma.mat'],'ArrayOffResp','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOnResp_cueResponsive_5sigma.mat'],'RewardOnResp','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOffResp_cueResponsive_5sigma.mat'],'RewardOffResp','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ElecInfo_cueResponsive_5sigma.mat'],'ElecInfo','-v7.3');
        end
    end
end

%% plot the cue-on and cue-away population responses 
% close all;
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; 
Time2=-1000:500;
areaName={'LIP','PUL','FEF'};
yplotRange1=[-0.1 0.15];
yplotRange2=[-0.1 0.15];
yplotRange3=[-0.1 0.2];
Nsmooth=50;
Selection={'','_5sigma','_cueResponsive_5sigma'};
iselect=3;
for igroup=2
    groupName=Group{igroup};
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        for isexo=1            
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOnResp' Selection{iselect} '.mat']);
            load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOffResp' Selection{iselect} '.mat']);
            figure;
            clr=[1 0 0; 0 0 1];
            for iarea=1:3
                subplot(3,3,iarea*3-2);
                if ~isempty(CueOnResp{iarea})
                    data1=zeros(size(CueOnResp{iarea}));
                    PeakVal=zeros(size(data1,1),1);
                    baseOn=zeros(size(data1,1),1);
                    for icell=1:size(CueOnResp{iarea},1)
                        data1(icell,:)=conv(CueOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOn(icell)=mean(data1(icell,Time<0&Time>-200));
                        data1(icell,:)=data1(icell,:)-baseOn(icell);
                        mResp=conv(ArrayOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        PeakVal(icell)=max(mResp(Time1>0&Time1<200));
                        if PeakVal(icell)<=0
                            PeakVal(icell)=1;
                        end
                        data1(icell,:)=data1(icell,:)*1000;
                        PeakVal(icell)=PeakVal(icell)*1000;
%                         data1(icell,:)=data1(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data1,1);
                    ee=std(data1,0,1)/sqrt(size(data1,1));
                    idx=find(Time<900&Time>-100);
                    TimePlt=Time(idx);
                    Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(1,:),'linewidth', 1.2); hold on;
                    
                    data2=zeros(size(CueOffResp{iarea}));
                    baseOff=zeros(size(data2,1),1);
                    for icell=1:size(CueOffResp{iarea},1)
                        data2(icell,:)=conv(CueOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOff(icell)=mean(data2(icell,Time<0&Time>-200));
                        data2(icell,:)=data2(icell,:)-baseOff(icell);
                        data2(icell,:)=data2(icell,:)*1000;
%                         data2(icell,:)=data2(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data2,1);
                    ee=std(data2,0,1)/sqrt(size(data2,1));
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(2,:),'linewidth', 1.2); hold on;
                    title([areaName{iarea} '-cue resp']);
                    sigDiff=zeros(numel(idx),1);
                    for jj=1:numel(idx)
                        sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                    end
                    plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.002,'.k','linewidth', 1); hold on;
                    axis tight;
                    %                 ylim(yplotRange1) ;
                    
                    subplot(3,3,iarea*3-1);
                    data1=zeros(size(ArrayOnResp{iarea}));
                    for icell=1:size(ArrayOnResp{iarea},1)
                        data1(icell,:)=conv(ArrayOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data1(icell,:)=data1(icell,:)-baseOn(icell);
                        data1(icell,:)=data1(icell,:)*1000;
%                         data1(icell,:)=data1(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data1,1);
                    ee=std(data1,0,1)/sqrt(size(data1,1));
                    idx=find(Time1<1100&Time1>-200);
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    TimePlt=Time1(idx);
                    Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(1,:),'linewidth', 1.2); hold on;
                    
                    data2=zeros(size(ArrayOffResp{iarea}));
                    for icell=1:size(ArrayOffResp{iarea},1)
                        data2(icell,:)=conv(ArrayOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data2(icell,:)=data2(icell,:)-baseOff(icell);
                        data2(icell,:)=data2(icell,:)*1000;
%                         data2(icell,:)=data2(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data2,1);
                    ee=std(data2,0,1)/sqrt(size(data2,1));
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(2,:),'linewidth', 1.2); hold on;
                    title([areaName{iarea} '-array resp']);
                    sigDiff=zeros(numel(idx),1);
                    for jj=1:numel(idx)
                        sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                    end
                    plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.002,'.k','linewidth', 1); hold on;
                    axis tight;
                    %                 ylim(yplotRange2) ;
                    
                    subplot(3,3,iarea*3);
                    data1=zeros(size(RewardOnResp{iarea}));
                    for icell=1:size(RewardOnResp{iarea},1)
                        data1(icell,:)=conv(RewardOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data1(icell,:)=data1(icell,:)-baseOn(icell);
                        data1(icell,:)=data1(icell,:)*1000;
%                         data1(icell,:)=data1(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data1,1);
                    ee=std(data1,0,1)/sqrt(size(data1,1));
                    idx=find(Time2<100&Time2>-800);
                    TimePlt=Time2(idx);
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(1,:),'linewidth', 1.2); hold on;
                    
                    data2=zeros(size(RewardOffResp{iarea}));
                    for icell=1:size(RewardOffResp{iarea},1)
                        data2(icell,:)=conv(RewardOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data2(icell,:)=data2(icell,:)-baseOff(icell);
                        data2(icell,:)=data2(icell,:)*1000;
%                         data2(icell,:)=data2(icell,:)/PeakVal(icell);
                    end
                    mm=mean(data2,1);
                    ee=std(data2,0,1)/sqrt(size(data2,1));
                    yy1=mm(idx)-ee(idx);
                    yy2=mm(idx)+ee(idx);
                    Vy=[yy1 fliplr(yy2) yy1(1)];
                    patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                    plot(TimePlt, mm(idx),'color',clr(2,:),'linewidth', 1.2); hold on;
                    title([areaName{iarea} '-reward resp']);
                    sigDiff=zeros(numel(idx),1);
                    for jj=1:numel(idx)
                        sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                    end
                    plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.002,'.k','linewidth', 1); hold on;
                    axis tight;
                    %                 ylim(yplotRange3) ;
                end
            end
        end
    end
end

%% converge the data from two monkeys together
close all;
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; 
Time2=-1000:500;
areaName={'LIP','PUL','FEF'};
yplotRange1=[-0.1 0.5];
yplotRange2=[-0.1 0.3];
yplotRange3=[-0.1 0.4];
Nsmooth=50;
Selection={'','_5sigma','_cueResponsive_5sigma'};
iselect=2;
for igroup=3
    groupName=Group{igroup};
    for isexo=1
        Source=cell(2,6);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
            Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
            Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
            Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
            Source{imonkey,5}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOnResp' Selection{iselect} '.mat']);
            Source{imonkey,6}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOffResp' Selection{iselect} '.mat']);
        end
        figure;
        clr=[1 0 0; 0 0 1];
        NN=zeros(1,3);
        for iarea=1:3
            subplot(3,3,iarea*3-2);
            iCndPlot=1;
            data1=[Source{1,1}.CueOnResp{iarea}; Source{2,1}.CueOnResp{iarea}];
            data2=[Source{1,2}.CueOffResp{iarea}; Source{2,2}.CueOffResp{iarea}];
            data3=[Source{1,3}.ArrayOnResp{iarea}; Source{2,3}.ArrayOnResp{iarea}];
            data4=[Source{1,4}.ArrayOffResp{iarea}; Source{2,4}.ArrayOffResp{iarea}];
            data5=[Source{1,5}.RewardOnResp{iarea}; Source{2,5}.RewardOnResp{iarea}];
            data6=[Source{1,6}.RewardOffResp{iarea}; Source{2,6}.RewardOffResp{iarea}];
            NN(iarea)=size(data1,1);
            if ~isempty(data1)
                PeakVal=zeros(size(data1,1),1);
                baseOn=zeros(size(data1,1),1);
                for icell=1:size(data1,1)
                    data1(icell,:)=conv(data1(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOn(icell)=mean(data1(icell,Time<0&Time>-300));
                    data1(icell,:)=data1(icell,:)-baseOn(icell);
                    mResp=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    PeakVal(icell)=max(mResp(Time1>0&Time1<200));
                    if PeakVal(icell)<=0
                        PeakVal(icell)=1;
                    end
                    Amplify=1000;
                    %                 data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
                
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));
                idx=find(Time<900&Time>-300);
                TimePlt=Time(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                baseOff=zeros(size(data2,1),1);
                for icell=1:size(data2,1)
                    data2(icell,:)=conv(data2(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOff(icell)=mean(data2(icell,Time<0&Time>-300));
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
                    %                 data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                mm=mean(data2,1);
                ee=std(data2,0,1)/sqrt(size(data2,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
                title([areaName{iarea} '-cue resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.001,'.k','linewidth', 1); hold on;
                axis tight;
                %                 ylim(yplotRange1) ;
                
                subplot(3,3,iarea*3-1);
                for icell=1:size(data3,1)
                    data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data3(icell,:)=data3(icell,:)-baseOn(icell);
                    %                 data3(icell,:)=data3(icell,:)/PeakVal(icell);
                end
                mm=mean(data3,1);
                ee=std(data3,0,1)/sqrt(size(data3,1));
                idx=find(Time1<1100&Time1>-300);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                TimePlt=Time1(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                for icell=1:size(data4,1)
                    data4(icell,:)=conv(data4(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data4(icell,:)=data4(icell,:)-baseOff(icell);
                    %                 data4(icell,:)=data4(icell,:)/PeakVal(icell);
                end
                mm=mean(data4,1);
                ee=std(data4,0,1)/sqrt(size(data4,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
                title([areaName{iarea} '-array resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data3(:, idx(jj)), data4(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.001,'.k','linewidth', 1); hold on;
                axis tight;
                %                 ylim(yplotRange2) ;
                
                subplot(3,3,iarea*3);
                for icell=1:size(data5,1)
                    data5(icell,:)=conv(data5(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data5(icell,:)=data5(icell,:)-baseOn(icell);
                    %                 data5(icell,:)=data5(icell,:)/PeakVal(icell);
                end
                mm=mean(data5,1);
                ee=std(data5,0,1)/sqrt(size(data5,1));
                idx=find(Time2<300&Time2>-800);
                TimePlt=Time2(idx);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                for icell=1:size(data6,1)
                    data6(icell,:)=conv(data6(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data6(icell,:)=data6(icell,:)-baseOff(icell);
                    %                 data6(icell,:)=data6(icell,:)/PeakVal(icell);
                end
                mm=mean(data6,1);
                ee=std(data6,0,1)/sqrt(size(data6,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
                title([areaName{iarea} '-reward resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data5(:, idx(jj)), data6(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*0.001,'.k','linewidth', 1); hold on;
                axis tight;
                %                 ylim(yplotRange3) ;
            end
        end
    end
end

%% get the d-prime and response latency for population responses 
close all;
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
igroup=4;
groupName=Group{igroup};
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; 
Time2=-1000:500;
areaName={'LIP','PUL','FEF'};
yplotRange1=[-0.1 0.15];
yplotRange2=[-0.1 0.15];
yplotRange3=[-0.1 0.2];
Nsmooth=50;

VisLatency=zeros(2,3);
FitQ=zeros(2,3);
DprimeTemp=cell(2,3);
Selection={'','_5sigma','_cueResponsive_5sigma'};
iselect=3;
for isexo=0:1
    figure;
    clr3=eye(3);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        subplot(1,2,imonkey);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOffResp' Selection{iselect} '.mat']);
        
        dataSource=ArrayOffResp;   %CueOnResp;   %
        TimePoint=Time1/1000;   % Time/1000;  %
        for iarea=1:3
            if ~isempty(dataSource{iarea})
%                 data1=zeros(size(dataSource{iarea}));
%                 PeakVal=zeros(size(data1,1),1);
%                 baseOn=zeros(size(data1,1),1);                 
%                 for icell=1:size(dataSource{iarea},1)
%                     data1(icell,:)=conv(dataSource{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
%                     baseOn(icell)=mean(CueOffResp{iarea}(icell,Time<0&Time>-200));                    
%                     data1(icell,:)=data1(icell,:)-baseOn(icell);
%                     PeakVal(icell)=max(ArrayOnResp{iarea}(icell,Time1>0&Time1<200))-baseOn(icell);
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);                 
%                 end
%                 [param,VisLatency(imonkey,iarea)]=LatencyGauss(TimePoint(TimePoint>-0.2&TimePoint<0.2),mean(data1(:,TimePoint>-0.2&TimePoint<0.2),1),0.33,1,clr3(iarea,:));
%                 FitQ(imonkey,iarea)=param(end);
%                 set(gca,'box','off');
%                 xlim([0 0.15])
%                 save(['Z:\RujiaChen\Results\VisLatency_ArrayOffResp' groupName '.mat'],'VisLatency','-v7.3') ;                   
                
%%%%%%%%%%%%%%%%%get the d-prime between cue-on and cue-away conditions
               data1=zeros(size(ArrayOnResp{iarea}));
                for icell=1:size(ArrayOnResp{iarea},1)
                    data1(icell,:)=conv(ArrayOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOn(icell)=mean(CueOnResp{iarea}(icell,Time<0&Time>-200));   
                    data1(icell,:)=data1(icell,:)-baseOn(icell);
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end               
                
                data2=zeros(size(ArrayOffResp{iarea}));
                baseOff=zeros(size(data1,1),1); 
                for icell=1:size(ArrayOffResp{iarea},1)
                    data2(icell,:)=conv(ArrayOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOff(icell)=mean(CueOffResp{iarea}(icell,Time<0&Time>-200));
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                [DprimeTemp{imonkey,iarea},~,T]=getDPrime(data1,data2,[0.025 0.001],1000);
                plot(T-0.505,DprimeTemp{imonkey,iarea},'color',clr3(iarea,:),'linewidth',1.5); hold on;
                
%                 mm=mean(data1-data2,1);
%                 ee=std(data1-data2,0,1)/sqrt(size(data1,1));
%                 AttResp=downsample(mm,5,2);
%                 AttSte=downsample(ee,5,2);
%                 TimeNew=downsample(Time1/1000,5,2);
%                 patchplot(TimeNew, AttResp,AttSte,  clr3(iarea,:)); hold on;

                line([-0.5 0.5],[0 0],'color','k'); hold on;
                xlim([-0.3 0.5]);
%                 ylim([-0.2 0.4])
                set(gca,'box','off');
                
            end
        end
       
    end
end


%% plot the cue response (0-200) versus attention modulation amplitude 
close all;
TimeCue=-600:1200;
TimeArray=-505:1305;
figure;
Selection={'','_5sigma','_cueResponsive_5sigma'};
iselect=3;
AreaName={'LIP','PUL','FEF'};
isexo=1;
for imonkey=1:2
    monkey=monkeyset{imonkey};
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
    for iarea=1:3
        base1=mean(CueOnResp{iarea}(:,TimeCue>-200&TimeCue<0),2);
        base2=mean(CueOffResp{iarea}(:,TimeCue>-200&TimeCue<0),2);
        peak=max(ArrayOnResp{iarea}(:,TimeArray>0&TimeArray<200),[],2);   %-base1
        mCueOn=(mean(CueOnResp{iarea}(:,TimeCue>0&TimeCue<200),2))-base1;
        mCueOff=(mean(CueOffResp{iarea}(:,TimeCue>0&TimeCue<200),2))-base2;
        CueIdx=(mCueOn-mCueOff)*1000;  %./(mCueOn+mCueOff);
        
        mArrayOn=(mean(ArrayOnResp{iarea}(:,TimeArray>250&TimeArray<500),2))-base1;
        mArrayOff=(mean(ArrayOffResp{iarea}(:,TimeArray>250&TimeArray<500),2))-base2;
        ArrayIdx=(mArrayOn-mArrayOff)*1000;  %./(mArrayOn+mArrayOff);
        
        p=polyfit(CueIdx, ArrayIdx, 1);
        xx=linspace(min(CueIdx), max(CueIdx));
        yy=polyval(p,xx);
        subplot(2,3,imonkey*3-3+iarea);        
        plot(CueIdx, ArrayIdx,'.'); hold on;
        plot(xx,yy,'r','linewidth',1.2); hold on;
        line([min(CueIdx)-5 max(CueIdx)+5], [0 0], 'color','k'); hold on;
        line([0 0], [min(ArrayIdx)-5 max(ArrayIdx)+5], 'color','k'); hold on;
        axis tight;
        xlabel('Relative cue response');
        ylabel('Attention effect');
        title(AreaName{iarea});
        set(gca,'box','off');
        
        %             figure;
        %             hh=histogram(ArrayIdx,-0.4:0.05:0.4); hold on;
        %             mValue=mean(ArrayIdx);
        %             xlabel('Attention effect');
        %             ylabel('Occurence');
        %             line([0 0],[0 max(hh.Values)+5],'color','k'); hold on;
        %             xlim([-0.4 0.4]);
        %             text(0.15, max(hh.Values),num2str(mValue,3));
    end
end

%% plot the responses of individual units
close all;
Time=-600:1200; 
Time1=-505:1305; 
Nsmooth=50;
clr=[1 0 0; 0 0 1];
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
igroup=2;
groupName=Group{igroup};
areaName={'LIP','PUL','FEF'};
for imonkey=2
    monkey=monkeyset{imonkey};
    for isexo=0:1
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp_cueResponsive_5sigma.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp_cueResponsive_5sigma.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp_cueResponsive_5sigma.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp_cueResponsive_5sigma.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ElecInfo_cueResponsive_5sigma.mat']);
        
    
        for iarea=3          
            figure;           
            if ~isempty(CueOnResp{iarea})  
                baseOn=zeros(1,size(CueOnResp{iarea},1));
                baseOff=zeros(1,size(CueOffResp{iarea},1));
                num=0;
                for icell=1:size(CueOnResp{iarea},1)
%                     if ElecInfo{iarea}(icell,3)>0  %==2||ElecInfo{iarea,1}(icell,3)==4
                        data1=conv(CueOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOn(icell)=mean(data1(Time<0&Time>-200));
                        data1=data1-baseOn(icell);
                        PeakCueOn=max(data1(Time>50&Time<200));
                        %                     if PeakCueOn>0.005
                        %                     PeakVal(icell)=max(ArrayOnResp{iarea}(icell,Time1>0&Time1<200))-baseOn(icell);
                        %                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                        
                        data2=conv(CueOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOff(icell)=mean(data2(Time<0&Time>-200));
                        data2=data2-baseOff(icell);
                        
                        num=num+1;
                        subplot(8,8,num);
                        idx=find(Time<900&Time>-100);
                        TimePlt=Time(idx);
                        plot(TimePlt, data1(idx)*1000,'color',clr(1,:),'linewidth', 1.2); hold on;
                        plot(TimePlt, data2(idx)*1000,'color',clr(2,:),'linewidth', 1.2); hold on;
                        title([num2str(ElecInfo{iarea}(icell,1)) '-' num2str(ElecInfo{iarea}(icell,2)) '-' num2str(ElecInfo{iarea}(icell,3))]);
%                     end
                end
                
                num=0;
                figure;
                for icell=1:size(CueOnResp{iarea},1)
%                     if ElecInfo{iarea}(icell,3)>0
                        data1=conv(ArrayOnResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data1=data1-baseOn(icell);
                        PeakCueOn=max(data1(Time>50&Time<200));
                        %                     if PeakCueOn>0.005
                        %                     PeakVal(icell)=max(ArrayOnResp{iarea}(icell,Time1>0&Time1<200))-baseOn(icell);
                        %                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                        
                        data2=conv(ArrayOffResp{iarea}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data2=data2-baseOff(icell);
                        
                        num=num+1;
                        subplot(8,8,num);
                        idx=find(Time1<800&Time1>-300);
                        TimePlt=Time1(idx);
                        plot(TimePlt, data1(idx)*1000,'color',clr(1,:),'linewidth', 1.2); hold on;
                        plot(TimePlt, data2(idx)*1000,'color',clr(2,:),'linewidth', 1.2); hold on;
                        title([num2str(ElecInfo{iarea}(icell,1)) '-' num2str(ElecInfo{iarea}(icell,2)) '-' num2str(ElecInfo{iarea}(icell,3))]);
%                     end
                end
            end                        
        end              
    end
end
%% plot the cue-on and cue-away response for individual units
close all;
iarea=3;
figure;
Nsmooth=50;
TimeCue=-600:1200;
TimeArray=-505:1305;
groupName='VisResp';
Selection={'','_5sigma','_cueResponsive_5sigma'};
iselect=3;
for imonkey=2
    monkey=monkeyset{imonkey};
    for isexo=1
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
%         load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ElecInfo' Selection{iselect} '.mat']);
        num=0;
        for ii=1:size(ArrayOnResp{iarea},1)
%             if ElecInfo{iarea}(ii,3)==1||ElecInfo{iarea}(ii,3)==3
                num=num+1;
                subplot(8,8, num);
                yy=conv(ArrayOnResp{iarea}(ii,:),ones(1,Nsmooth)/Nsmooth,'same');
                zz=conv(CueOnResp{iarea}(ii,:),ones(1,Nsmooth)/Nsmooth,'same');                
%                 idx=TimeArray>-300&TimeArray<600;
%                 yy=yy-mean(zz(TimeCue>-200&TimeCue<0));   %%%% for array response
                %     plot(TimeArray(idx), yy(idx),'r'); hold on;
                idx=TimeCue>-200&TimeCue<600;
                zz=zz-mean(zz(TimeCue>-200&TimeCue<0));   %%% for cue response
                plot(TimeCue(idx), zz(idx),'r'); hold on;
                
                yy=conv(ArrayOffResp{iarea}(ii,:),ones(1,Nsmooth)/Nsmooth,'same');
                zz=conv(CueOffResp{iarea}(ii,:),ones(1,Nsmooth)/Nsmooth,'same');
                yy=yy-mean(zz(TimeCue>-200&TimeCue<0));  % for array response
                %     plot(TimeArray(idx), yy(idx),'b'); hold on;
                zz=zz-mean(zz(TimeCue>-200&TimeCue<0));  %% for cue response
                plot(TimeCue(idx), zz(idx),'b'); hold on;
                axis tight;
                title(num2str(ElecInfo{iarea}(ii,3)));
%             end
            
        end
    end
end

%% scatter plot of spikes in one channel
close all;
date='111817';
load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
Time=-300:1200;
plotRange=[-200 500];
figure;
for ich=10
    idxCh=find(SpikeTrain.channel==ich);
    data=squeeze(SpikeTrain.arrayAlign(:,idxCh,:));
%     subplot(3, 4, ich);
    for ii=1:size(data,1)  %numel(iTrl)  %size(data) %%% scatter plot for all trials in one condition
        mResp=squeeze(data(ii,:));
        idx1=find(Time>plotRange(1)&Time<plotRange(2)&mResp==1);
        plot(Time(idx1), ones(1,numel(idx1))*ii,'k.','linewidth', 1 ); hold on;
    end
    line([0 0],[0 ii+1],'linestyle','--','color','r', 'linewidth', 1.5);
    xlabel('Time relative to array onset (ms)')
    ylabel('Trial number');
end

%% compare the neural responses (spikes) in cue-on and cue-away conditions 
close all;
clear all;
folder='Z:\RujiaChen\Results\';
cueTy='exo';
monkey='Vasco';  %'Mikey';
Unit='ArrayResp_unSorted';  % for spike
Unit1='ArrayResp_corrected_by_N-1';  % for LFP
load([folder  monkey 'CueSpk_' cueTy '_' Unit '.mat']);
load([folder  monkey 'ArraySpk_' cueTy '_' Unit '.mat']);
TimeCue=-600:1200; 
TimeArray=-900:900;
TimeCoh=[-300 800];

%% plot the mean spike activities of units in the same area for each session on each condition
Nsession=size(CueSpk,1);
% Nsession=size(ArraySpk,1);
Source=CueSpk;
figure;
RFloc=cell(size(Source,1), size(Source, 2));
for idate=1:Nsession
    for iarea=1:3
%         subplot(Nsession, 3, (idate-1)*3+iarea);
%         for icnd=1:size(Source{idate,iarea},1)
%             mm=zeros(size(Source{idate,iarea},2), size(Source{idate,iarea}{1},1));
%             for icell=1:size(Source{idate,iarea},2)
%                 mm(icell,:)=mean(Source{idate,iarea}{icnd,icell}, 2);                
%             end
%             mResp=mean(mm);
%             plot(1:length(mResp)-19, mResp(10:end-10)); hold on;
%         end

%%%%%% identify the location with strongest response
        if ~isempty(Source{idate, iarea})
            for icell=1:size(Source{idate, iarea}, 2)
                mResp=zeros(1,4);
                for icnd=1:size(Source{idate, iarea}, 1)
                    mm=mean(Source{idate, iarea}{icnd, icell}, 2);
                    idxT=TimeCue>0&TimeCue<200;
                    idxBase=TimeCue>-200&TimeCue<0;
                    mResp(icnd)=mean(mm(idxT))-mean(mm(idxBase));
                end                
                [~ , RFloc{idate, iarea}(icell)]=max(mResp);
            end
        end
    end    
end

save([folder 'RFloc_' monkey '.mat'], 'RFloc');

%% get the d-prime of each cell at each location

TimeCue=-600:1200; % for sessions after 092718
TimeArray=-900:900;
TimeCoh=[-300 800];
DataSource= ArraySpk;   %CueSpk;  
Nsession=size(DataSource,1);
for idate=1:Nsession
    TimeIdx=TimeCue>-100&TimeCue<=600;
%         TimeIdx=TimeArray>-100&TimeArray<=600;
    for iarea=1:3
        if ~isempty(DataSource{idate,iarea})
            for ipos=1:size(DataSource{idate,iarea})/2
                for icell=1:size(DataSource{idate,iarea},2)
                    if ~isempty(DataSource{idate,iarea}{ipos*2-1, icell})&&~isempty(DataSource{idate,iarea}{ipos*2, icell})
                        data1=transpose(DataSource{idate,iarea}{ipos*2-1, icell});   %CueSpk
                        data2=transpose(DataSource{idate,iarea}{ipos*2, icell});
                        [dprime{idate, iarea}(ipos, icell,:), ~, tt]= getDPrime(data1(:,TimeIdx),data2(:,TimeIdx),[0.03 0.01],1000);
                    end
                end
            end
        end
    end    
end
%%
close all;
clear dprimeAll;
AreaName={'LIP', 'PUL', 'FEF'};
dateOrder=[10,8,2,5,6,7,3, 1,4,9, 11:29 ];
   figure;
for iarea=1:3
    subplot(1,3, iarea);
%    for ipos=1:2
%        subplot(1,2,ipos);
        dprimeAll=[];
        for iorder=1:size(dprime,1)
            idate=dateOrder(iorder);
            if ~isempty(dprime{idate, iarea})
                [n1, n2, n3]= size(dprime{idate, iarea});
                %                 dprimeAll=[dprimeAll; reshape(-dprime{idate, iarea}(ipos,:,:), n2, n3)];
                nn=squeeze(mean(dprime{idate, iarea},3));
                mm=reshape(nn,n1,n2);
                [~, maxPos]=max(mm, [], 1);
                MaxDprime=zeros(n2, n3);
                for icell=1:numel(maxPos)
                    MaxDprime(icell,:)=squeeze(dprime{idate, iarea}(maxPos(icell),icell,:));
                end
                dprimeAll=[dprimeAll; MaxDprime];
            end
            if iorder==10
                border = size(dprimeAll, 1) ;
            end
        end
        nWin=10;
        SmoothKernal=ones(nWin)/nWin^2;
        dprimeAll=conv2(dprimeAll, SmoothKernal, 'same');
        imagesc(tt-0.1, 1:size(dprimeAll, 1), dprimeAll);
        line([min(tt)-0.1 max(tt)-0.1], [border border], 'color', 'k'); hold on;
        title(AreaName{iarea});
        caxis([0 1]);
        colorbar;
        %    end
end

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

%%  plot the distribution of Cue distinguished responsive units and attentional modulated units
close all;
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};
dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey 
folder='Z:\RujiaChen\Results\';
figure;
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
load('Z:\RujiaChen\Results\Mikey_RecordingDepth.mat');
% load([folder 'RelativeDepth.mat']);
yy=1:96;
areaName={'LIP','PUL','FEF'};
count=zeros(1,3);
count1=zeros(1,3);
count2=zeros(1,3);
for ifile=1:numel(dateUsed)
    date=dateUsed{ifile};
    load([folder 'CueOn_' date '.mat']);
    load([folder 'bVisResp_' date '.mat']);
    load([folder  'Endo_AttentionEffect_' date '.mat']);   %ArrayDelay_
    load([folder 'Endo_CueFlag_' date '.mat']);   %ArrayDelay_
    load([folder 'bArrayResp_' date '.mat']);
    isession=find(strcmp(RawInfo(1,:), date)==1);
    
%     %%%%% plot the distribution for two positions seperately
%     for ii=1:2
%         xx=(3*ifile-3+ii)*ones(1,96);
% %         xx=ifile*ones(1,96);
% %         idx=CueOn(:,ii)==0;
%         idx=AttentionEffect(:,ii)==0;
%         plot(xx(idx), yy(idx),'ko','markersize', 5, 'linewidth',0.5); hold on;
% 
%         idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii);  %CueOn(:,ii)>0&
%         plot(xx(idx), yy(idx),'ro','markerfacecolor','r','markersize', 5, 'linewidth',0.5); hold on;
%         idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)~=CueFlag(:,ii);  %CueOn(:,ii)>0&
%         plot(xx(idx), yy(idx),'bo','markerfacecolor','b','markersize', 5, 'linewidth',0.5); hold on;
% %         idx=CueOn(:,ii)>0;  %&AttentionEffect(:,ii)==0;
% %         plot(xx(idx), yy(idx),'ko','markerfacecolor','k','markersize', 5, 'linewidth',0.5); hold on;            
%     end        
    
     %%%% plot the distribution pool two positions together   
     bResp=sum(CueOn,2);
     bArrResp=sum(bArrayResp,2);
     Facilitate=zeros(numel(bResp),2);
    for ii=1:2
        idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii);
        Facilitate(idx,ii)=1;             
    end
    bFacilitate=sum(Facilitate, 2);
    yProbe=[4.65:-0.15:0; 4.65:-0.15:0; 2.325:-0.075:0];
    nUnits=[numel(RawInfo{2,isession}), numel(RawInfo{3,isession}), numel(RawInfo{4,isession})];  % for Mikey
    unitID=[0, nUnits(1), nUnits(1)+nUnits(2)];
    for iarea=1:3
        subplot(1,3,iarea);
%         xx=ifile*ones(1,96);
%         yy= yProbe(iarea,:)-RelativeDepth(iarea,ifile);   %depth(iarea,ifile)  % for Vasco
%         idx=bResp(iarea*32-31:iarea*32)==0;
        
        xx=ifile*ones(1,nUnits(iarea));
        yy=[nUnits(iarea):-1:1]*0.15-0.15;
        yy=yy-Depth{iarea+1, isession};
        idxx=1:nUnits(iarea);  % for Mikey
        idx=bArrResp(idxx+unitID(iarea))==0;  %&bArrResp(iarea*32-31:iarea*32)>0          bResp
%         idx=bFacilitate(idxx+unitID(iarea))==0;  
        plot(xx(idx), yy(idx),'ks','markersize', 6, 'linewidth',0.5); hold on;

%         idx=bResp(iarea*32-31:iarea*32)>0;
        idx=bArrResp(idxx+unitID(iarea))>0;    %bResp
%         idx=bFacilitate(idxx+unitID(iarea))>0;    
        count(iarea)=count(iarea)+sum(idx);
        plot(xx(idx), yy(idx),'ks','markerfacecolor','k','markersize', 6, 'linewidth',0.5); hold on;   
        axis tight;
%         xlim([0 numel(dateUsed)+1]);
%         ylim([0 33]);
        title(areaName{iarea});
        set(gca,'box', 'off');
    end   
end

% %%%%%% plot the overlap of exo_ and endo_ facilitated units
% for ifile=1:numel(dateUsed)
%     date=dateUsed{ifile};
%     load([folder 'CueOn_' date '.mat']);
%     load([folder  'Endo_AttentionEffect_' date '.mat']);   %ArrayDelay_
%     load([folder  'Endo_CueFlag_' date '.mat']);   %ArrayDelay_
%     isession=find(strcmp(RawInfo(1,:), date)==1);
%     nUnits=[numel(RawInfo{2,isession}), numel(RawInfo{3,isession}), numel(RawInfo{4,isession})];  % for Mikey
%     unitID=[0, nUnits(1), nUnits(1)+nUnits(2)];
%     
%      %%%% plot the distribution pool two positions together   
%      Facilitate=zeros(size(AttentionEffect,1),2);
%      for ii=1:2
%          idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii);
%          Facilitate(idx,ii)=1;
%      end
%      bFacilitateEndo=sum(Facilitate, 2);
%      
%      for iarea=1:3
%          subplot(1,3,iarea);
%          xx=ifile*ones(1,nUnits(iarea));
%          yy=[nUnits(iarea):-1:1]*0.15-0.15;
%          yy=yy-Depth{iarea+1, isession};
%          idxx=1:nUnits(iarea);  % for Mikey
%          idx=bFacilitateEndo(idxx+unitID(iarea))==0;
%          plot(xx(idx), yy(idx),'ks','markersize', 6, 'linewidth',0.5); hold on;    
%          idx=bFacilitateEndo(idxx+unitID(iarea))>0;
%          plot(xx(idx), yy(idx),'rs','markerfacecolor','r','markersize', 6, 'linewidth',0.5); hold on;
%      end
%      
%      load([folder  'AttentionEffect_' date '.mat']);   %ArrayDelay_
%      load([folder  'CueFlag_' date '.mat']);   %ArrayDelay_
%      Facilitate=zeros(size(AttentionEffect,1),2);
%      for ii=1:2
%          idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii);
%          Facilitate(idx,ii)=1;
%      end
%      bFacilitateExo=sum(Facilitate, 2);
%      
%      for iarea=1:3
%          subplot(1,3,iarea);
%          xx=ifile*ones(1,nUnits(iarea));
%          yy=[nUnits(iarea):-1:1]*0.15-0.15;
%          yy=yy-Depth{iarea+1, isession};
%          idxx=1:nUnits(iarea);  % for Mikey
% %          yy=1:nUnits(iarea);  % for Mikey
%          
%          idx=bFacilitateExo(idxx+unitID(iarea))>0&bFacilitateEndo(idxx+unitID(iarea))==0;
%          plot(xx(idx), yy(idx),'bs','markerfacecolor','b','markersize', 6, 'linewidth',0.5); hold on;
%          idx=bFacilitateExo(idxx+unitID(iarea))>0&bFacilitateEndo(idxx+unitID(iarea))>0;
%          count(iarea)=count(iarea)+sum(idx);
%          plot(xx(idx), yy(idx),'ks','markerfacecolor','k','markersize', 6, 'linewidth',0.5); hold on;
%          
% %          xlim([0 numel(dateUsed)+1]);
% %          ylim([0 33]);
%          title(areaName{iarea});
%          set(gca,'box', 'off');
%      end
% end


% %%%%%% merge the facilated units during both cue delay and array delay
% %%%%%% period together
% Facilitate=zeros(numel(dateUsed), 2, 96);
% CuePosition=zeros(numel(dateUsed), 2, 96);
% for ifile = 1:numel(dateUsed)
%     date=dateUsed{ifile};
%     load([folder  'ArrayDelay_AttentionEffect_' date '.mat']);
%     load([folder 'ArrayDelay_CueFlag_' date '.mat']);
%     isession=find(strcmp(RawInfo(1,:), date)==1);
%     
%     for ii=1:2
%         idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii);
%         Facilitate(ifile, ii, idx)=1;
%         CuePosition(ifile, ii, idx)=CueFlag(idx,ii);
%         idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)~=CueFlag(:,ii);
%         Facilitate(ifile, ii, idx)=-1;
%         CuePosition(ifile, ii, idx)=CueFlag(idx,ii);
%     end
%     
%     load([folder  'AttentionEffect_' date '.mat']);
%     load([folder 'CueFlag_' date '.mat']);
%      for ii=1:2
%          idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)==CueFlag(:,ii); 
%          Facilitate(ifile, ii, idx)=Facilitate(ifile, ii, idx)+1;
%          CuePosition(ifile, ii, idx)=CueFlag(idx, ii);
%          idx=AttentionEffect(:,ii)>0&AttentionEffect(:,ii)~=CueFlag(:,ii); 
%          Facilitate(ifile, ii, idx)=Facilitate(ifile, ii, idx)-1;
%          CuePosition(ifile, ii, idx)=CueFlag(idx, ii);
%          
%          xx=(3*ifile-3+ii)*ones(1,96);
%          idx=squeeze(Facilitate(ifile, ii, :))==0;
%          plot(xx(idx), yy(idx),'ks','markersize', 5, 'linewidth',0.5); hold on;
%          idx=squeeze(Facilitate(ifile, ii, :))>0;
%          plot(xx(idx), yy(idx),'rs', 'markerfacecolor','r','markersize',5, 'linewidth',0.5); hold on;
%          idx=squeeze(Facilitate(ifile, ii, :))<0;
%          plot(xx(idx), yy(idx),'bs','markerfacecolor','b','markersize', 5, 'linewidth',0.5); hold on; 
%      end
% end
% save([folder 'CuePosition.mat'], 'CuePosition', '-v7.3');
% save([folder 'Facilitate_units.mat'], 'Facilitate', '-v7.3');
%% plot the distribution of array responsive units
% close all;
figure;
folder='Z:\RujiaChen\Results\';
load([folder 'RecordingDepth.mat']);
load('Z:\RujiaChen\Results\Regress_location_depth_params.mat');
depth=zeros(3,numel(dateUsed));
AreaName={'LIP','PUL','FEF'};
RelativeDepth=zeros(3,numel(dateUsed));
for iarea=1:3
    for ii=1:numel(dateUsed)
        for jj=iarea:3:numel(dataRaw{1})
            if strcmp(dateUsed{ii},dataRaw{1}{jj})&&strcmp(dataRaw{2}{jj},AreaName{iarea})
                depth(iarea, ii)=dataRaw{3}(jj);
                xloc=dataRaw{4}(jj);
                yloc=dataRaw{5}(jj);
                zz=params(iarea,1) + params(iarea,2)*xloc + params(iarea,3)*yloc ;
                RelativeDepth(iarea, ii)=depth(iarea, ii)-zz;
            end
        end
    end
end
clr4=colormap(jet(4));
count=zeros(1,3);
yProbe=[4.65:-0.15:0; 4.65:-0.15:0; 2.325:-0.075:0];
for ifile=[1:5 7:9]    %numel(dateUsed)
    date=dateUsed{ifile};
    load([folder 'CueOn_' date '.mat']);
    load([folder 'bVisResp_' date '.mat']);
    load([folder 'bArrayResp_' date '.mat']);
    load([folder 'Resp2cue_' date '.mat']);  % cueMeanResp:  response to cue stimuli at the flanker task
    
    bResp=sum(bArrayResp,2)>0;
    for iarea=1:3
        subplot(1,3,iarea);
        load([folder 'EstimateRange_' AreaName{iarea} '.mat']);
        xx=ifile*ones(1,32);
%         yy=1:32;
%         yy=31:-1:0;
        yy=yProbe(iarea,:)-RelativeDepth(iarea,ifile);   %depth(iarea,ifile);   %
        idx=bResp(iarea*32-31:iarea*32)==0&EstimateRange(ifile,:)'==1;
        plot(xx(idx), yy(idx),'ks','markersize', 6, 'linewidth',0.5); hold on;
       
        RespRatio{iarea}(:, ifile)=max(cueMeanResp(iarea*32-31:iarea*32,:),[],2)./mean(cueMeanResp(iarea*32-31:iarea*32,:),2);
        RespRatio{iarea}(idx, ifile)=1;
        idx=EstimateRange(ifile,:)'==1;  %bResp(iarea*32-31:iarea*32)>0&
        count(iarea)=count(iarea)+sum(idx);
        plot(xx(idx), yy(idx),'ks','markerfacecolor','k','markersize', 6, 'linewidth',0.5); hold on;
        
        
%         [~,maxIdx1]=max(cueMeanResp(iarea*32-31:iarea*32,:),[],2); hold on;
%         for iRf=1:4
%             idx=bResp(iarea*32-31:iarea*32)>0&maxIdx1==iRf;           
%             plot(xx(idx), yy(idx),'s','markeredgecolor',clr4(iRf, :), 'markerfacecolor',clr4(iRf, :),'markersize', 6, 'linewidth',0.5); hold on; 
%         end
        
        title(AreaName{iarea});
        axis tight;
        xlim([0 11]);
        
    end
end
% figure;
% for iarea=1:3
%     subplot(1,3,iarea);
%     imagesc(1:10, 31:-1:0, RespRatio{iarea});
%     axis xy;
% end

save([folder 'RelativeDepth.mat'],'RelativeDepth','-v7.3');



%% plot the population averaged responses
close all;
clear all;
plotRange=[-150 600];
TargetPosY=[200,-200];
Time=-300:1200;
clr=[1 0 0; 0 0 1];
area={'LIP','PUL','FEF'};
figure;
for iarea=1:3
    load(['Z:\RujiaChen\Results\MUA_exoCueResp_' area{iarea} '.mat']);
%     figure;  
    for ipos=1:2
        subplot(3,2,iarea*2+ipos-2);
        CueOn=[];
        CueOff=[];
        ArrayOn=[];
        ArrayOff=[];
%         subplot(1,2,ipos);
        for ifile=1:size(MUA.cueOn,1)
            if ~isempty(MUA.cueOn{ifile,ipos})
                CueOn  = [CueOn;squeeze(MUA.cueOnEndo{ifile,ipos})]; % population response to array / cue stimulus
                ArrayOn  = [ArrayOn;squeeze(MUA.arrayOnEndo{ifile,ipos})]; % population response to array / cue stimulus
%                 ArrayOn  = [ArrayOn;squeeze(MUA.arrayOn{ifile,ipos})]; % population response to array / cue stimulus
%                 CueOn  = [CueOn;squeeze(MUA.cueOn{ifile,ipos})]; % population response to array / cue stimulus
            end
            if ~isempty(MUA.cueOff{ifile,ipos})
                CueOff = [CueOff;squeeze(MUA.cueOffEndo{ifile,ipos})];
                ArrayOff = [ArrayOff;squeeze(MUA.arrayOffEndo{ifile,ipos})];
%                  CueOff = [CueOff;squeeze(MUA.cueOff{ifile,ipos})];
%                  ArrayOff = [ArrayOff;squeeze(MUA.arrayOff{ifile,ipos})];
            end
        end
        NormOn=zeros(size(CueOn));
        NormOff=zeros(size(CueOn));
        for icell=1:size(ArrayOn,1)
%             NormOn(icell,:)=ArrayOn(icell,:)-mean(CueOn(icell,Time>-250&Time<0)); % all condition subtract their own baseline(-250 to 0 ms, before cue onset)
            NormOn(icell,:)=CueOn(icell,:)-mean(CueOn(icell,Time>-250&Time<0));
            mm=conv(NormOn(icell,:),ones(1,20)/20,'same');
            NormOn(icell,:)= mm/max(mm(300:700));
            
%             NormOff(icell,:)=ArrayOff(icell,:)-mean(CueOff(icell,Time>-250&Time<0));
            NormOff(icell,:)=CueOff(icell,:)-mean(CueOff(icell,Time>-250&Time<0));
            NormOff(icell,:)=conv(NormOff(icell,:),ones(1,20)/20,'same');
            NormOff(icell,:)= NormOff(icell,:)/max( mm(300:700));
        end
        mPSTH{1}=mean(NormOn,1);
        mPSTH{2}=mean(NormOff,1);
        eePSTH{1}=std(NormOn,0,1)/sqrt(size(NormOn,1));
        eePSTH{2}=std(NormOff,0,1)/sqrt(size(NormOff,1));
        
        for ii=1:2
            idx1=find(ismember(Time,plotRange(1):5:plotRange(2)));
            TimePlt=Time(idx1);
            %         hh=errorbar(Time(idx1), mPSTH{ii}(idx1),eePSTH{ii}(idx1),'color',clr(ii,:),'linewidth', 1.2); hold on;
            Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
            yy1=mPSTH{ii}(idx1)-eePSTH{1}(idx1);
            yy2=mPSTH{ii}(idx1)+eePSTH{1}(idx1);
            Vy=[yy1 fliplr(yy2) yy1(1)];
            patch(Vx,Vy,clr(ii,:),'edgecolor',clr(ii,:),'facealpha',0.5,'edgealpha',0.5); hold on;
            
            hh(ii)=plot(TimePlt, mPSTH{ii}(idx1),'color',clr(ii,:),'linewidth', 1.2); hold on;
%             if ipos==2
%                 set(hh(ii),'linestyle','--');
%             else
%                 set(hh(ii),'linestyle','-');
%             end
        end
        SigResp=zeros(1,length(idx1));
        for itime=1:numel(idx1)
            SigResp(itime)= ttest(NormOn(:,idx1(itime)), NormOff(:,idx1(itime)));
        end
        flag=0*ones(1,length(TimePlt));
        plot(TimePlt(SigResp==1),flag(SigResp==1),'ko','markerfacecolor','k','markersize',3); hold on;
        legend(hh,'On','Off','location','northeast');
%         title(['(200,' num2str(TargetPosY(ipos)) ')']);
        title(area{iarea});
        axis tight ;
        %     ylim([-0.01 0.04]);
    end
end








