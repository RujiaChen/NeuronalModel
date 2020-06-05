%% get the mean response for different groups
clear all;
close all;
folder='Z:\Rujia\Results\';
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
load([folder 'Mikey_RecordingInfo.mat']);
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
CueCondition={'Endo','Exo'};
TimeArray=-505:1305; 
for igroup=4
    groupName=Group{igroup};   
    for imonkey=2  
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end
        
        for isexo=1
            CueOnResp=cell(3,2);
            ArrayOnResp=cell(3,2);
            RewardOnResp=cell(3,2);
            CueOffResp=cell(3,2);
            ArrayOffResp=cell(3,2);
            RewardOffResp=cell(3,2);
            ElecInfo=cell(3,2);  
            SaccadeOnResp=cell(3,1);
            SaccadeOffResp=cell(3,1);
            
            for iarea=1:3
                num=0;
                count=[0 0];
                for ifile = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{ifile};
                    load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                    load([folder 'CorrectTrialParam_' date '.mat']);
                    load([folder 'bVisResp_' date '_SingleUnit.mat']);
                    load([folder 'bMotorResp_' date '_SingleUnit.mat']);
                    load([folder 'bVisMotor_' date '_SingleUnit.mat']);
                    load([folder 'CueOn_' date '_SingleUnit.mat']);
                    load([folder 'RFOn_' date '_SingleUnit.mat']);
                    load([folder 'SaccadeOn_' date '_SingleUnit.mat']);
                    load([folder 'SaccadeTime_ToArray_' CueCondition{isexo+1} '_' date '.mat']);
                    if imonkey==1
                        isession=find(strcmp(RawInfo(1,:), date)==1);   
                        channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
                    else
                        channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
                    end
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

                            idx1=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
                            SpikeTrain_group.CueOn{iarea}{num}=squeeze(SpikeTrain.cueAlign(idx1,ich,:));
                            SpikeTrain_group.ArrayOn{iarea}{num}=squeeze(SpikeTrain.arrayAlign(idx1,ich,:));
                            
                            idx2=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
                            SpikeTrain_group.CueOff{iarea}{num}=squeeze(SpikeTrain.cueAlign(idx2,ich,:));
                            SpikeTrain_group.ArrayOff{iarea}{num}=squeeze(SpikeTrain.arrayAlign(idx2,ich,:));  
                                                        
%                             CueOnResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx1,ich,:),1));   
%                             ArrayOnResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx1,ich,:),1));  
%                             RewardOnResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx1,ich,:),1)); 
%                             ntrl=0;
%                             motorResp=[];
%                             if numel(idx1)>=10 
%                                 count(1)=count(1)+1;
%                                 for itrl=1:numel(idx1)
%                                     if SaccadeTime{ipos, icue}(itrl)>400&&SaccadeTime{ipos, icue}(itrl)<=1100
%                                         ntrl=ntrl+1;
%                                         idxT=TimeArray-SaccadeTime{ipos, icue}(itrl)>=-800&TimeArray-SaccadeTime{ipos, icue}(itrl)<=200;
%                                         motorResp(ntrl,:)=squeeze(SpikeTrain.arrayAlign(idx1(itrl),ich,idxT));  
%                                     end
%                                 end
%                             end
%                             SaccadeOnResp{iarea}(count(1),:)=mean(motorResp,1);
                            
%                             idx2=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
%                             CueOffResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx2,ich,:),1));
%                             ArrayOffResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx2,ich,:),1));  
%                             RewardOffResp{iarea,1}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx2,ich,:),1));  
%                             ntrl=0;
%                             motorResp=[];
%                             if numel(idx2)>=10
%                                 count(2)=count(2)+1;
%                                 for itrl=1:numel(idx2)
%                                     if SaccadeTime{ipos, 3-icue}(itrl)>400&&SaccadeTime{ipos, 3-icue}(itrl)<=1100
%                                         ntrl=ntrl+1;
%                                         idxT=TimeArray-SaccadeTime{ipos, 3-icue}(itrl)>=-800&TimeArray-SaccadeTime{ipos, 3-icue}(itrl)<=200;
%                                         motorResp(ntrl,:)=squeeze(SpikeTrain.arrayAlign(idx2(itrl),ich,idxT));
%                                     end
%                                 end
%                             end
%                             SaccadeOffResp{iarea}(count(2),:)=mean(motorResp,1);
%                             ElecInfo{iarea,1}(num,:)=[SpikeTrain.channel(ich), ifile, RFOn(ich)];
%                             %%%%%%%% another position
%                             jpos=3-ipos;
%                             icue=imonkey;
%                             idx1=CorrectTrialParam.rfyi==TargetPosY(jpos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
%                             CueOnResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx1,ich,:),1));
%                             ArrayOnResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx1,ich,:),1));
%                             RewardOnResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx1,ich,:),1));
% 
%                             idx2=CorrectTrialParam.rfyi==TargetPosY(jpos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
%                             CueOffResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.cueAlign(idx2,ich,:),1));
%                             ArrayOffResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.arrayAlign(idx2,ich,:),1));
%                             RewardOffResp{iarea,2}(num,:)=squeeze(mean(SpikeTrain.rewardAlign(idx2,ich,:),1));
                        end
                    end
                end
            end
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeTrain_cueResponsive_SU.mat'],'SpikeTrain_group','-v7.3');
            clear SpikeTrain_group;
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp_cueResponsive_SU.mat'],'CueOnResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp_cueResponsive_SU.mat'],'CueOffResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp_cueResponsive_SU.mat'],'ArrayOnResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp_cueResponsive_SU.mat'],'ArrayOffResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOnResp_cueResponsive_SU.mat'],'RewardOnResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_RewardOffResp_cueResponsive_SU.mat'],'RewardOffResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOnResp_cueResponsive_SU.mat'],'SaccadeOnResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOffResp_cueResponsive_SU.mat'],'SaccadeOffResp','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_ElecInfo_cueResponsive_SU.mat'],'ElecInfo','-v7.3');
        end
    end
end

%% plot the cue-on and cue-away population responses 
close all;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %CueOnResp;
% Time1=-905:905;
Time1=-505:1305;  %ArrayOnResp
Time2=-1000:500; % RewardOnResp
Time3=-800:200; % RewardOnResp
areaName={'LIP','PUL','FEF'};
yplotRange=[-1 10];

Nsmooth=50;
Selection={'','_SU','_cueResponsive_SU'};
iselect=3;
igroup=2;
groupName=Group{igroup};
for imonkey=2
    monkey=monkeyset{imonkey};
    for isexo=1
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOnResp' Selection{iselect} '.mat'])
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOffResp' Selection{iselect} '.mat'])
        figure;
        clr=[1 0 0; 0 0 1];
        for iarea=1            
            iCndPlot=1;
            if ~isempty(CueOnResp{iarea,iCndPlot})
                data1=zeros(size(CueOnResp{iarea,iCndPlot}));
                PeakVal=zeros(size(data1,1),1);
                baseOn=zeros(size(data1,1),1);
                for icell=1:size(CueOnResp{iarea,iCndPlot},1)
                    data1(icell,:)=conv(CueOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOn(icell)=mean(data1(icell,Time<0&Time>-300));
                    data1(icell,:)=data1(icell,:)-baseOn(icell);
                    mResp=conv(ArrayOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    PeakVal(icell)=max(mResp(Time1>0&Time1<200));
                    if PeakVal(icell)<=0
                        PeakVal(icell)=1;
                    end
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
             
                subplot(1,3,iarea*3-2);
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));
                idx=find(Time<=600&Time>=-300);
                TimePlt=Time(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Amplify=1000; 
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                data2=zeros(size(CueOffResp{iarea,iCndPlot}));
                baseOff=zeros(size(data2,1),1);
                for icell=1:size(CueOffResp{iarea,iCndPlot},1)
                    data2(icell,:)=conv(CueOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOff(icell)=mean(data2(icell,Time<0&Time>-300));
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                mm=mean(data2,1);
                ee=std(data2,0,1)/sqrt(size(data2,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
%                 title([areaName{iarea} '-cue resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                end
%                 plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.001,'.k','linewidth', 1); hold on;
                axis tight;
                yticks([0 5 10]);
                xticks(-200:200:600);
                ylim(yplotRange) ;
                set(gca,'linewidth',2,'fontsize',20);                
                
                subplot(1,3,iarea*3-1);
                data1=zeros(size(ArrayOnResp{iarea,iCndPlot}));
                for icell=1:size(ArrayOnResp{iarea,iCndPlot},1)
                    data1(icell,:)=conv(ArrayOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data1(icell,:)=data1(icell,:)-baseOn(icell);
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));                
                idx=find(Time1<=400&Time1>=-300);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                TimePlt=Time1(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                data2=zeros(size(ArrayOffResp{iarea,iCndPlot}));
                for icell=1:size(ArrayOffResp{iarea,iCndPlot},1)
                    data2(icell,:)=conv(ArrayOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                mm=mean(data2,1);
                ee=std(data2,0,1)/sqrt(size(data2,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
%                 title([areaName{iarea} '-array resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                end
%                 plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-0.001,'.k','linewidth', 1); hold on;
                axis tight;
                yticks([0 5 10]);
                xticks(-200:200:400);
                ylim(yplotRange) ;
                set(gca,'linewidth',2,'fontsize',20);
                
                subplot(1,3,iarea*3);
                data1=zeros(size(SaccadeOnResp{iarea,iCndPlot}));
                for icell=1:size(SaccadeOnResp{iarea,iCndPlot},1)
                    data1(icell,:)=conv(SaccadeOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data1(icell,:)=data1(icell,:)-baseOn(icell);
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));
                idx=find(Time3>=-400&Time3<=160);
                TimePlt=Time3(idx);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                data2=zeros(size(SaccadeOffResp{iarea,iCndPlot}));
                for icell=1:size(SaccadeOffResp{iarea,iCndPlot},1)
                    data2(icell,:)=conv(SaccadeOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                mm=mean(data2,1);
                ee=std(data2,0,1)/sqrt(size(data2,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
%                 title([areaName{iarea} '-reward resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                end
%                 plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*0.001,'.k','linewidth', 1); hold on;
%                 axis tight;
                yticks([0 5 10]);
                xticks(-400:200:200);
                ylim(yplotRange);
                set(gca,'linewidth',2,'fontsize',20);
            end
        end
    end
end

%% converge the data from two monkeys together
close all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; 
Time2=-1000:500;
Time3=-800:200; % RewardOnResp
areaName={'LIP','PUL','FEF'};
yplotRange=[-1 5.5;-1 9; -2 8];
% yplotRange2=[-2 7];
% yplotRange3=[-2 7];
Nsmooth=50;
Selection={'_','_SU','_cueResponsive_SU'};
for igroup=4
    groupName=Group{igroup};
    if igroup==3
        iselect=2;
    else
        iselect=3;
    end
    iselect1=3;
    for isexo=1
        Source=cell(2,6);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
            Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
            Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
            Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
            Source{imonkey,5}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOnResp' Selection{iselect1} '.mat']);
            Source{imonkey,6}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOffResp' Selection{iselect1} '.mat']);
        end
%         figure;
        clr=[1 0 0; 0 0 1];
        NN=zeros(1,3);
        for iarea=1
            subplot(3,3,iarea*3-2);
            iCndPlot=1;
            data1=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];
            data2=[Source{1,2}.CueOffResp{iarea,iCndPlot}; Source{2,2}.CueOffResp{iarea,iCndPlot}];
            data3=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}; Source{2,3}.ArrayOnResp{iarea,iCndPlot}];
            data4=[Source{1,4}.ArrayOffResp{iarea,iCndPlot}; Source{2,4}.ArrayOffResp{iarea,iCndPlot}];
            data5=[Source{1,5}.SaccadeOnResp{iarea,iCndPlot}; Source{2,5}.SaccadeOnResp{iarea,iCndPlot}];
            data6=[Source{1,6}.SaccadeOffResp{iarea,iCndPlot}; Source{2,6}.SaccadeOffResp{iarea,iCndPlot}];
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
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
                
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));
                idx=find(Time<600&Time>-300);
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
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
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
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:)) ;
                
                subplot(3,3,iarea*3-1);
                for icell=1:size(data3,1)
                    data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data3(icell,:)=data3(icell,:)-baseOn(icell);
%                     data3(icell,:)=data3(icell,:)/PeakVal(icell);
                end
                mm=mean(data3,1);
                ee=std(data3,0,1)/sqrt(size(data3,1));
                idx=find(Time1<600&Time1>-400);
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
%                     data4(icell,:)=data4(icell,:)/PeakVal(icell);
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
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:)) ;
                
                subplot(3,3,iarea*3);
                for icell=1:size(data5,1)
                    data5(icell,:)=conv(data5(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data5(icell,:)=data5(icell,:)-baseOn(icell);
%                     data5(icell,:)=data5(icell,:)/PeakVal(icell);
                end
                mm=mean(data5,1);
                ee=std(data5,0,1)/sqrt(size(data5,1));
                idx=find(Time3>=-500&Time3<=170);
                TimePlt=Time3(idx);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', 1.2); hold on;
                
                for icell=1:size(data6,1)
                    data6(icell,:)=conv(data6(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data6(icell,:)=data6(icell,:)-baseOff(icell);
%                     data6(icell,:)=data6(icell,:)/PeakVal(icell);
                end
                mm=mean(data6,1);
                ee=std(data6,0,1)/sqrt(size(data6,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', 1.2); hold on;
                title([areaName{iarea} '-saccade resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data5(:, idx(jj)), data6(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*0.001,'.k','linewidth', 1); hold on;
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:)) ;
            end
        end                
    end
end

%% converge the data from two monkeys together, and plot response of different groups for each area
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; 
Time2=-1000:500;
Time3=-800:200; % RewardOnResp
areaName={'LIP','PUL','FEF'};
yplotRange=[-2 8;-2 12; -4 8];

% yplotRange2=[-2 7];
% yplotRange3=[-2 7];
Nsmooth=50;
Selection={'_','_SU','_cueResponsive_SU'};
for iarea=1
    if iarea==1
        clr=[0 0.5 0; 0 0 0];
    elseif iarea==2
        clr=[1 0 0; 0 0 0];  % [0 0.5 0]  LIP; [1 0 0]  PUL; FEF [0 0.33 0.83]
    else
        clr=[0 0.33 0.83; 0 0 0];
    end
    for isexo=0:1
        figure;
        NN=zeros(1,3); 
        for igroup=1:3
            groupName=Group{igroup};
            if igroup==3
                iselect=2;
            else
                iselect=3;
            end
            iselect1=3;
            Source=cell(2,7);
            for imonkey=1:2
                monkey=monkeyset{imonkey};
                Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
                Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
                Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
                Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
                Source{imonkey,5}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOnResp' Selection{iselect1} '.mat']);
                Source{imonkey,6}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOffResp' Selection{iselect1} '.mat']);
                Source{imonkey,7}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
            end            
            
            MainlineWidth=1.5;                                               
            iCndPlot=1;  % 1 is cue in RF          
            
            SpkType=2;  % 1 is for narrow and 2 is for broad
            idxSpkType1=Source{1,7}.SpikeShape.CellType{iarea}==SpkType;  %&Source{1,7}.SpikeShape.ValueRight{iarea}>=0.3;
            idxSpkType2=Source{2,7}.SpikeShape.CellType{iarea}==SpkType;  %&Source{2,7}.SpikeShape.ValueRight{iarea}>=0.3;
            
            data1=[Source{1,1}.CueOnResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,1}.CueOnResp{iarea,iCndPlot}(idxSpkType2,:)];
            data2=[Source{1,2}.CueOffResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,2}.CueOffResp{iarea,iCndPlot}(idxSpkType2,:)];
            data3=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,3}.ArrayOnResp{iarea,iCndPlot}(idxSpkType2,:)];
            data4=[Source{1,4}.ArrayOffResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,4}.ArrayOffResp{iarea,iCndPlot}(idxSpkType2,:)];
            data5=[Source{1,5}.SaccadeOnResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,5}.SaccadeOnResp{iarea,iCndPlot}(idxSpkType2,:)];
            data6=[Source{1,6}.SaccadeOffResp{iarea,iCndPlot}(idxSpkType1,:); Source{2,6}.SaccadeOffResp{iarea,iCndPlot}(idxSpkType2,:)];
            
            NN(igroup)=size(data1,1);
            if ~isempty(data1)
                subplot(3,3,igroup*3-2);                
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
%                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                end
                
                mm=mean(data1,1);
                ee=std(data1,0,1)/sqrt(size(data1,1));
                idx=find(Time<=500&Time>=-200);
                TimePlt=Time(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                shade_vx=[250 500 500 250 250];
                shade_vy=[-3 -3 8 8 -3];
                patch(shade_vx, shade_vy,[1 1 1]*0.7, 'edgecolor',[1 1 1]*0.7, 'facealpha',0.5); hold on;
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', MainlineWidth); hold on;
                
                baseOff=zeros(size(data2,1),1);
                for icell=1:size(data2,1)
                    data2(icell,:)=conv(data2(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    baseOff(icell)=mean(data2(icell,Time<0&Time>-300));
                    data2(icell,:)=data2(icell,:)-baseOff(icell);
%                     data2(icell,:)=data2(icell,:)/PeakVal(icell);
                end
                mm=mean(data2,1);
                ee=std(data2,0,1)/sqrt(size(data2,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', MainlineWidth); hold on;
%                 title([areaName{iarea} '-cue resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data1(:, idx(jj)), data2(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-1,'.k','linewidth', 1); hold on;
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:)) ;
                yticks(-4:4:12);
                xticks(-400:400:600);
                set(gca,'linewidth',2,'fontsize',10);
                
                
                subplot(3,3,igroup*3-1);
                for icell=1:size(data3,1)
                    data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data3(icell,:)=data3(icell,:)-baseOn(icell);
%                     data3(icell,:)=data3(icell,:)/PeakVal(icell);
                end
                mm=mean(data3,1);
                ee=std(data3,0,1)/sqrt(size(data3,1));
                idx=find(Time1<=600&Time1>=-300);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                TimePlt=Time1(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                
                patch(shade_vx, shade_vy,[1 1 1]*0.7, 'edgecolor',[1 1 1]*0.7, 'facealpha',0.5); hold on;
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', MainlineWidth); hold on;
                
                for icell=1:size(data4,1)
                    data4(icell,:)=conv(data4(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data4(icell,:)=data4(icell,:)-baseOff(icell);
%                     data4(icell,:)=data4(icell,:)/PeakVal(icell);
                end
                mm=mean(data4,1);
                ee=std(data4,0,1)/sqrt(size(data4,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', MainlineWidth); hold on;
%                 title([areaName{iarea} '-array resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data3(:, idx(jj)), data4(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-1,'.k','linewidth', 1); hold on;
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:));
                yticks([-4:4:12]);
                xticks(-400:400:600);
                set(gca,'linewidth',2,'fontsize',10);
                
                subplot(3,3,igroup*3);
                for icell=1:size(data5,1)
                    data5(icell,:)=conv(data5(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data5(icell,:)=data5(icell,:)-baseOn(icell);
%                     data5(icell,:)=data5(icell,:)/PeakVal(icell);
                end
                mm=mean(data5,1);
                ee=std(data5,0,1)/sqrt(size(data5,1));
                idx=find(Time3>=-400&Time3<=170);
                TimePlt=Time3(idx);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];                
                patch(Vx,Vy,clr(1,:),'edgecolor',clr(1,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(1,:),'linewidth', MainlineWidth); hold on;
                
                for icell=1:size(data6,1)
                    data6(icell,:)=conv(data6(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                    data6(icell,:)=data6(icell,:)-baseOff(icell);
%                     data6(icell,:)=data6(icell,:)/PeakVal(icell);
                end
                mm=mean(data6,1);
                ee=std(data6,0,1)/sqrt(size(data6,1));
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr(2,:),'edgecolor',clr(2,:),'facealpha',0.5,'edgealpha',0.5); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr(2,:),'linewidth', MainlineWidth); hold on;
%                 title([areaName{iarea} '-saccade resp']);
                sigDiff=zeros(numel(idx),1);
                for jj=1:numel(idx)
                    sigDiff(jj)=ttest(data5(:, idx(jj)), data6(:,idx(jj)),'tail','right');
                end
                plot(TimePlt(sigDiff==1),  sigDiff(sigDiff==1)*-1,'.k','linewidth', 1); hold on;
                line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
                axis tight;
                ylim(yplotRange(igroup,:));
                yticks([-4:4:12]);
                xticks(-400:400:600);
                set(gca,'linewidth',2,'fontsize',10);
            end
        end                
    end
end

%% plot PSTH of different area together to compare the visual response latency
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; %ArrayOnResp
areaName={'LIP','PUL','FEF'};
Nsmooth=25;

xx=-30:30;
sigma=5;
yy=exp(-xx.^2/sigma^2/2);
kernal=yy/sum(yy);

isexo = 1;
Amplify=1000;
Selection={'_','_SU','_cueResponsive_SU'};
iselect=3; 
MainlineWidth=1.5;           
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
figure;
RespLat=zeros(2,2,3);
igroup=4;
border = 0.15;  % for spike waveform amplitude
iCndPlot=1;  % 1 is cue in RF
single_lat=cell(2,3);
param=cell(2,3);
for SpkType=1:2 % 1 is for narrow and 2 is for broad
    if SpkType==1
        yplotRange=[-2 12;-2 17; -4 8;-2 16];
    else
        yplotRange=[-2 12;-2 17; -4 8;-2 16]/2;
    end    
    groupName=Group{igroup};
    Source=cell(2,5);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);    
%         Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);        
        Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);                
%         Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);                
        Source{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end
        
    for iarea=1:3
%         figure;
        
        SpikeType = [Source{1,5}.SpikeShape.CellType{iarea}, Source{2,5}.SpikeShape.CellType{iarea}];
        WaveformClass = [Source{1,5}.SpikeShape.WaveformClass{iarea}, Source{2,5}.SpikeShape.WaveformClass{iarea}];
        
        idx0=SpikeType == SpkType;        
        cueon=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];  
        arrayon=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}; Source{2,3}.ArrayOnResp{iarea,iCndPlot}];      
                
        data1=cueon(idx0,:);
        data3=arrayon(idx0,:);
        nfig=1;
        if ~isempty(data1)    
%             RespLat(1,SpkType,iarea)=LatencySTD(Time, data1, 0);
%             RespLat(2,SpkType,iarea)=LatencySTD(Time1, data3, 0);
            PeakVal=zeros(size(data1,1),1);
            baseOn=zeros(size(data1,1),1);
            idx=find(Time<=200&Time>=-100);   
            for icell=1:size(data1,1)                
                baseOn(icell)=mean(data1(icell,Time<0&Time>-200));
                data1(icell,:)=data1(icell,:)-baseOn(icell);                                
                data1(icell,:)=conv(data1(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                [param{SpkType,iarea}(icell,:),single_lat{SpkType,iarea}(1,icell)]=LatencyGauss(Time(idx)/1000,data1(icell,idx),0.33,0,clr3(iarea,:));
%                 if param{SpkType,iarea}(icell,6)>0.7 && single_lat{SpkType,iarea}(1,icell)<0.03
%                     subplot(4,4,nfig);
%                     [pp,ll] = LatencyGauss(Time(idx)/1000,data1(icell,idx),0.33,1,clr3(iarea,:));   
%                     nfig=nfig+1;                    
%                     title(['L=' num2str(ll*1000) ',F=' num2str(pp(6),2)]);
%                 end                
            end
%             subplot(2,2,SpkType*2-1);            
%             mm=mean(data1,1)*Amplify;            
%             ee=std(data1,0,1)/sqrt(size(data1,1))*Amplify;                                             
% %             [~,RespLat(1,SpkType,iarea)]=LatencyGauss(Time(idx)/1000,mm(idx),0.33,1,clr3(iarea,:));           
%             TimePlt=Time(idx);
% %             idxPlt=1:10:length(TimePlt);
%             patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
%             line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
%             line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
%             axis tight;
%             ylim(yplotRange(igroup,:)) ;
%             xlim([-50 200])
%             set(gca,'linewidth',2,'fontsize',10);            
                        
%             idx=find(Time1<=200&Time1>=-100);
%             for icell=1:size(data3,1)
%                 data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
%                 data3(icell,:)=data3(icell,:)-baseOn(icell);                
%                 [~,single_lat{SpkType,iarea}(2,icell)]=LatencyGauss(Time1(idx)/1000,data3(icell,idx),0.33,0,clr3(iarea,:));              
%             end
%             subplot(2,2,SpkType*2);                        
%             mm=mean(data3,1)*Amplify;            
%             ee=std(data3,0,1)/sqrt(size(data3,1))*Amplify;            
% %             [~,RespLat(2,SpkType,iarea)]=LatencyGauss(Time1(ildx)/1000,mm(idx),0.33,1,clr3(iarea,:));                        
%             TimePlt=Time1(idx);
% %             idxPlt=1:10:length(TimePlt);            
%             pp(iarea)=patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
%             line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
%             line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
%             axis tight;
%             ylim(yplotRange(igroup,:));
%             xlim([-100 200])
%             set(gca,'linewidth',2,'fontsize',10);                
        end        
    end
%     legend(pp,areaName)
end
save([folder 'Visual_latency_fitting_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat'],'single_lat','-v7.3');
save([folder 'Visual_latency_fitting_param_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat'],'param','-v7.3');

%% plot the distribution of latency of individual cells
close all;
figure;
Nsmooth=50;
load([folder 'Visual_latency_fitting_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat']);
load([folder 'Visual_latency_fitting_param_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat']);
Ncell=zeros(2,3);
for itype=1:2
    subplot(1,2,itype);
    label=cell(1,3);
    for iarea=1:3   
        idx = param{itype,iarea}(:,6)>0.7;
        Ncell(itype,iarea)=sum(idx);
        lat = single_lat{itype,iarea}(1,idx)*1000;
        [count,center]=hist(lat,30);
        cum = cumsum(count)/numel(lat);
        if itype==2
            plot(center,cum,'color',clr3(iarea,:),'linewidth',2,'linestyle','--'); hold on;
        else
            plot(center,cum,'color',clr3(iarea,:),'linewidth',2); hold on;
        end
        label{iarea}=[areaName{iarea} ':' num2str(Ncell(itype,iarea))];
        xlabel('Latency (ms)');
        ylabel('Proportion of units');                
    end
    legend(label,'location','southeast');
    xlim([-20 200]);
    ylim([0 1.1]);
end
%%%%%% converge both narrow and broad together
figure;
Nall=zeros(1,3);
for iarea=1:3
    fitness = [param{1,iarea}(:,6); param{2,iarea}(:,6)];
    lat_total=[single_lat{1,iarea}(1,:), single_lat{2,iarea}(1,:)];
    idx = fitness>0.7;
    Nall(iarea)=sum(idx);
    lat = lat_total(idx)*1000;
    [count,center]=hist(lat,40);
    cum = cumsum(count)/numel(lat);
    plot(center,cum,'color',clr3(iarea,:),'linewidth',2); hold on;
    xlabel('Latency (ms)');
    ylabel('Proportion of units');    
    label{iarea}=[areaName{iarea} ':' num2str(Nall(iarea))];
    xlim([-20 200]);
    ylim([0 1.1]);
end
legend(label,'location','southeast');

%% plot PSTH of cue-on and cue-away together to get the attention modulation latency
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; %ArrayOnResp
Time3=-800:200; % RewardOnResp
areaName={'LIP','PUL','FEF'};

xx=-30:30;
sigma=5;
yy=exp(-xx.^2/sigma^2/2);
kernal=yy/sum(yy);
% Nsmooth=30;
% kernal=ones(1,Nsmooth)/Nsmooth;
isexo = 1;
Amplify=1000;
Selection={'_','_SU','_cueResponsive_SU'};
iselect=3;          
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
figure;
RespLat=zeros(2,2,3);
igroup=4;
border = 0.15;  % for spike waveform amplitude
iCndPlot=1;  % 1 is cue in RF
single_lat=cell(2,3);
param=cell(2,3);

for SpkType=1:2 % 1 is for narrow and 2 is for broad
    if SpkType==1
        yplotRange=[-4 12;-2 17; -4 8;-2 16];
    else
        yplotRange=[-4 12;-2 17; -4 8;-2 16]/2;
    end    
    groupName=Group{igroup};
    Source=cell(2,5);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);    
        Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);        
        Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);                
        Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);                
        Source{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end
        
    for iarea=1:3
        AmpRight = [Source{1,5}.SpikeShape.ValueRight{iarea}, Source{2,5}.SpikeShape.ValueRight{iarea}];        
        AmpLeft = [Source{1,5}.SpikeShape.ValueLeft{iarea}, Source{2,5}.SpikeShape.ValueLeft{iarea}];  
        SpikeType = [Source{1,5}.SpikeShape.CellType{iarea}, Source{2,5}.SpikeShape.CellType{iarea}];
        SpikeWave = [Source{1,5}.SpikeShape.WaveShape{iarea}; Source{2,5}.SpikeShape.WaveShape{iarea}];        
        Valley = abs(min(SpikeWave, [],2));
        NormWave = SpikeWave./repmat(Valley, 1, size(SpikeWave,2));
        NormLeft=AmpLeft./Valley';
        NormRight=AmpRight./Valley';
        idx0=SpikeType == SpkType;
        idx_cell{1} = NormLeft<=border&NormRight<=border&idx0;
        idx_cell{2} = NormLeft<=border&NormRight>border&idx0;
        idx_cell{3} = NormLeft>border&NormRight<=border&idx0;
        idx_cell{4} = NormLeft>border&NormRight>border&idx0; 
        WaveformClass=zeros(1,length(idx0));
        for ii =1:4
            WaveformClass(idx_cell{ii})=ii;
        end
        cueon=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];  
        cueoff=[Source{1,2}.CueOffResp{iarea,iCndPlot}; Source{2,2}.CueOffResp{iarea,iCndPlot}];  
        arrayon=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}; Source{2,3}.ArrayOnResp{iarea,iCndPlot}];     
        arrayoff=[Source{1,4}.ArrayOffResp{iarea,iCndPlot}; Source{2,4}.ArrayOffResp{iarea,iCndPlot}]; 
                
        idxValid=idx_cell{2};  %%% select the interested units 
        data1=cueon(idxValid,:);
        data2=cueoff(idxValid,:);
        data3=arrayon(idxValid,:);
        data4=arrayoff(idxValid,:);
        nfig=1;
        if ~isempty(data1)    
%             RespLat(1,SpkType,iarea)=LatencySTD(Time, data1, 0);
%             RespLat(2,SpkType,iarea)=LatencySTD(Time1, data3, 0);            
            baseOn=zeros(size(data1,1),1);
            baseOff=zeros(size(data2,1),1);
            
            for icell=1:size(data1,1)                
                baseOn(icell)=mean(data1(icell,Time<0&Time>-200));
                baseOff(icell)=mean(data2(icell,Time<0&Time>-200));
                data1(icell,:)=data1(icell,:)-baseOn(icell);                                
                data1(icell,:)=conv(data1(icell,:),kernal,'same');
                data2(icell,:)=data2(icell,:)-baseOff(icell);  
                data2(icell,:)=conv(data2(icell,:),kernal,'same');
                data3(icell,:)=data3(icell,:)-baseOn(icell);                                
                data3(icell,:)=conv(data3(icell,:),kernal,'same');
                data4(icell,:)=data4(icell,:)-baseOff(icell);  
                data4(icell,:)=conv(data4(icell,:),kernal,'same');
            end
            subplot(2,3,iarea+(SpkType-1)*3); 
            plot_cue=0;
            if plot_cue==1
                idx=find(Time<=500&Time>=-200);   
                TimePlt=Time(idx);            
                yy1=data1(:,idx)*Amplify;
                yy2=data2(:,idx)*Amplify;
            else
                idx=find(Time1<=500&Time1>=-300);   
                TimePlt=Time1(idx);            
                yy1=data3(:,idx)*Amplify;
                yy2=data4(:,idx)*Amplify;
            end  
            sig=zeros(1,size(yy1,2));
            for itime=1:size(yy1,2)
               sig(itime)=ttest(yy1(:,itime),yy2(:,itime)) ;                
            end
            mm=mean(yy1,1);
            ee=std(yy1,0,1)/sqrt(size(yy1,1));
            mm2=mean(yy2,1);
            ee2=std(yy2,0,1)/sqrt(size(yy2,1));            
%             [~,RespLat(1,SpkType,iarea)]=LatencyGauss(Time(idx)/1000,mm(idx),0.33,1,clr3(iarea,:));   
            step=1;
            idxT=1:step:length(TimePlt); % downsample timepoint
            xx=TimePlt(idxT);
            patchplot(xx,mm(idxT),ee(idxT),clr3(iarea,:));            
            patchplot(xx,mm2(idxT),ee2(idxT),[1 1 1]*0.4);
            sig_step =  sig(idxT);
            plot(xx(sig_step==1),sig_step(sig_step==1)*-2,'k.'); hold on;  
                      
%             for icell =1:size(yy1,1)
%                 diff=yy1(icell,:)-yy2(icell,:);
%                 plot(TimePlt,diff); hold on;
%                 
%             end
            line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
            line([-300 600],[0 0],'color','k','linestyle','--'); hold on;
            axis tight;
%             ylim(yplotRange(igroup,:)) ;
            xlim([min(TimePlt) max(TimePlt)]);
            set(gca,'linewidth',2,'fontsize',10);                        
        end        
    end    
end

%% plot PSTH of different cell group (spike waveform) in each area
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};   
Time=-600:1200;  %-300:1200;
% Time1=-905:905;
Time1=-505:1305; %ArrayOnResp
areaName={'LIP','PUL','FEF'};

Nsmooth=50;
isexo = 0;
Selection={'_','_SU','_cueResponsive_SU'};
iselect=3; iselect1=3;
MainlineWidth=2;           
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
Amplify=1000;
figure;
RespLat=zeros(2,2,3);
igroup=4;
clr4 = colormap(hsv(4));
iCndPlot=1;  % 1 is cue in RF
for SpkType=1  % 1 is for narrow and 2 is for broad
    if SpkType==1
        yplotRange=[-2 12;-2 17; -4 8;-2 10];
    else
        yplotRange=[-2 12;-2 17; -4 8;-2 16]/2;
    end    
    groupName=Group{igroup};
    Source=cell(2,3);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);        
        Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);                
        Source{imonkey,3}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end
    border = 0.15;
    for iarea=1:3
        AmpRight = [Source{1,3}.SpikeShape.ValueRight{iarea}, Source{2,3}.SpikeShape.ValueRight{iarea}];
        PeakRight = [Source{1,3}.SpikeShape.PeakRight{iarea}, Source{2,3}.SpikeShape.PeakRight{iarea}];
        AmpLeft = [Source{1,3}.SpikeShape.ValueLeft{iarea}, Source{2,3}.SpikeShape.ValueLeft{iarea}];       
        SpikeWave = [Source{1,3}.SpikeShape.WaveShape{iarea}; Source{2,3}.SpikeShape.WaveShape{iarea}];
        Valley = abs(min(SpikeWave, [],2));
        NormWave = SpikeWave./repmat(Valley, 1, size(SpikeWave,2));
        NormLeft=AmpLeft./Valley';
        NormRight=AmpRight./Valley';
        idx0=PeakRight<=25;
        idx_cell{1} = NormLeft<=border&NormRight<=border&idx0;
        idx_cell{2} = NormLeft<=border&NormRight>border&idx0;
        idx_cell{3} = NormLeft>border&NormRight<=border&idx0;
        idx_cell{4} = NormLeft>border&NormRight>border&idx0;                
        
        data1=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];        
        data2=[Source{1,2}.ArrayOnResp{iarea,iCndPlot}; Source{2,2}.ArrayOnResp{iarea,iCndPlot}];                
                
        if ~isempty(data1)                        
            baseOn=zeros(size(data1,1),1);
            for icell=1:size(data1,1)
                data1(icell,:)=conv(data1(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                baseOn(icell)=mean(data1(icell,Time<0&Time>-300));
                data1(icell,:)=data1(icell,:)-baseOn(icell);                               
            end
            subplot(3,2,iarea*2-1); 
            for itype=1:4
                mm=mean(data1(idx_cell{itype},:),1);
                ee=std(data1(idx_cell{itype},:),0,1)/sqrt(sum(idx_cell{itype}));
                idx=find(Time<=200&Time>=-100);                
                TimePlt=Time(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr4(itype,:),'edgecolor',clr4(itype,:),'facealpha',0.3,'edgealpha',0.3); hold on;
                plot(TimePlt, mm(idx)*Amplify,'color',clr4(itype,:),'linewidth', MainlineWidth); hold on;
            end
            line([30 30],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
            line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
            axis tight;
            ylim(yplotRange(igroup,:)) ;
            xlim([-50 200])
            set(gca,'linewidth',2,'fontsize',10);            
            
            subplot(3,2,iarea*2);            
            for icell=1:size(data2,1)
                data2(icell,:)=conv(data2(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                data2(icell,:)=data2(icell,:)-baseOn(icell);
            end
            for itype=1:4
                mm=mean(data2(idx_cell{itype},:),1);
                ee=std(data2(idx_cell{itype},:),0,1)/sqrt(sum(idx_cell{itype}));
                idx=find(Time1<=200&Time1>=-100);
                yy1=(mm(idx)-ee(idx))*Amplify;
                yy2=(mm(idx)+ee(idx))*Amplify;
                TimePlt=Time1(idx);
                Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
                Vy=[yy1 fliplr(yy2) yy1(1)];
                patch(Vx,Vy,clr4(itype,:),'edgecolor',clr4(itype,:),'facealpha',0.3,'edgealpha',0.3); hold on;
                pp(iarea)=plot(TimePlt, mm(idx)*Amplify,'color',clr4(itype,:),'linewidth', MainlineWidth); hold on;
            end
            line([30 30],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
            line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
            axis tight;
            ylim(yplotRange(igroup,:));
            xlim([-100 200])
            set(gca,'linewidth',2,'fontsize',10);                
        end        
    end    
end


%% plot the responses of individual units
close all;
Time=-600:1200; 
Time1=-505:1305; 
Time2=-800:200;
Nsmooth=80;
clr=[1 0 0; 0 0 1];
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
igroup=2;
groupName=Group{igroup};
areaName={'LIP','PUL','FEF'};
folder='Z:\RujiaChen\Results\';
for imonkey=1
    monkey=monkeyset{imonkey};
    for isexo=1
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ElecInfo_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOnResp_cueResponsive_SU.mat']);
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SaccadeOffResp_cueResponsive_SU.mat']);
        iCndPlot=1;
        for iarea=3          
            figure;           
            if ~isempty(CueOnResp{iarea,iCndPlot})  
                baseOn=zeros(1,size(CueOnResp{iarea,iCndPlot},1));
                baseOff=zeros(1,size(CueOffResp{iarea,iCndPlot},1));
                num=0;
                for icell=26  %1:size(CueOnResp{iarea,iCndPlot},1)
                    if ElecInfo{iarea,1}(icell,3)>0  %==2||ElecInfo{iarea,1}(icell,3)==4
                        data1=conv(CueOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOn(icell)=mean(data1(Time<0&Time>-200));
                        data1=data1-baseOn(icell);
                        PeakCueOn=max(data1(Time>50&Time<200));
                        %                     if PeakCueOn>0.005
                        %                     PeakVal(icell)=max(ArrayOnResp{iarea}(icell,Time1>0&Time1<200))-baseOn(icell);
                        %                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
                        
                        data2=conv(CueOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        baseOff(icell)=mean(data2(Time<0&Time>-200));
                        data2=data2-baseOff(icell);
                        
                        num=num+1;
                        CellID(num)=icell;
%                         subplot(8,8,num);
                        subplot(1,2,1);
                        idx=find(Time<=400&Time>=-200);
                        TimePlt=Time(idx);
                        plot(TimePlt, data1(idx)*1000,'color','k','linewidth', 2); hold on;
%                         plot(TimePlt, data2(idx)*1000,'color',clr(2,:),'linewidth', 1.5); hold on;
                        Vx=[0 200 200 0 0];                       
                        Vy=[-1 -1 8 8 -1];
                        patch(Vx,Vy,[0.7 0.7 0.7],'edgecolor',[1 1 1],'facealpha',0.5,'edgealpha',0); hold on;
                        axis tight;
                        xticks([-200 0 200 400]);
                        ylim([-1 8]);
                        set(gca,'box','off','linewidth',2.0,'fontsize',20);
%                         title([num2str(ElecInfo{iarea}(icell,1)) '-' num2str(ElecInfo{iarea}(icell,2)) '-' num2str(ElecInfo{iarea}(icell,3))]);
                    end
                end
                
%                 num=0;
%                 figure;
%                 for icell=1:size(ArrayOnResp{iarea},1)
%                     if ElecInfo{iarea,1}(icell,3)>0
%                         data1=conv(ArrayOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
%                         data1=data1-baseOn(icell);
%                         PeakCueOn=max(data1(Time>50&Time<200));
%                         %                     if PeakCueOn>0.005
%                         %                     PeakVal(icell)=max(ArrayOnResp{iarea}(icell,Time1>0&Time1<200))-baseOn(icell);
%                         %                     data1(icell,:)=data1(icell,:)/PeakVal(icell);
%                         
%                         data2=conv(ArrayOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
%                         data2=data2-baseOff(icell);
%                         
%                         num=num+1;
%                         subplot(8,8,num);
%                         idx=find(Time1<=500&Time1>=-400);
%                         TimePlt=Time1(idx);
%                         plot(TimePlt, data1(idx)*1000,'color',clr(1,:),'linewidth', 1.2); hold on;
%                         plot(TimePlt, data2(idx)*1000,'color',clr(2,:),'linewidth', 1.2); hold on;
%                         title([num2str(ElecInfo{iarea}(icell,1)) '-' num2str(ElecInfo{iarea}(icell,2)) '-' num2str(ElecInfo{iarea}(icell,3))]);
%                     end
%                 end
%                 
%                 num=0;
%                 figure;
                for icell=26  %1:size(SaccadeOnResp{iarea},1)
                    if ElecInfo{iarea,1}(icell,3)>0
                        data1=conv(SaccadeOnResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data1=data1-baseOn(icell);
                        
                        data2=conv(SaccadeOffResp{iarea,iCndPlot}(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                        data2=data2-baseOff(icell);
                        
                        num=num+1;
%                         subplot(8,8,num);
                        subplot(1,2,2);
                        idx=find(Time2<=150&Time2>=-300);
                        TimePlt=Time2(idx);
                        plot(TimePlt, data1(idx)*1000,'color','k','linewidth', 2); hold on;
%                         plot(TimePlt, data2(idx)*1000,'color',clr(2,:),'linewidth', 1.5); hold on;
                        Vx=[-100 100 100 -100 -100];
                        Vy=[-1 -1 8 8 -1];
                        patch(Vx,Vy,[0.7 0.7 0.7],'edgecolor',[1 1 1],'facealpha',0.5,'edgealpha',0); hold on;
                        axis tight;
                        ylim([-1 8]);
                        xticks([-200 0 200 400]);
                        set(gca,'box','off','linewidth',2.0,'fontsize',20);
%                         title([num2str(ElecInfo{iarea}(icell,1)) '-' num2str(ElecInfo{iarea}(icell,2)) '-' num2str(ElecInfo{iarea}(icell,3))]);
                    end
                end
            end                        
        end              
    end
end

%% get the d-prime for different groups 
clear all; 
close all;
folder='Z:\Rujia\Results\';
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
load('Z:\Rujia\Results\Mikey_RecordingInfo.mat');
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    
Time=-600:1200; 
Time1=-505:1305;
window=[0.05 0.025];
Fs=1000;
for igroup=4 %1:3
    groupName=Group{igroup};
    for imonkey=1 
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end        
        for isexo=0:1
            DprimeCue=cell(3,1);
            DprimeArray=cell(3,1);            
%             ElecInfo=cell(3,2);
            DprimeCueDelay=cell(3,1);
            DprimeArrayDelay=cell(3,1);
            pValueCueDelay=cell(3,1);
            pValueArrayDelay=cell(3,1);
            for iarea=1:3
                num=0;
                for ifile = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{ifile};
                    load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                    load([folder 'CorrectTrialParam_' date '.mat']);
                    load([folder 'bVisResp_' date '_SingleUnit.mat']);
                    load([folder 'bMotorResp_' date '_SingleUnit.mat']);
                    load([folder 'bVisMotor_' date '_SingleUnit.mat']);                    
                    load([folder 'CueOn_' date '_SingleUnit.mat']);
                    load([folder 'RFOn_' date '_SingleUnit.mat']);
                    load([folder 'SaccadeOn_' date '_SingleUnit.mat']);                    
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
                            
                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            idxCueT=Time>=-300&Time<=1000;
                            idxArrayT=Time1>=-300&Time1<=1000;
                            CueOnResp=squeeze(SpikeTrain.cueAlign(idx1,ich,idxCueT)); 
                            ArrayOnResp=squeeze(SpikeTrain.arrayAlign(idx1,ich,idxArrayT));  
                            BaseOn=repmat(squeeze(mean(SpikeTrain.cueAlign(idx1,ich,Time>=-300&Time<=0),3)),1,size(CueOnResp,2));
                            CueOnResp=CueOnResp-BaseOn;
                            ArrayOnResp=ArrayOnResp-BaseOn;
                             
                            idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            CueOffResp=squeeze(SpikeTrain.cueAlign(idx2,ich,idxCueT));
                            ArrayOffResp=squeeze(SpikeTrain.arrayAlign(idx2,ich,idxArrayT));
                            BaseOff=repmat(squeeze(mean(SpikeTrain.cueAlign(idx2,ich,Time>=-300&Time<=0),3)),1,size(CueOffResp,2));
                            CueOffResp=CueOffResp-BaseOff;
                            ArrayOffResp=ArrayOffResp-BaseOff;
                            
%                             [DprimeCue{iarea}(num,:),~,~]=getDPrime(CueOnResp,CueOffResp,window,Fs,0);
%                             [DprimeArray{iarea}(num,:),~,T]=getDPrime(ArrayOnResp,ArrayOffResp,window,Fs,0);

%                             ElecInfo{iarea,1}(num,:)=[SpikeTrain.channel(ich), ifile, RFOn(ich)];
                            
                            [DprimeCueDelay{iarea}(num),pValueCueDelay{iarea}(num)]=getDPrime(mean(CueOnResp(:,550:800),2),mean(CueOffResp(:,550:800),2));
                            [DprimeArrayDelay{iarea}(num),pValueArrayDelay{iarea}(num)]=getDPrime(mean(ArrayOnResp(:,550:800),2),mean(ArrayOffResp(:,550:800),2));
                        end
                    end
                end
            end
            
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCue_cueResponsive_SU.mat'],'DprimeCue','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArray_cueResponsive_SU.mat'],'DprimeArray','-v7.3'); 
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeT.mat'],'T','-v7.3');
            
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat'],'DprimeCueDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat'],'DprimeArrayDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat'],'pValueCueDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat'],'pValueArrayDelay','-v7.3');
        end
    end
end

%% plot the averaged dprime during cue delay and array delay, compare dprime for narrow-spike units and broad-spike units
close all;
% clear all;
folder='Z:\Rujia\Results\';
AreaName = {'LIP', 'PUL', 'FEF'};
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'}; 
SpkRange=[5 10; 12 25];
clr4 = colormap(hsv(4));
border = 0.15;
for igroup=4
    groupName=Group{igroup};
%     figure;
    for isexo=0:1
        figure;
        clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];  % [0 0.5 0]  LIP; [1 0 0]  PUL
        SourceData=cell(2,4);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            SourceData{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat']);  %
            SourceData{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
        end
        mDpCue=zeros(2,3);
        mDpArray=zeros(2,3);
        eeDpCue=zeros(2,3);
        eeDpArray=zeros(2,3);
        Ncell=zeros(2,3);
        CellNum=zeros(1, 3);
        for iarea=1:3
            data1=[SourceData{1,1}.DprimeCueDelay{iarea}, SourceData{2,1}.DprimeCueDelay{iarea}];
            data2=[SourceData{1,2}.DprimeArrayDelay{iarea}, SourceData{2,2}.DprimeArrayDelay{iarea}];
            pvalue1=[SourceData{1,3}.pValueCueDelay{iarea}, SourceData{2,3}.pValueCueDelay{iarea}];
            pvalue2=[SourceData{1,4}.pValueArrayDelay{iarea}, SourceData{2,4}.pValueArrayDelay{iarea}];
              
            AmpRight = [SourceData{1,5}.SpikeShape.ValueRight{iarea}, SourceData{2,5}.SpikeShape.ValueRight{iarea}];
            PeakRight = [SourceData{1,5}.SpikeShape.PeakRight{iarea}, SourceData{2,5}.SpikeShape.PeakRight{iarea}];
            AmpLeft = [SourceData{1,5}.SpikeShape.ValueLeft{iarea}, SourceData{2,5}.SpikeShape.ValueLeft{iarea}];            
            SpikeWave = [SourceData{1,5}.SpikeShape.WaveShape{iarea}; SourceData{2,5}.SpikeShape.WaveShape{iarea}];
            Valley = abs(min(SpikeWave, [],2));            
            NormLeft=AmpLeft./Valley';
            NormRight=AmpRight./Valley';
            
            idx0 = ~isnan(data1)&~isnan(data2);
            idx{1} = NormLeft<=border&NormRight<=border&idx0;
            idx{2} = NormLeft<=border&NormRight>border&idx0;
            idx{3} = NormLeft>border&NormRight<=border&idx0;
            idx{4} = NormLeft>border&NormRight>border&idx0;
            subplot(3,2,2*iarea-1);
            for kk=4
                scatter(PeakRight(idx{kk})/40*1000, data1(idx{kk}), 10 ,'markerfacecolor',clr4(kk,:),'markeredgecolor',clr4(kk,:)); hold on;
            end
            line([300 300],[-1 2], 'color', 'k', 'linestyle','--');
%             line([250 250],[-1 1.5], 'color', 'k', 'linestyle','--');
            line([100 1000],[0 0], 'color', 'k', 'linestyle','--');
            title([AreaName{iarea} '-Cue']);
            xlim([100 650])
            ylim([-1 1.5])
            
            subplot(3,2,2*iarea);
            for kk=4
                scatter(PeakRight(idx{kk})/40*1000, data2(idx{kk}),10,'markerfacecolor',clr4(kk,:),'markeredgecolor',clr4(kk,:)); hold on;
            end
%             line([250 250],[-1 1.5], 'color', 'k', 'linestyle','--');
            line([300 300],[-1 2], 'color', 'k', 'linestyle','--');
            line([100 1000],[0 0], 'color', 'k', 'linestyle','--');
            title([AreaName{iarea} '-Array']);
            xlim([100 650])
            ylim([-1 1.5])
            
%             iSig = pvalue1<1|pvalue2<1; % select units show significant d'~=0 
%             idxCell{1}= iSig & idx0& PeakRight>= SpkRange(1,1) & PeakRight<= SpkRange(1,2) & idx4;
%             idxCell{2}= iSig & idx0& PeakRight>= SpkRange(2,1) & PeakRight<= SpkRange(2,2) & idx4;                       
%             for icelltype=1:2
%                 ss=[0 0];                            
%                 Ncell(icelltype,iarea)=sum(idxCell{icelltype});                 
%                 mDpCue(icelltype,iarea)=mean(data1(idxCell{icelltype}));
%                 eeDpCue(icelltype,iarea)=std(data1(idxCell{icelltype}))/sqrt(Ncell(icelltype,iarea));               
%                 ss(1)=ttest(data1(idxCell{icelltype}));
% %                 idx1=idx0&pvalue2<0.05;
%                 mDpArray(icelltype,iarea)=mean(data2(idxCell{icelltype}));
%                 eeDpArray(icelltype,iarea)=std(data2(idxCell{icelltype}))/sqrt(Ncell(icelltype,iarea));                               
%                 ss(2)=ttest(data2(idxCell{icelltype}));
%                 subplot(3,2,iarea*2-1);
%                 if icelltype==1  
%                     Barlegend(1)=bar(icelltype+3*isexo, mDpCue(icelltype,iarea),'facecolor','w','edgecolor',clr3(iarea,:),'linewidth',2);hold on;                    
%                 else
%                     Barlegend(2)=bar(icelltype+3*isexo, mDpCue(icelltype,iarea),'facecolor',clr3(iarea,:),'edgecolor',clr3(iarea,:),'linewidth',2);hold on;
%                 end
%                 if ss(1)==1
%                     plot(icelltype+3*isexo,-0.05, 'k*'); hold on;
%                 end
%                 errorbar(icelltype+isexo*3, mDpCue(icelltype,iarea),eeDpCue(icelltype,iarea), 'color', 'k','linewidth',1.5); hold on;
%                 set(gca,'box','off','linewidth',2,'fontsize',10);
%                 xticks([1.5 4.5]);   
%                 yticks(-0.4:0.4:0.8);
%                 xticklabels({'Endo','Exo'});
%                 ylim([-0.1 0.4]);
% 
%                 subplot(3,2,iarea*2);
%                 if icelltype==1  %isexo==0
%                     bar(icelltype+3*isexo, mDpArray(icelltype,iarea),'facecolor','w','edgecolor',clr3(iarea,:),'linewidth',2); hold on;
%                 else
%                     bb(iarea)=bar(icelltype+3*isexo, mDpArray(icelltype,iarea),'facecolor',clr3(iarea,:),'edgecolor',clr3(iarea,:),'linewidth',2); hold on;
%                 end
%                 if ss(2)==1
%                     plot(icelltype+3*isexo,-0.05, 'k*'); hold on;
%                 end
%                 errorbar(icelltype+isexo*3, mDpArray(icelltype,iarea),eeDpArray(icelltype,iarea), 'color', 'k','linewidth',1.5); hold on;               
%                 set(gca,'box','off','linewidth',2,'fontsize',10);
%                 xticks([1.5 4.5]);
%                 yticks(-0.4:0.4:0.8);
%                 xticklabels({'Endo','Exo'});
%                 ylim([-0.1 0.4]);
%             end
%             legend(Barlegend,{'NS','BS'},'location','northwest');
        end
    end
end

%% plot the averaged dprime during cue delay and array delay, compare dprime for narrow-spike units and broad-spike units
%%%%%% balance the cell number by resampling from the raw dateset 
% close all;
% clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'}; 
SpkRange=[5 10; 12 25];
border = 0.15;
for igroup=2
    
    groupName=Group{igroup};
    figure;
    for isexo=0:1
        clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];  % [0 0.5 0]  LIP; [1 0 0]  PUL
        SourceData=cell(2,4);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            SourceData{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat']);  %
            SourceData{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
        end
        mDpCue=zeros(2,3);
        mDpArray=zeros(2,3);
        eeDpCue=zeros(2,3);
        eeDpArray=zeros(2,3);
        Ncell=zeros(2,3);
        CellNum=zeros(1, 3);
        for iarea=1:3
            data1=[SourceData{1,1}.DprimeCueDelay{iarea}, SourceData{2,1}.DprimeCueDelay{iarea}];
            data2=[SourceData{1,2}.DprimeArrayDelay{iarea}, SourceData{2,2}.DprimeArrayDelay{iarea}];
            pvalue1=[SourceData{1,3}.pValueCueDelay{iarea}, SourceData{2,3}.pValueCueDelay{iarea}];
            pvalue2=[SourceData{1,4}.pValueArrayDelay{iarea}, SourceData{2,4}.pValueArrayDelay{iarea}];          
              
            AmpRight = [SourceData{1,5}.SpikeShape.ValueRight{iarea}, SourceData{2,5}.SpikeShape.ValueRight{iarea}];
            PeakRight = [SourceData{1,5}.SpikeShape.PeakRight{iarea}, SourceData{2,5}.SpikeShape.PeakRight{iarea}];
            AmpLeft = [SourceData{1,5}.SpikeShape.ValueLeft{iarea}, SourceData{2,5}.SpikeShape.ValueLeft{iarea}];            
            SpikeWave = [SourceData{1,5}.SpikeShape.WaveShape{iarea}; SourceData{2,5}.SpikeShape.WaveShape{iarea}];
            Valley = abs(min(SpikeWave, [],2));            
            NormLeft=AmpLeft./Valley';
            NormRight=AmpRight./Valley';
                   
            idx2 = NormLeft<=border&NormRight>border;
    
            iSig = pvalue1<0.05|pvalue2<0.05; % select units show significant d'~=0 
            idxCell{1}= iSig & ~isnan(data1)&~isnan(data2)&PeakRight>= SpkRange(1,1) & PeakRight<= SpkRange(1,2);  % & idx2;
            idxCell{2}= iSig & ~isnan(data1)&~isnan(data2)&PeakRight>= SpkRange(2,1) & PeakRight<= SpkRange(2,2) ;% & idx2;
            CellNum(iarea) = min([sum(idxCell{1}), sum(idxCell{2})]);
            
            for icelltype=1:2
                ss=[0 0];                            
                Ncell(icelltype,iarea)=sum(idxCell{icelltype}); 
                idxRaw = find(idxCell{icelltype}==1);
                NN=1000;
                mDpCue_new=zeros(1, NN);
                mDpArray_new = zeros(1, NN);
                for itime = 1:NN
                    idxNew=randi(numel(idxRaw),1,CellNum(iarea));
                    idx1=idxRaw(idxNew); 
                    mDpCue_new(itime) = mean(data1(idx1));
                    mDpArray_new(itime)=mean(data2(idx1)) ;                               
                end
                mDpCue(icelltype,iarea)=mean(mDpCue_new);
                eeDpCue(icelltype,iarea)=std(mDpCue_new);  %/sqrt(100);
                ss(1)=sum(mDpCue_new>0)/NN<0.05 | sum(mDpCue_new>0)/NN>0.95;
                
%                 idx1=idx0&pvalue2<0.05;
                mDpArray(icelltype,iarea)=mean(mDpArray_new);
                eeDpArray(icelltype,iarea)=std(mDpArray_new);  %/sqrt(100);               
                ss(2)=sum(mDpArray_new>0)/NN<0.05 | sum(mDpArray_new>0)/NN>0.95;
                
                subplot(3,2,iarea*2-1);
                if icelltype==1  %isexo==0
                    Barlegend(1)=bar(icelltype+3*isexo, mDpCue(icelltype,iarea),'facecolor','w','edgecolor',clr3(iarea,:),'linewidth',2);hold on;                    
                else
                    Barlegend(2)=bar(icelltype+3*isexo, mDpCue(icelltype,iarea),'facecolor',clr3(iarea,:),'edgecolor',clr3(iarea,:),'linewidth',2);hold on;
                end
                if ss(1)==1
                    plot(icelltype+3*isexo,-0.05, 'k*'); hold on;
                end
                errorbar(icelltype+isexo*3, mDpCue(icelltype,iarea),eeDpCue(icelltype,iarea), 'color', 'k','linewidth',1.5); hold on;
                set(gca,'box','off','linewidth',2,'fontsize',10);
                xticks([1.5 4.5]);   
                yticks(-0.4:0.4:0.8);
%                 yticks(-0.4:0.2:0.8);
                xticklabels({'Endo','Exo'});
                ylim([-0.2 0.4]);

                subplot(3,2,iarea*2);
                if icelltype==1  %isexo==0
                    bar(icelltype+3*isexo, mDpArray(icelltype,iarea),'facecolor','w','edgecolor',clr3(iarea,:),'linewidth',2); hold on;
                else
                    bb(iarea)=bar(icelltype+3*isexo, mDpArray(icelltype,iarea),'facecolor',clr3(iarea,:),'edgecolor',clr3(iarea,:),'linewidth',2); hold on;
                end
                if ss(2)==1
                    plot(icelltype+3*isexo,-0.05, 'k*'); hold on;
                end
                errorbar(icelltype+isexo*3, mDpArray(icelltype,iarea),eeDpArray(icelltype,iarea), 'color', 'k','linewidth',1.5); hold on;               
                set(gca,'box','off','linewidth',2,'fontsize',10);
                xticks([1.5 4.5]);
                yticks(-0.4:0.4:0.8);
%                 yticks(-0.4:0.2:0.8);
                xticklabels({'Endo','Exo'});
                ylim([-0.2 0.4]);
            end
            legend(Barlegend,{'NS','BS'},'location','northwest');
        end
    end
%     ll=legend(bb,{'LIP','PUL','FEF'});
%     set(ll,'box','off');
end

%% plot the distribution of d-prime on different conditions after balance the cell number
close all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'}; 
SpkRange=[5 10; 12 25];
igroup=4;
groupName=Group{igroup};
AreaName = {'LIP', 'PUL', 'FEF'};
figure;
lt= {'-', '--'};
Ncell=cell(1,2);
for isexo=0:1
    clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];  % [0 0.5 0]  LIP; [1 0 0]  PUL
    SourceData=cell(2,4);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        SourceData{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat']);  %
        SourceData{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat']);
        SourceData{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat']);
        SourceData{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat']);
        SourceData{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end    
    Ncell{isexo+1}=zeros(2,3);
    CellNum=zeros(1, 3);
    for iarea=1:3
        data1=[SourceData{1,1}.DprimeCueDelay{iarea}, SourceData{2,1}.DprimeCueDelay{iarea}];
        data2=[SourceData{1,2}.DprimeArrayDelay{iarea}, SourceData{2,2}.DprimeArrayDelay{iarea}];
        pvalue1=[SourceData{1,3}.pValueCueDelay{iarea}, SourceData{2,3}.pValueCueDelay{iarea}];
        pvalue2=[SourceData{1,4}.pValueArrayDelay{iarea}, SourceData{2,4}.pValueArrayDelay{iarea}];
        WaveWidth=[SourceData{1,5}.SpikeShape.PeakRight{iarea}, SourceData{2,5}.SpikeShape.PeakRight{iarea}];
        WaveAmp=[SourceData{1,5}.SpikeShape.ValueRight{iarea}, SourceData{2,5}.SpikeShape.ValueRight{iarea}];
             
        iSig = pvalue1<0.2|pvalue2<0.2;
        idxCell1= iSig &~isnan(data1)&~isnan(data2)&WaveWidth>= SpkRange(1,1) & WaveWidth<= SpkRange(1,2);  
        idxCell2= iSig &~isnan(data1)&~isnan(data2)&WaveWidth>= SpkRange(2,1) & WaveWidth<= SpkRange(2,2);   %iSig &
        CellNum(iarea) = min([sum(idxCell1), sum(idxCell2)]);
        
        subplot(3, 2, iarea*2-1+isexo);
        clr2 = [1 1 1; clr3(iarea,:)];
        for icelltype=1:2
            ss=[0 0];
            idx0=~isnan(data1)&~isnan(data2)& WaveWidth>= SpkRange(icelltype,1) & WaveWidth<= SpkRange(icelltype,2)  ;   %& WaveAmp>=0.3            
            idx2=find(idx0 & iSig);  %
%             idx2= find(idx0==1);
            Ncell{isexo+1}(icelltype,iarea)=numel(idx2);
        
            %%%%%% plot the histogram of raw dprime
%             [h, edges] = histcounts(data1(idx2), 10);   % for cue delay
%             [h, edges] = histcounts(data2(idx2), 10);    % for array delay
%             BinCenter = (edges(1:end-1)+edges(2:end))/2;
%             plot(BinCenter, h, 'color', clr3(iarea, :), 'linestyle', lt{icelltype}, 'linewidth', 2.0); hold on;
            
%%%%%%% resample from the raw dprime values and plot the distribution of
%%%%%%% averaged surrogate d-prime 
            mDpCue_new=zeros(1, 200);
            mDpArray_new = zeros(1, 200);
            for itime = 1:200
                idxNew=randi(numel(idx2),1,CellNum(iarea));
                idx1=idx2(idxNew);
                mDpCue_new(itime) = mean(data1(idx1));
                mDpArray_new(itime)=mean(data2(idx1)) ;
            end            
            [h, edges] = histcounts(mDpCue_new, 10);   % for cue delay
%             [h, edges] = histcounts(mDpArray_new, 10);    % for array delay
            BinCenter = (edges(1:end-1)+edges(2:end))/2;
            plot(BinCenter, h, 'color', clr3(iarea, :), 'linestyle', lt{icelltype}, 'linewidth', 2.0); hold on;
        end
        line([0 0], [0 max(h)+5], 'color', 'k'); hold on;
        legend('NS', 'WS', 'location', 'best');
        set(gca, 'box', 'off');
        xlim([-0.1 0.3+isexo*0.2])
        title([AreaName{iarea} '-' cuetype{isexo+1}]);
    end
end

%%  plot the temporal dprime for narrow- and broad- spikes units
close all;
folder='Z:\Rujia\Results\';
monkeyset= {'Mikey','Vasco'}; 
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
igroup=2;
groupName=Group{igroup};
SpkRange=[5 10; 12 25];
SpkType=1;
cuetype={'endo','exo'};
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];  % [0 0.5 0]  LIP; [1 0 0]  PUL  ;0.83 0.33 0 FEF
figure;
for isexo=0:1
    DataSource=cell(2,5);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        DataSource{imonkey,1} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCue_cueResponsive_SU.mat']);
        DataSource{imonkey,2} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArray_cueResponsive_SU.mat'])   ;     
        DataSource{imonkey,3} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_cueResponsive_SU.mat']);
        DataSource{imonkey,4} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_cueResponsive_SU.mat']);
        DataSource{imonkey,5} = load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end 
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeT.mat']);
%     load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_cueResponsive_SU.mat']);
%     load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_cueResponsive_SU.mat']);
        
    subplot(2,2,isexo*2+1);    
    for iarea=1:3
        data1=[DataSource{1,1}.DprimeCue{iarea}; DataSource{2,1}.DprimeCue{iarea}];
        pvalue1=[DataSource{1,3}.pValueCueDelay{iarea}, DataSource{2,3}.pValueCueDelay{iarea}];      
        WaveWidth=[DataSource{1,5}.SpikeShape.PeakRight{iarea}, DataSource{2,5}.SpikeShape.PeakRight{iarea}];
        WaveAmp=[DataSource{1,5}.SpikeShape.ValueRight{iarea}, DataSource{2,5}.SpikeShape.ValueRight{iarea}];
        idxCell = WaveWidth>= SpkRange(SpkType,1) & WaveWidth<= SpkRange(SpkType,2);
        idx=~isnan(mean(data1,2)) & pvalue1'<0.05; %& idxCell'
        Ncell(1,iarea)=sum(idx);
        mm=mean(data1(idx,:),1);
        ee=std(data1(idx,:),0,1)/sqrt(sum(idx));
        TimePlt=T-0.3;
        pp(iarea)=patchplot(TimePlt, mm, ee, clr3(iarea,:)); hold on;      
        line([min(TimePlt) max(TimePlt)],[0 0],'color','k','linestyle','--'); hold on;
        line([0 0],[-0.1 0.6],'color','k','linestyle','-.'); hold on;
    end
%     ll=legend(pp,{'LIP','PUL','FEF'});
%     set(ll,'box','off');
    title(['Cue-' cuetype{isexo+1}]);
    xlim([-0.2 0.5]);
%     ylim([-0.2 0.7]);
    set(gca,'box','off','linewidth', 2.0,'fontsize',10);
    
    subplot(2,2,isexo*2+2);
    for iarea=1:3
        data1=[DataSource{1,2}.DprimeArray{iarea}; DataSource{2,2}.DprimeArray{iarea}];        
        pvalue2=[DataSource{1,4}.pValueArrayDelay{iarea}, DataSource{2,4}.pValueArrayDelay{iarea}];
        WaveWidth=[DataSource{1,5}.SpikeShape.PeakRight{iarea}, DataSource{2,5}.SpikeShape.PeakRight{iarea}];
        WaveAmp=[DataSource{1,5}.SpikeShape.ValueRight{iarea}, DataSource{2,5}.SpikeShape.ValueRight{iarea}];
        idxCell = WaveWidth>= SpkRange(SpkType,1) & WaveWidth<= SpkRange(SpkType,2);        
        idx=~isnan(mean(data1,2))& pvalue2'<0.05;    %& idxCell'
        Ncell(2,iarea)=sum(idx);
        mm=mean(data1(idx,:),1);
        ee=std(data1(idx,:),0,1)/sqrt(sum(idx));
        TimePlt=T-0.3;
        patchplot(TimePlt, mm, ee, clr3(iarea,:)); hold on;      
        line([min(TimePlt) max(TimePlt)],[0 0],'color','k','linestyle','--'); hold on;
        line([0 0],[-0.1 0.3],'color','k','linestyle','-.'); hold on;
    end    
    title(['Array-' cuetype{isexo+1}]);
    xlim([-0.2 0.5]);
%     ylim([-0.1 0.3]);
    set(gca,'box','off','linewidth', 2.0,'fontsize',10);
end
%% plot the cell ratio, averaged d-prime and distribution of d-prime for different groups
close all;
clear all;
folder='Z:\Rujia\Results\';
cuetype={'endo','exo'};
monkeyset={'Mikey','Vasco'};
AreaName={'LIP','PUL','FEF'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'}; 
for igroup=4    
    groupName=Group{igroup};
    figure;
    clr3=[0 0.5 0;0.5 0 0; 0 0.33 0.84];  % [0 0.5 0]  LIP; [1 0 0]  PUL
    for isexo=0:1       
        SourceData=cell(2,4);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            SourceData{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat']);
        end
              
        for iarea=1:3
            data1{iarea}(isexo+1,:)=[SourceData{1,1}.DprimeCueDelay{iarea}, SourceData{2,1}.DprimeCueDelay{iarea}];
            data2{iarea}(isexo+1,:)=[SourceData{1,2}.DprimeArrayDelay{iarea}, SourceData{2,2}.DprimeArrayDelay{iarea}];
            pvalue1{iarea}(isexo+1,:)=[SourceData{1,3}.pValueCueDelay{iarea}, SourceData{2,3}.pValueCueDelay{iarea}];
            pvalue2{iarea}(isexo+1,:)=[SourceData{1,4}.pValueArrayDelay{iarea}, SourceData{2,4}.pValueArrayDelay{iarea}];            
        end
    end    
    %%%%%%%%%%%%%%% plot the distribution of units showing significant
    %%%%%%%%%%%%%%% attention modulation d'>0 (p<0.05)
      
    for isexo=0:1
        ratioDp=zeros(3,4);
        for iarea=1:3
            idx0=~isnan(data1{iarea}(isexo+1,:))&~isnan(data2{iarea}(isexo+1,:));
            idx1=idx0&pvalue1{iarea}(isexo+1,:)<0.05;    % cue delay
            idx2=idx0&pvalue2{iarea}(isexo+1,:)<0.05;    %  array delay
            idx3=idx1&idx2;
            idx4=idx0&idx1==0&idx2==0;
            ratioDp(iarea,1)=sum(idx1)/sum(idx0)-sum(idx3)/sum(idx0); % cue only
            ratioDp(iarea,2)=sum(idx3)/sum(idx0);  % both
            ratioDp(iarea,3)=sum(idx2)/sum(idx0)-sum(idx3)/sum(idx0); %array only
            ratioDp(iarea,4)=sum(idx4)/sum(idx0); % none
            subplot(2,3,iarea+isexo*3);
            pp=pie(ratioDp(iarea,:),{'Cue','Both','Array','None'}); hold on;
            set(pp(2:2:8),'fontsize',12)
            set(pp(7),'facecolor',[1 1 1]*0.8);
            set(gca,'box','off','linewidth',2.0,'fontsize',12);
        end
    end
%%%%%%%% plot the mean d-prime averaged across all task-related units
    figure;
    Ncell=zeros(2,3);
    ss=cell(1,2);
    for isexo=0:1
        mDprime=cell(1,2);
        eeDprime=cell(1,2);
        for iarea=1:3
            idx0=~isnan(data1{iarea}(isexo+1,:)) & ~isnan(data2{iarea}(isexo+1,:));
            idx1=pvalue1{iarea}(isexo+1,:)<0.05 | pvalue2{iarea}(isexo+1,:)<0.05; %either d' during cue delay or array delay is significant
            idx2=idx0&idx1;
            Ncell(isexo+1,iarea)=sum(idx2);
            
            mDprime{1}(iarea)=mean(data1{iarea}(isexo+1,idx2));
            ss{isexo+1}(1,iarea)=ttest(data1{iarea}(isexo+1,idx2));
            eeDprime{1}(iarea)=std(data1{iarea}(isexo+1,idx2))/sqrt(sum(idx2));
            
            mDprime{2}(iarea)=mean(data2{iarea}(isexo+1,idx2));
            ss{isexo+1}(2,iarea)=ttest(data2{iarea}(isexo+1,idx2));
            eeDprime{2}(iarea)=std(data2{iarea}(isexo+1,idx2))/sqrt(sum(idx2));
                  
            for idelay=1:2
                subplot(2,1,idelay);
                if isexo==0
                    bar(iarea+isexo*4, mDprime{idelay}(iarea),'facecolor','w','edgecolor',clr3(iarea,:),'linewidth',2);hold on;
                else
                    bar(iarea+isexo*4, mDprime{idelay}(iarea),'facecolor',clr3(iarea,:),'edgecolor',clr3(iarea,:),'linewidth',2);hold on;
                end
                errorbar(iarea+isexo*4, mDprime{idelay}(iarea),eeDprime{idelay}(iarea), 'color', 'k','linewidth',1.5); hold on;
                if ss{isexo+1}(idelay,iarea)==1
                    plot(iarea+isexo*4,ss{isexo+1}(idelay,iarea)*0.45,'k*');
                end
                set(gca,'box','off','linewidth',2,'fontsize',12);
                xticks([2 6]);
                xticklabels({'Endo','Exo'});
                ylim([-0.1 0.5]);
            end          
        end
    end
     
      %%%%%% plot the d' distribution for endo and exo
      figure;
      for isexo=0:1
          for iarea=1:3
              idx0=~isnan(data1{iarea}(isexo+1,:)) & ~isnan(data2{iarea}(isexo+1,:));
              idx1=pvalue1{iarea}(isexo+1,:)<0.05 | pvalue2{iarea}(isexo+1,:)<0.05; %either d' during cue delay or array delay is significant
              idx2=idx0&idx1;
              subplot(2,3,iarea+3*isexo);
              data=transpose([data1{iarea}(isexo+1,idx2); data2{iarea}(isexo+1,idx2)]);
              hist(data);
              line([0 0],[0 55],'color','k','linestyle','--','linewidth',1.5); hold on;
              if iarea==3
                  ll=legend('Cue','Array');
                  set(ll,'box','off');
              end
              set(gca,'box','off','linewidth',2.0,'fontsize',12);
              xlim([-1.5 1.5]);              
          end
      end
end

%%  plot the temporal dprime
close all;
folder='Z:\Rujia\Results\';
monkeyset= {'Mikey','Vasco'}; 
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
igroup=4;
groupName=Group{igroup};
cuetype={'endo','exo'};
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];  % [0 0.5 0]  LIP; [1 0 0]  PUL  ;0.83 0.33 0 FEF
figure;
for isexo=0:1
    DataSource=cell(2,4);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        DataSource{imonkey,1} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCue_cueResponsive_SU.mat']);
        DataSource{imonkey,2} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArray_cueResponsive_SU.mat'])   ;     
        DataSource{imonkey,3} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_cueResponsive_SU.mat']);
        DataSource{imonkey,4} = load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_cueResponsive_SU.mat']);
    end 
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeT.mat']);
%     load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_cueResponsive_SU.mat']);
%     load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_cueResponsive_SU.mat']);
        
    subplot(2,2,isexo+1);    
    for iarea=1:3
        data1=[DataSource{1,1}.DprimeCue{iarea}; DataSource{2,1}.DprimeCue{iarea}];
        pvalue1=[DataSource{1,3}.pValueCueDelay{iarea}, DataSource{2,3}.pValueCueDelay{iarea}];
        idx=~isnan(mean(data1,2))&pvalue1'<0.05;        
        mm=mean(data1(idx,:),1);
        ee=std(data1(idx,:),0,1)/sqrt(sum(idx));
        TimePlt=T-0.3;
        pp(iarea)=patchplot(TimePlt, mm, ee, clr3(iarea,:)); hold on;      
        line([min(TimePlt) max(TimePlt)],[0 0],'color','k','linestyle','--'); hold on;
        line([0 0],[-0.1 0.6],'color','k','linestyle','-.'); hold on;
    end
%     ll=legend(pp,{'LIP','PUL','FEF'});
%     set(ll,'box','off');
    title(['Cue-' cuetype{isexo+1}]);
    xlim([-0.2 0.5]);
    ylim([-0.1 0.6]);
    set(gca,'box','off','linewidth', 2.0,'fontsize',20);
%     ylim([-0.2 0.5]);
    
    subplot(2,2,isexo+3);
    for iarea=1:3
        data1=[DataSource{1,2}.DprimeArray{iarea}; DataSource{2,2}.DprimeArray{iarea}];
        pvalue1=[DataSource{1,4}.pValueArrayDelay{iarea}, DataSource{2,4}.pValueArrayDelay{iarea}];
        idx=~isnan(mean(data1,2))&pvalue1'<0.05;        
        mm=mean(data1(idx,:),1);
        ee=std(data1(idx,:),0,1)/sqrt(sum(idx));
        TimePlt=T-0.3;
        patchplot(TimePlt, mm, ee, clr3(iarea,:)); hold on;      
        line([min(TimePlt) max(TimePlt)],[0 0],'color','k','linestyle','--'); hold on;
        line([0 0],[-0.1 0.3],'color','k','linestyle','-.'); hold on;
    end    
    title(['Array-' cuetype{isexo+1}]);
    xlim([-0.2 0.5]);
    ylim([-0.1 0.3]);
    set(gca,'box','off','linewidth', 2.0,'fontsize',20);
%     ylim([-0.1 0.2])
end

%% get the attention modulation index for different groups 
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
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    
Time=-600:1200; 
Time1=-505:1305;
window=[0.05 0.025];
Fs=1000;
SmoothWin=50;
for igroup=[2 4] %1:3
    groupName=Group{igroup};
    for imonkey=1:2  
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end        
        for isexo=0:1
            AMICue=cell(3,1);
            AMIArray=cell(3,1);            
%             ElecInfo=cell(3,2);
%             DprimeCueDelay=cell(3,1);
%             DprimeArrayDelay=cell(3,1);
%             pValueCueDelay=cell(3,1);
%             pValueArrayDelay=cell(3,1);
            for iarea=1:3
                num=0;
                for ifile = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{ifile};
                    load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                    load([folder 'CorrectTrialParam_' date '.mat']);
                    load([folder 'bVisResp_' date '_SingleUnit.mat']);
                    load([folder 'bMotorResp_' date '_SingleUnit.mat']);
                    load([folder 'bVisMotor_' date '_SingleUnit.mat']);                    
                    load([folder 'CueOn_' date '_SingleUnit.mat']);
                    load([folder 'RFOn_' date '_SingleUnit.mat']);
                    load([folder 'SaccadeOn_' date '_SingleUnit.mat']);                    
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
                            %%%%%% prime position
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
                            idxCueT=Time>=-300&Time<=1000;
                            idxArrayT=Time1>=-300&Time1<=1000;
                            CueOnResp=mean(squeeze(SpikeTrain.cueAlign(idx1,ich,idxCueT)),1); 
                            ArrayOnResp=mean(squeeze(SpikeTrain.arrayAlign(idx1,ich,idxArrayT)),1);  
                            BaseOn=squeeze(mean(mean(SpikeTrain.cueAlign(idx1,ich,Time>=-300&Time<=0),3),1));
%                             CueOnResp=CueOnResp-BaseOn;                            
%                             ArrayOnResp=ArrayOnResp-BaseOn;
                            CueOnResp=conv(CueOnResp, ones(1,SmoothWin)/SmoothWin,'same');
                            ArrayOnResp=conv(ArrayOnResp, ones(1,SmoothWin)/SmoothWin,'same');

                            idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            CueOffResp=mean(squeeze(SpikeTrain.cueAlign(idx2,ich,idxCueT)),1);
                            ArrayOffResp=mean(squeeze(SpikeTrain.arrayAlign(idx2,ich,idxArrayT)),1);
                            BaseOff=squeeze(mean(mean(SpikeTrain.cueAlign(idx2,ich,Time>=-300&Time<=0),3),1));
%                             CueOffResp=CueOffResp-BaseOff;
%                             ArrayOffResp=ArrayOffResp-BaseOff;
                            CueOffResp=conv(CueOffResp, ones(1,SmoothWin)/SmoothWin,'same');
                            ArrayOffResp=conv(ArrayOffResp, ones(1,SmoothWin)/SmoothWin,'same');
                            
                            AMICue{iarea}(num,:)=(CueOnResp-CueOffResp)./(CueOnResp+CueOffResp);
                            AMIArray{iarea}(num,:)=(ArrayOnResp-ArrayOffResp)./(ArrayOnResp+ArrayOffResp);
%                             ElecInfo{iarea,1}(num,:)=[SpikeTrain.channel(ich), ifile, RFOn(ich)];
                            
%                             [DprimeCueDelay{iarea}(num),pValueCueDelay{iarea}(num)]=getDprime(mean(CueOnResp(:,500:800),2),mean(CueOffResp(:,500:800),2));
%                             [DprimeArrayDelay{iarea}(num),pValueArrayDelay{iarea}(num)]=getDprime(mean(ArrayOnResp(:,500:800),2),mean(ArrayOffResp(:,500:800),2));
                        end
                    end
                end
            end
            
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_AMICue_cueResponsive_SU.mat'],'AMICue','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_AMIArray_cueResponsive_SU.mat'],'AMIArray','-v7.3'); 
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeT.mat'],'T','-v7.3');
%             
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_cueResponsive_SU.mat'],'DprimeCueDelay','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_cueResponsive_SU.mat'],'DprimeArrayDelay','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_cueResponsive_SU.mat'],'pValueCueDelay','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_cueResponsive_SU.mat'],'pValueArrayDelay','-v7.3');
        end
    end
end

%% plot the time course of AMI (attention modulation index)
close all;
Time=-300:1000;
figure;
clr3=eye(3);
clr2=[1 0 0 ; 0 0 1];
Group={'VisResp','VisMotor','MotorResp','AllVisual'};   
igroup=4;
groupName=Group{igroup};
for isexo=0:1
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_AMICue_cueResponsive_SU.mat']);
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_AMIArray_cueResponsive_SU.mat']);
    load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    for iarea=1:3
        mCue=zeros(1,length(Time));
        mArray=zeros(1,length(Time));
        for itime=1:size(AMICue{iarea},2)
            idx=~isnan(AMICue{iarea}(:,itime));
            mCue(itime)=mean(AMICue{iarea}(idx,itime),1);
            idx=~isnan(AMIArray{iarea}(:,itime));
            mArray(itime)=mean(AMIArray{iarea}(idx,itime),1);
        end
                
        subplot(3,2,iarea*2-1);
        plot(Time,mCue,'color',clr2(isexo+1,:)); hold on;
        line([-300 1000],[0 0],'linestyle','--');
        xlim([-200 600]);
        %     ylim([-1 1]);
        subplot(3,2,iarea*2);
        plot(Time,mArray,'color',clr2(isexo+1,:)); hold on;
        line([-300 1000],[0 0],'linestyle','--');
        xlim([-300 600]);
        %     ylim([-1 1]);
    end
end






