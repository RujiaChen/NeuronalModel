%% plot PSTH of different area together to compare the visual response latency
close all;
clear all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    %'ArrayResp';  %'CueOn';   %'CueRespUnit';  %'cuedelay_facilitation';
Time=-600:1200;  %-300:1200;
Time1=-505:1305; %ArrayOnResp
areaName={'LIP','PUL','FEF'};

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
groupName=Group{igroup};
border = 0.15;  % for spike waveform amplitude
Nsmooth=50;
iCndPlot=1;  % 1 is cue in RF
for SpkType=1:2 % 1 is for narrow and 2 is for broad
    if SpkType==1
        yplotRange=[-2 12;-2 17; -4 8;-2 16];
    else
        yplotRange=[-2 12;-2 17; -4 8;-2 16]/2;
    end        
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
        Latency=[Lat{iarea,1}, Lat{iarea,2}];
        fitness = [param{iarea,1}, param{iarea,2}];
        idx=fitness>0.4&Latency<100;
        idx0=SpikeType == SpkType&idx;                
        cueon=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];  
        arrayon=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}; Source{2,3}.ArrayOnResp{iarea,iCndPlot}];      
                
        data1=cueon(idx0,:);
        data3=arrayon(idx0,:);
        nfig=1;
        if ~isempty(data1)    
            PeakVal=zeros(size(data1,1),1);
            baseOn=zeros(size(data1,1),1);
            idx=find(Time<=200&Time>=-100);   
            for icell=1:size(data1,1)                
                baseOn(icell)=mean(data1(icell,Time<0&Time>-200));
                data1(icell,:)=data1(icell,:)-baseOn(icell);                                
                data1(icell,:)=conv(data1(icell,:),ones(1,Nsmooth)/Nsmooth,'same');

            end
            subplot(2,2,SpkType*2-1);            
            mm=mean(data1,1)*Amplify;            
            ee=std(data1,0,1)/sqrt(size(data1,1))*Amplify;                                                      
            TimePlt=Time(idx);
            patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
            line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
            line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
            axis tight;
            ylim(yplotRange(igroup,:)) ;
            xlim([-50 200])
            set(gca,'linewidth',2,'fontsize',10);            
                        
            idx=find(Time1<=200&Time1>=-100);
            for icell=1:size(data3,1)
                data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
                data3(icell,:)=data3(icell,:)-baseOn(icell);
            end
            subplot(2,2,SpkType*2);
            mm=mean(data3,1)*Amplify;
            ee=std(data3,0,1)/sqrt(size(data3,1))*Amplify;            
            TimePlt=Time1(idx);           
            pp(iarea)=patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
            line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
            line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
            axis tight;
            ylim(yplotRange(igroup,:));
            xlim([-100 200])
            set(gca,'linewidth',2,'fontsize',10);
        end        
    end
    legend(pp,areaName)
end

%% plot the neural responses on cue-on/cue-off conditions for narrow- and broad-spike units, divided into early and late groups based on visual latency
close all;
Nsmooth=50;
load([folder 'Visual_latency_two_monkey_fitting_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat']);
load([folder 'Visual_latency_two_monkey_fitting_param_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat']);
isexo=1; % 0 for endo and 1 for exo
iarea=2;
Ncell=zeros(1,2);
figure;
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
LatRange=[0, 80; 80, 200];
irange=1;
for SpkType=1:2 % 1 is for narrow and 2 is for broad
    if SpkType==1
        yplotRange=[-2 12;-2 20; -4 8;-2 25];
    else
        yplotRange=[-2 12;-1 6; -4 8;-1 6];
    end
    Source=cell(2,5);
    for imonkey=1:2
        monkey=monkeyset{imonkey};
        Source{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOnResp' Selection{iselect} '.mat']);
        Source{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_CueOffResp' Selection{iselect} '.mat']);
        Source{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOnResp' Selection{iselect} '.mat']);
        Source{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_ArrayOffResp' Selection{iselect} '.mat']);
        Source{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
    end
    
    SpikeType = [Source{1,5}.SpikeShape.CellType{iarea}, Source{2,5}.SpikeShape.CellType{iarea}];
    WaveformClass = [Source{1,5}.SpikeShape.WaveformClass{iarea}, Source{2,5}.SpikeShape.WaveformClass{iarea}];
    Latency=[Lat{iarea,1}, Lat{iarea,2}];
    fitness = [param{iarea,1}, param{iarea,2}];    
    idx1=fitness>0.4&Latency>LatRange(irange,1)&Latency<=LatRange(irange,2);
    
    idx0=SpikeType == SpkType&idx1;
    cueon=[Source{1,1}.CueOnResp{iarea,iCndPlot}; Source{2,1}.CueOnResp{iarea,iCndPlot}];
    cueoff=[Source{1,2}.CueOffResp{iarea,iCndPlot}; Source{2,2}.CueOffResp{iarea,iCndPlot}];
    arrayon=[Source{1,3}.ArrayOnResp{iarea,iCndPlot}; Source{2,3}.ArrayOnResp{iarea,iCndPlot}];
    arrayoff=[Source{1,4}.ArrayOffResp{iarea,iCndPlot}; Source{2,4}.ArrayOffResp{iarea,iCndPlot}];
    Ncell(SpkType)=sum(idx0);
    
    data1=cueon(idx0,:);
    data2=cueoff(idx0,:);
    data3=arrayon(idx0,:);
    data4=arrayoff(idx0,:);
    nfig=1;
    if ~isempty(data1)
        PeakVal=zeros(size(data1,1),1);
        baseOn=zeros(size(data1,1),1);
        baseOff=zeros(size(data2,1),1);
        idx=find(Time<=500&Time>=-100);
        for icell=1:size(data1,1)
            baseOn(icell)=mean(data1(icell,Time<0&Time>-200));
            baseOff(icell)=mean(data2(icell,Time<0&Time>-200));
            data1(icell,:)=data1(icell,:)-baseOn(icell);
            data1(icell,:)=conv(data1(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
            data2(icell,:)=data2(icell,:)-baseOff(icell);
            data2(icell,:)=conv(data2(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
        end
        subplot(2,2,SpkType*2-1);
        mm=mean(data1,1)*Amplify;
        ee=std(data1,0,1)/sqrt(size(data1,1))*Amplify;
        TimePlt=Time(idx);
        patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
        
        nn=mean(data2,1)*Amplify;
        ee=std(data2,0,1)/sqrt(size(data2,1))*Amplify;        
        patchplot(TimePlt,nn(idx),ee(idx),[0 0 0]);
        text(300,yplotRange(igroup,2)*3/4,['N=' num2str(Ncell(SpkType))]);
        line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
        line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
        axis tight;
        ylim(yplotRange(igroup,:)) ;
        xlim([-100 500])
        set(gca,'linewidth',2,'fontsize',10);
        
        idx=find(Time1<=500&Time1>=-100);
        for icell=1:size(data3,1)
            data3(icell,:)=conv(data3(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
            data3(icell,:)=data3(icell,:)-baseOn(icell);
            data4(icell,:)=conv(data4(icell,:),ones(1,Nsmooth)/Nsmooth,'same');
            data4(icell,:)=data4(icell,:)-baseOff(icell);
        end
        subplot(2,2,SpkType*2);
        mm=mean(data3,1)*Amplify;
        ee=std(data3,0,1)/sqrt(size(data3,1))*Amplify;
        TimePlt=Time1(idx);
        patchplot(TimePlt,mm(idx),ee(idx),clr3(iarea,:));
        
        nn=mean(data4,1)*Amplify;
        ee=std(data4,0,1)/sqrt(size(data4,1))*Amplify;        
        patchplot(TimePlt,nn(idx),ee(idx),[0 0 0]);
        
        line([0 0],yplotRange(igroup,:),'color','k','linestyle','--'); hold on;
        line([-300 400],[0 0],'color','k','linestyle','--'); hold on;
        axis tight;
        ylim(yplotRange(igroup,:));
        xlim([-100 500])
        set(gca,'linewidth',2,'fontsize',10);    
    end
end

%% plot the averaged dprime during cue delay and array delay, compare dprime for narrow-spike units and broad-spike units
close all;
folder='Z:\Rujia\Results\';
AreaName = {'LIP', 'PUL', 'FEF'};
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'}; 
SpkRange=[5 10; 12 25];
clr4 = colormap(hsv(4));
Nsmooth=50;
load([folder 'Visual_latency_two_monkey_fitting_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat']);
load([folder 'Visual_latency_two_monkey_fitting_param_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat']);
LatRange=[0, 80;80, 200];
DelayName={'Cue','Array'};
for igroup=4
    groupName=Group{igroup};
%     figure;
    for isexo=0
        figure;
        clr3=[0 0.5 0;0.8 0 0; 0 0.33 0.84];  % [0 0.5 0]  LIP; [1 0 0]  PUL
        SourceData=cell(2,4);
        for imonkey=1:2
            monkey=monkeyset{imonkey};
            SourceData{imonkey,1}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat']);  %
            SourceData{imonkey,2}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,3}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,4}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat']);
            SourceData{imonkey,5}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);            
        end
        
        for iarea=2
            data{1}=[SourceData{1,1}.DprimeCueDelay{iarea}, SourceData{2,1}.DprimeCueDelay{iarea}];
            data{2}=[SourceData{1,2}.DprimeArrayDelay{iarea}, SourceData{2,2}.DprimeArrayDelay{iarea}];
            pvalue{1}=[SourceData{1,3}.pValueCueDelay{iarea}, SourceData{2,3}.pValueCueDelay{iarea}];
            pvalue{2}=[SourceData{1,4}.pValueArrayDelay{iarea}, SourceData{2,4}.pValueArrayDelay{iarea}];
            SpikeType=[SourceData{1,5}.SpikeShape.CellType{iarea}, SourceData{2,5}.SpikeShape.CellType{iarea}];
            WaveClass = [SourceData{1,5}.SpikeShape.WaveformClass{iarea},SourceData{2,5}.SpikeShape.WaveformClass{iarea}];
            
            Latency=[Lat{iarea,1}, Lat{iarea,2}];
            fitness = [param{iarea,1}, param{iarea,2}];
        
            for SpkType=1:2           
                idx0 = ~isnan(data{1})&~isnan(data{2})&SpikeType==SpkType&fitness>0.4&Latency>0;                                                
                %%%%% plot the dprime against latency
                for idelay=1:2
                    subplot(2,2,SpkType*2-2+idelay);
                    dp = data{idelay}(idx0);
                    pv=pvalue{idelay}(idx0);
                    ll=Latency(idx0);
                    Ncell(SpkType,iarea)=sum(idx0);
                    scatter(ll,dp, 15 ,'markeredgecolor',clr3(iarea,:)); hold on;%  'markerfacecolor',clr3(iarea,:)
                    idx2=pv<0.1;
                    scatter(ll(idx2),dp(idx2), 15 ,'markeredgecolor',clr3(iarea,:),'markerfacecolor',clr3(iarea,:)); hold on;%  
                    line([0 200],[0 0], 'color', 'k', 'linestyle','--');
                    line([LatRange(1,2) LatRange(1,2)],[-1 1.2], 'color', 'k', 'linestyle','--');
                    ylim([-1 1.2]);
                    xlim([0 200]);
                    title([DelayName{idelay} ' delay, N = ' num2str(Ncell(SpkType,iarea))]);
                end
                
                %%%%%% plot the mean dprime for early and late groups
%                 for idelay=1:2
%                     subplot(2,2,SpkType*2-2+idelay);
%                     for irange=1:2
%                         idx1=idx0&Latency>LatRange(irange,1)&Latency<=LatRange(irange,2);
%                         Ncell(SpkType,irange)=sum(idx1);
%                         dp=data{idelay}(idx1);
%                         mm=mean(dp);
%                         ee=std(dp)/sqrt(length(dp));
%                         bar(irange,mm,0.5); hold on;
%                         errorbar(irange,mm,ee,'linewidth',1.5); hold on;
%                         [ss,pp]=ttest(dp,0,'tail','right');
%                         if ss>0
%                             plot(irange,-0.1,'*'); hold on;
%                         end
%                         text(irange-0.3,0.25,{['N=' num2str(sum(idx1))],['p=' num2str(pp,2)]});
%                     end
%                     title([DelayName{idelay} ' delay']);
%                     xticks([1,2]);
%                     xticklabels({'Early','Late'});
%                     xlim([0 3]);
%                     ylim([-0.15 0.3]);
%                     set(gca,'box','off');
%                 end
            end            
        end
    end
end


%% get the visual response latency of individual units by t-test methods 
xx=-30:30;
sigma=10;
yy=exp(-xx.^2/sigma^2/2);
kernal=yy/sum(yy);
groupName=Group{igroup};
Source=cell(2,1);
S=cell(1,2);
for imonkey=1:2
    monkey=monkeyset{imonkey};
    S{imonkey}=load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeTrain_cueResponsive_SU.mat']);
    Source{imonkey}=load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
end
for SpkType=1:2 % 1 is for narrow and 2 is for broad
    for iarea=1:3
        nfig=1;
        count=0;
        for imonkey=1:2
            idx0 = find(Source{imonkey}.SpikeShape.CellType{iarea} == SpkType);
            idx1 = find(ismember(Source{imonkey}.SpikeShape.WaveformClass{iarea},[2,4]));
            for icell = idx0
                count=count+1;
                data1=S{imonkey}.SpikeTrain_group.CueOn{iarea}{icell};
                for itrl=1:size(data1,1)
                   data1(itrl,:)=conv(data1(itrl,:),ones(1,25)/25,'same') ;                    
                end
                bin=10;
                RespLatency{iarea,SpkType}(count)=LatencySTD(Time, data1,bin, 0);
            end
        end
    end
end

% save([folder 'Visual_latency_ttest_' cuetype{isexo+1} '_' groupName '.mat'],'RespLatency','-v7.3');

%% get the visual response latency by curve fitting method for each session
isexo=1;
Nsmooth=25;
TimeCue=-600:1000;
groupName = 'AllVisual';
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
for imonkey=1:2
    monkey=monkeyset{imonkey};    
    load([folder monkey '_exo_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']);  %  _5sigma    
    RespLatency=cell(size(SpikeLFPInfo.ElecInfo,1),3);
    idx=find(TimeCue<=200&TimeCue>=-100);
    for iarea=1:3        
        for idate=1:size(SpikeLFPInfo.ElecInfo,1)
            if ~isempty(SpikeLFPInfo.ElecInfo{idate, iarea})
                NN=size(SpikeLFPInfo.ElecInfo{idate, iarea},1);                                
                RespLatency{idate,iarea}=zeros(2,NN);
                for icell =1:NN
                    cueResp=mean(SpikeLFPInfo.CueOnSpike{idate,iarea}{icell},1);
                    cueResp=cueResp-mean(cueResp(TimeCue>-200&TimeCue<0));
                    cueResp=conv(cueResp,ones(1,Nsmooth)/Nsmooth,'same');
                    [param,single_lat]=LatencyGauss(TimeCue(idx)/1000,cueResp(idx),0.33,0,clr3(iarea,:));
                    
                    RespLatency{idate,iarea}(1,icell)=single_lat*1000;
                    RespLatency{idate,iarea}(2,icell)=param(6);                     
                end
            end
        end
    end
    save([folder 'Visual_latency_fitting_' num2str(Nsmooth)  'ms_smooth_exo_cue_' monkey '_' groupName '.mat'],'RespLatency','-v7.3');
end

%% plot the cumulative distribution of individual latencies
close all;
monkeyset={'Mikey','Vasco'};
groupName = 'AllVisual';
Nsmooth=50;
lines={'-','--'};
figure;
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
for itype =1:2
    subplot(1,2,itype);
    param=cell(3,2);
    Lat=cell(3,2);
    for iarea=1:3
        Lat_all=[];        
        for imonkey=1:2            
            monkey=monkeyset{imonkey};  %'Mikey';
            load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
            load([folder 'Visual_latency_fitting_' num2str(Nsmooth)  'ms_smooth_exo_cue_' monkey '_' groupName '.mat']);
            %     load([folder 'Visual_latency_ttest_' cuetype{isexo+1} '_' groupName '.mat']);
            Nfile=size(RespLatency,1);
            for ifile = 1:Nfile
                if ~isempty(RespLatency{ifile,iarea})
                    Lat{iarea,imonkey} = [Lat{iarea,imonkey}, RespLatency{ifile,iarea}(1,:)];
                    param{iarea,imonkey} = [param{iarea,imonkey}, RespLatency{ifile,iarea}(2,:)];
                end
            end
            idx = Lat{iarea,imonkey}>0&param{iarea,imonkey}>0.4&SpikeShape.CellType{iarea}==itype;
            Lat_all=[Lat_all, Lat{iarea,imonkey}(idx)];
        end
        cc=0:10:200;
        [count,cc] = hist(Lat_all,30);
        cum = cumsum(count)/numel(Lat_all);
        label{iarea}=[areaName{iarea} ':' num2str(numel(Lat_all))];        
        pp(iarea)=plot(cc, cum, 'linewidth',2.0,'color',clr3(iarea,:),'linestyle',lines{itype}); hold on;
        if iarea==2&&itype==2
            [~,idx]=min(abs(cum-0.5));
            border = cc(idx);
            line([0 200], [0.5, 0.5], 'linestyle','--','color','k'); hold on;
            line([border, border], [0, 1], 'linestyle','--','color','k'); hold on;
        end
    end
    legend(pp,label,'location','southeast');
    xlabel('Latency (ms)');
    ylabel('Proportion of units');
%     save([folder 'Visual_latency_two_monkey_fitting_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat'],'Lat');
%     save([folder 'Visual_latency_two_monkey_fitting_param_' num2str(Nsmooth)  'ms_smooth_cue_' groupName '.mat'],'param');
end

%% compare the visual latency computed by two different methods
close all;
groupName=Group{igroup};
Nsmooth = 25;
load([folder 'Visual_latency_ttest_' cuetype{isexo+1} '_' groupName '.mat']);
load([folder 'Visual_latency_fitting_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat']);
load([folder 'Visual_latency_fitting_param_' num2str(Nsmooth) 'ms_' cuetype{isexo+1} '_' groupName '.mat']);
figure;
Ncell=zeros(3,2);
for iarea=1:3
    for itype=1:2
        subplot(2,3,(itype-1)*3+iarea);
        idx1=RespLatency{iarea,itype}>0&RespLatency{iarea,itype}<400;  
        idx2=transpose(param{itype,iarea}(:,6)>0.4);  %&RespLatency{iarea,itype}<200;
        idx3=idx1&idx2;
        Ncell(iarea, itype) = sum(idx1)+sum(idx2)-sum(idx3);
        scatter(single_lat{itype,iarea}(idx1)*1000,RespLatency{iarea,itype}(idx1),12,'markeredgecolor','r'); hold on;
        scatter(single_lat{itype,iarea}(idx2)*1000,RespLatency{iarea,itype}(idx2),12,'markeredgecolor','b'); hold on;
        scatter(single_lat{itype,iarea}(idx1&idx2)*1000,RespLatency{iarea,itype}(idx1&idx2),12,'markerfacecolor','b'); hold on;
        line([0, 500], [0, 500], 'linestyle', '--','color','k');     
        title(['N=' num2str(Ncell(iarea,itype))]);
        xlim([0 400]);
        ylim([0 400]);
    end
end

%% plot the distribution of latency of individual cells based on t-test method
close all;
figure;
Ncell=zeros(2,3);
clr3=[0 0.5 0;1 0 0;0 0.33 0.83];
load([folder 'Visual_latency_ttest_exo_' groupName '.mat']);
for itype=1:2
    subplot(1,2,itype);
    label=cell(1,3);
    for iarea=1:3   
        idx = RespLatency{iarea,itype}>0;
        Ncell(itype,iarea)=sum(idx);
        lat = RespLatency{iarea,itype}(idx);
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
    xlim([0 200]);
    ylim([0 1.1]);
end
%%%%%%%% converge both narrow and broad together
figure;
Nall=zeros(1,3);
for iarea=1:3   
    lat_total=[RespLatency{iarea,1}, RespLatency{iarea,2}];  %
    idx = lat_total>0;
    Nall(iarea)=sum(idx);
    lat = lat_total(idx);
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




