%% compare the local spike-field coherence between early and late cell groups within pulvinar
close all;
clear all;
% folder='Z:\Rujia\Results\';
folder='C:\Data\';
CellTypeName={'NarrowUnits','BroadUnits'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
iarea=2;
irange=1;
LatRange=[20,66 ;66, 200];
for iCellType=1
    celltype=CellTypeName{iCellType};
    figure;
    for isexo=0:1             
        C1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_despike_linear_SU.mat']);  %
        C2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_despike_linear_SU.mat']);
        V1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_despike_linear_SU.mat']);
        V2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_despike_linear_SU.mat']);                
        Vars=fieldnames(C1.C);                        
        subplot(2,2,isexo*2+1);
        ivar=3; 
        Var=Vars{ivar};
        coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});
        latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
        param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];        
        
        t_cue=C1.C.t_cue*1000;
        idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
        idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
        NN=sum(idx&idx0);
        mm=transpose(squeeze(mean(coh(idx&idx0,:,:),1)));
        t=t_cue;
        baseline=repmat(mean(mm(:,t>-200&t<-50),2),1,size(mm,2));
        
        mm=(mm-baseline)./baseline*100;
        imagesc(t,C1.C.f,mm); hold on;
        line([0,0],[0 100],'color','k','linestyle','--'); hold on;
        line([-200, 500],[15, 15],'color','k','linestyle','--'); hold on;
        line([-200, 500],[30, 30],'color','k','linestyle','--'); hold on;
        axis xy;
        xlim([-200, 500]);
        title([cuetype{isexo+1} ', cue onset' ]);
        caxis([-30,30]);
%         caxis([-0.05 0.05]); % for raw SFC
        colorbar;    
        
        subplot(2,2,isexo*2+2);
        ivar=4;                
        Var=Vars{ivar};
        coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});
        latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
        param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];        
       
        idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
        idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
        NN=sum(idx&idx0);
        mm=transpose(squeeze(mean(coh(idx&idx0,:,:),1)));
        
        t=C1.C.t_array*1000;
        coh1=cat(1,C1.C.(Vars{ivar-1}){iarea},C2.C.(Vars{ivar-1}){iarea});
        latency1=[squeeze(V1.Vislatency.(Vars{ivar-1}){iarea}(:,1)); squeeze(V2.Vislatency.(Vars{ivar-1}){iarea}(:,1))];
        param1=[squeeze(V1.Vislatency.(Vars{ivar-1}){iarea}(:,2)); squeeze(V2.Vislatency.(Vars{ivar-1}){iarea}(:,2))];
        
        idx2=latency1<=LatRange(irange,2)&latency1>LatRange(irange,1);
        idx3=~isnan(squeeze(mean(mean(coh1,2),3)))&param1>0.4;
        nn=transpose(squeeze(mean(coh1(idx2&idx3,:,:),1)));
        baseline = repmat(mean(nn(:,t_cue>-200&t_cue<-50),2),1,size(mm,2));
        mm=(mm-baseline)./baseline*100;
        imagesc(t,C1.C.f,mm); hold on;
        line([0,0],[0 100],'color','k','linestyle','--'); hold on;
        line([-200, 500],[15, 15],'color','k','linestyle','--'); hold on;
        line([-200, 500],[30, 30],'color','k','linestyle','--'); hold on;
        axis xy;
        xlim([-200, 500]);
        title([cuetype{isexo+1} ', array onset' ]);
        caxis([-30,30]);
%         caxis([-0.05 0.05]);  % for raw SFC
        colorbar;
    end
end

%% Compare the difference of SFC during some frequency bands (time periods) between early and late cell groups within pulvinar 
close all;
clear all;
% folder='Z:\Rujia\Results\';
folder='C:\Data\';
CellTypeName={'NarrowUnits','BroadUnits'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
iarea=2;
LatRange=[0,66 ;66, 200];
fband=[40 80];
for irange=1:2
for iCellType=1
    celltype=CellTypeName{iCellType};
    Ncell=zeros(2);
    figure;
    for isexo=0:1             
        C1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_despike_linear_SU.mat']);  %
        C2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_despike_linear_SU.mat']);
        V1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_despike_linear_SU.mat']);
        V2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_despike_linear_SU.mat']);                
        Vars=fieldnames(C1.C); 
        f=C1.C.f;
        subplot(2,2,isexo*2+1);
        clr2=[1,0,0;0,0,1];
        icnd=0;
        baseline=zeros(1,2);
        mSFC=cell(1,2);
        for ivar=[1 3]
            Var=Vars{ivar};
            coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});
            latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
            param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];
            icnd=icnd+1;
            t=C1.C.t_cue*1000;
            idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
            idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
            Ncell(isexo+1,1)=sum(idx&idx0);
            idx_f=f>=fband(1)&f<fband(2);
            mm=squeeze(mean(coh(idx&idx0,:,idx_f),3));            
            baseline(icnd)=mean(mean(mm(:,t>-200&t<-50),2),1);
            mm=(mm-baseline(icnd))/baseline(icnd)*100;
            mSFC{icnd}=mm;
            mm_f=mean(mm,1);
            ee=std(mm,0,1)/sqrt(size(mm,1));                        
            patchplot(t,mm_f,ee,clr2(icnd,:)); hold on;
            line([-300 500],[0 0],'linestyle','--','linewidth',1.5,'color','k'); hold on;
            title([cuetype{isexo+1} ', cue delay' ]);             
        end
        ss=zeros(1,size(mSFC{1},2));
        pp=zeros(1,size(mSFC{1},2));
        for itime=1:size(mSFC{1},2)
            [pp(itime),ss(itime)]=ranksum(mSFC{1}(:,itime),mSFC{2}(:,itime),'alpha',0.1);            
%             [ss(itime),~]=ttest2(mSFC{1}(:,itime),mSFC{2}(:,itime));            
        end
        plot(t(ss==1),ss(ss==1)*-20,'.k'); hold on;
        xlim([-300 500]);
        ylim([-30 50]);
        
        subplot(2,2,isexo*2+2);
        icnd=0;
        for ivar=[2 4]
            Var=Vars{ivar};
            coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});
            latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
            param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];            
            idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
            idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
            Ncell(isexo+1,2)=sum(idx&idx0);            
            t=C1.C.t_array*1000;
            icnd=icnd+1;            
            mm=squeeze(mean(coh(idx&idx0,:,idx_f),3));                        
            mm=(mm-baseline(icnd))/baseline(icnd)*100;
            mSFC{icnd}=mm;
            mm_f=mean(mm,1);            
            ee=std(mm,0,1)/sqrt(size(mm,1));            
            patchplot(t,mm_f,ee,clr2(icnd,:)); hold on;  
            line([-300 500],[0 0],'linestyle','--','linewidth',1.5); hold on;
            title([cuetype{isexo+1} ', array delay' ]);             
        end                
        ss=zeros(1,size(mSFC{1},2));
        pp=zeros(1,size(mSFC{1},2));
        for itime=1:size(mSFC{1},2)
            [pp(itime),ss(itime)]=ranksum(mSFC{1}(:,itime),mSFC{2}(:,itime),'alpha',0.1);  
%             [ss(itime),~]=ttest2(mSFC{1}(:,itime),mSFC{2}(:,itime));            
        end
        plot(t(ss==1),ss(ss==1)*-20,'.k'); hold on;  
        xlim([-300 500]);
        ylim([-30 50]);
    end
end
end
%% converge the spike_field_coherence between a pair of areas of two monkeys together
close all;
clear all;
AreaName={'LIP','PUL','FEF'};
% folder='Z:\Rujia\Results\';
folder='C:\Data\';
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
irange=1;
LatRange=[0,200; 66, 200];
CellTypeName={'NarrowUnits','BroadUnits'};
iarea=3;
jarea=1;
Nunit=zeros(2);
for icelltype=1:2
    celltype=CellTypeName{icelltype};    
    figure;        
    for isexo=0:1
        C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);  %cueResponsive
        C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);     %AttentionModulated
        V1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);
        V2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);
        Vars=fieldnames(C1.C);
        ivar=1;
        Var=Vars{ivar};
        Time1=C1.C.t_cue;                
        data1=cat(1, C1.C.(Var){iarea,jarea}, C2.C.(Var){iarea,jarea});
        param1 = [V1.Vislatency.(Var){iarea,jarea}(:,2);V2.Vislatency.(Var){iarea,jarea}(:,2)];
        Latency1 = [V1.Vislatency.(Var){iarea,jarea}(:,1);V2.Vislatency.(Var){iarea,jarea}(:,1)];
        
        idx0=param1>0.4&Latency1>=LatRange(irange,1)&Latency1<=LatRange(irange,2);
        idx1=~isnan(squeeze(sum(sum(data1,2),3)))&idx0;
        Nunit(icelltype,1)=sum(idx1);
        aa=transpose(squeeze(mean(data1(idx1,:,:),1)));                
        Time=C1.C.t_cue*1000;
        mbaseline=mean(aa(:, Time1>=-0.25& Time1<=-0.1),2);
        baseline=repmat(mbaseline,1,size(aa,2));
       
        aa=(aa-baseline)./baseline*100;
        mm=conv2(aa,ones(4)/16,'same');        
        subplot(2,2,isexo*2+1);
        idxT=Time>=-300&Time<=500;  
        imagesc(Time(idxT), C1.C.f, mm(:,idxT)); hold on;  % -cc(:,idxT)
        line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
        line([-300 500], [15 15],'color','k','linestyle','--' ); hold on;
        line([-300 500], [30 30],'color','k','linestyle','--' ); hold on;
        caxis([-20 20]);
        ylim([5 100]);
        axis xy;
        title([cuetype{isexo+1} ',cue onset']);
        set(gca,'box','off','linewidth',2.0);  %,'fontsize',10
        colorbar;
        
        subplot(2,2,isexo*2+2);
        ivar=2;
        Var=Vars{ivar};                      
        data1=cat(1, C1.C.(Var){iarea,jarea}, C2.C.(Var){iarea,jarea});
        param1 = [V1.Vislatency.(Var){iarea,jarea}(:,2);V2.Vislatency.(Var){iarea,jarea}(:,2)];
        Latency1 = [V1.Vislatency.(Var){iarea,jarea}(:,1);V2.Vislatency.(Var){iarea,jarea}(:,1)];
        
        idx0=param1>0.4&Latency1>=LatRange(irange,1)&Latency1<=LatRange(irange,2);
        idx1=~isnan(squeeze(sum(sum(data1,2),3)))&idx0;
        Nunit(icelltype,2)=sum(idx1);
        aa=transpose(squeeze(mean(data1(idx1,:,:),1)));
        
        Time=C1.C.t_array*1000;        
        baseline=repmat(mbaseline,1,size(aa,2));                
        aa=(aa-baseline)./baseline*100;
        mm=conv2(aa,ones(4)/16,'same');                
        idxT=Time>=-300&Time<=500;  
        imagesc(Time(idxT), C1.C.f, mm(:,idxT)); hold on;  % -cc(:,idxT)
        line([0 0], [0 100],'color','k','linestyle','--' ); hold on;
        line([-300 500], [15 15],'color','k','linestyle','--' ); hold on;
        line([-300 500], [30 30],'color','k','linestyle','--' ); hold on;
        caxis([-20 20]);
        ylim([5 100]);
        axis xy;
        title([cuetype{isexo+1} ',array onset']);
        set(gca,'box','off','linewidth',2.0,'fontsize',10);
        colorbar;
    end
end

%% Compare the difference of SFC during some frequency bands between early and late cell groups among pulvinar and other regions
close all;
clear all;
% folder='Z:\Rujia\Results\';
folder='C:\Data\';
CellTypeName={'NarrowUnits','BroadUnits'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};
iarea=1;
jarea=3;
LatRange=[0,200 ;66, 200];
fband=[40 80];

for iCellType=1:2
    celltype=CellTypeName{iCellType};
    Ncell=zeros(2);
    for irange=1
        figure;
        for isexo=0:1
            C1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);  %cueResponsive
            C2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_between_areas_AllChannelReference_' celltype '_SU.mat']);     %AttentionModulated
            V1=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);
            V2=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_between_area_' celltype '_SU.mat']);
            Vars=fieldnames(C1.C);
            f=C1.C.f;
            subplot(2,2,isexo*2+1);
            clr2=[1,0,0;0,0,1];
            icnd=0;
            baseline=zeros(1,2);
            mSFC=cell(1,2);
            for ivar=[1 3]
                Var=Vars{ivar};
                coh=cat(1,C1.C.(Var){iarea,jarea},C2.C.(Var){iarea,jarea});
                latency=[squeeze(V1.Vislatency.(Var){iarea,jarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea,jarea}(:,1))];
                param=[squeeze(V1.Vislatency.(Var){iarea,jarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea,jarea}(:,2))];
                icnd=icnd+1;
                t=C1.C.t_cue*1000;
                idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
                idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
                Ncell(irange,isexo+1)=sum(idx&idx0);
                idx_f=f>=fband(1)&f<fband(2);
                mm=squeeze(mean(coh(idx&idx0,:,idx_f),3));
                baseline(icnd)=mean(mean(mm(:,t>-200&t<-50),2),1);
                mm=(mm-baseline(icnd))/baseline(icnd)*100;
                mSFC{icnd}=mm;
                mm_f=mean(mm,1);
                ee=std(mm,0,1)/sqrt(size(mm,1));
                patchplot(t,mm_f,ee,clr2(icnd,:)); hold on;
                line([-300 500],[0 0],'linestyle','--','linewidth',1.5,'color','k'); hold on;
                title([cuetype{isexo+1} ', cue delay' ]);
            end
            ss=zeros(1,size(mSFC{1},2));
            pp=zeros(1,size(mSFC{1},2));
            for itime=1:size(mSFC{1},2)
                [pp(itime),ss(itime)]=ranksum(mSFC{1}(:,itime),mSFC{2}(:,itime),'alpha',0.1);
                %             [ss(itime),~]=ttest2(mSFC{1}(:,itime),mSFC{2}(:,itime));
            end
            plot(t(ss==1),ss(ss==1)*-10,'.k'); hold on;
            xlim([-300 500]);
            ylim([-20 30]);
            
            subplot(2,2,isexo*2+2);
            icnd=0;
            for ivar=[2 4]
                Var=Vars{ivar};
                coh=cat(1,C1.C.(Var){iarea,jarea},C2.C.(Var){iarea,jarea});
                latency=[squeeze(V1.Vislatency.(Var){iarea,jarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea,jarea}(:,1))];
                param=[squeeze(V1.Vislatency.(Var){iarea,jarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea,jarea}(:,2))];
                idx=latency<=LatRange(irange,2)&latency>LatRange(irange,1);
                idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4;
                t=C1.C.t_array*1000;
                icnd=icnd+1;
                mm=squeeze(mean(coh(idx&idx0,:,idx_f),3));
                mm=(mm-baseline(icnd))/baseline(icnd)*100;
                mSFC{icnd}=mm;
                mm_f=mean(mm,1);
                ee=std(mm,0,1)/sqrt(size(mm,1));
                patchplot(t,mm_f,ee,clr2(icnd,:)); hold on;
                line([-300 500],[0 0],'linestyle','--','linewidth',1.5); hold on;
                title([cuetype{isexo+1} ', array delay' ]);
            end
            ss=zeros(1,size(mSFC{1},2));
            pp=zeros(1,size(mSFC{1},2));
            for itime=1:size(mSFC{1},2)
                [pp(itime),ss(itime)]=ranksum(mSFC{1}(:,itime),mSFC{2}(:,itime),'alpha',0.1);
                %             [ss(itime),~]=ttest2(mSFC{1}(:,itime),mSFC{2}(:,itime));
            end
            plot(t(ss==1),ss(ss==1)*-10,'.k'); hold on;
            xlim([-300 500]);
            ylim([-20 30]);
        end
    end
end


%% scatter plot the averaged SFC in different frequency band for individual unit
CellTypeName={'NarrowUnits','BroadUnits'};
iCellType=1;
celltype=CellTypeName{iCellType}; 
window=[0,200];
Gamma_f=[30,40;50,70;80,100];
close all;
band=[30,100;16,25;8,15];
for isexo=0
    figure;
    Ncell=zeros(1,3);
    C1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_SU.mat']);
    C2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_spike-field-coherence_within_area_' celltype '_SU.mat']);
    V1=load([folder 'Vasco_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_SU.mat']);
    V2=load([folder 'Mikey_' cuetype{isexo+1} '_' groupName '_Visual_latency_for_SFC_within_area_' celltype '_SU.mat']);
    
    rr=zeros(3);
    pp=zeros(3);
    for iarea=1:3
        Var='CueOn';  %'ArrayOn';  %
        coh=cat(1,C1.C.(Var){iarea},C2.C.(Var){iarea});        
        latency=[squeeze(V1.Vislatency.(Var){iarea}(:,1)); squeeze(V2.Vislatency.(Var){iarea}(:,1))];
        param=[squeeze(V1.Vislatency.(Var){iarea}(:,2)); squeeze(V2.Vislatency.(Var){iarea}(:,2))];
        Ncell(iarea)=length(latency);
        t_cue=C1.C.t_cue*1000;
        mm_f_band=cell(1,3);
        idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4&latency>20;        
        
        for iband=1:3            
            if iband ==1                
                idx_f=C1.C.f>=Gamma_f(iarea,1) &C1.C.f<=Gamma_f(iarea,2);
            else
                idx_f=C1.C.f>=band(iband,1) &C1.C.f<=band(iband,2);
            end         
            mm_f_band{iband}=squeeze(mean(coh(idx0,:,idx_f),3));                       
            if strcmp(Var,'ArrayOn')
                t=C1.C.t_array*1000;
                coh1=cat(1,C1.C.CueOn{iarea},C2.C.CueOn{iarea});
                latency1=[squeeze(V1.Vislatency.CueOn{iarea}(:,1)); squeeze(V2.Vislatency.CueOn{iarea}(:,1))];
                param1=[squeeze(V1.Vislatency.CueOn{iarea}(:,2)); squeeze(V2.Vislatency.CueOn{iarea}(:,2))];
                
                idx1=~isnan(squeeze(mean(mean(coh1,2),3)))&param1>0.4&latency1>20;                
                nn_f=squeeze(mean(coh1(idx1,:,idx_f),3));
                baseline = repmat(mean(nn_f(:,t_cue>-200&t_cue<-50),2),1,size(mm_f_band{iband},2));                
            else
                t=t_cue;                
                baseline=repmat(mean(mm_f_band{iband}(:,t>-200&t<-50),2),1,size(mm_f_band{iband},2));                
            end
            Lat_valid=latency(idx0);
            mm_f_band{iband}=(mm_f_band{iband}-baseline)./baseline*100;            
            
            subplot(3,3,iarea+(iband-1)*3);
            mc_tf=mean(mm_f_band{iband}(:,t>=window(1)&t<=window(2)),2);
            [rr(iband,iarea),pp(iband,iarea)]=corr(mc_tf,Lat_valid,'type','Pearson');
            if pp(iband,iarea)<=0.05
                clr=[1,0,0];
            else
                clr=[0,0,1];
            end
            num=1;
            t_bin=30:10:190;
            mm_bin=zeros(1,length(t_bin));
            for icenter=30:10:190
                idx=Lat_valid>=icenter-10&Lat_valid<=icenter+10;
                if sum(idx)>0
                    mm_bin(num)=mean(mc_tf(idx));
                end
                num=num+1;
            end
            
%             bar(t_bin,mm_bin); hold on;
%             ylim([0.05 0.2])
            plot(Lat_valid, mc_tf, '.','color',clr); hold on;            
            line([-100, 250],[0 0],'linestyle','--','color','k'); hold on;
            line([100 100],[-100 200],'linestyle','--','color','k'); hold on;
            axis tight;
            xlim([0 220]);
%             ylim([-100, 200]);
            text(150,150,{['r= ' num2str(rr(iband,iarea),2)], ['p= ' num2str(pp(iband,iarea),2)]});
            if iband ==1
                title([AreaName{iarea} ':' num2str(Gamma_f(iarea,1)) '-' num2str(Gamma_f(iarea,2)) 'Hz']);
            else
                title([AreaName{iarea} ':' num2str(band(iband,1)) '-' num2str(band(iband,2)) 'Hz']);
            end
        end        
    end
end

%% plot the averaged SFC  during the delay window

CellTypeName={'NarrowUnits','BroadUnits'};
iCellType=2;
celltype=CellTypeName{iCellType}; 
close all;
window=[0,200;300,500];
lat_band=[20,70; 71, 200];
for isexo=0
    figure;
    clr2=[1 0 0; 0 0 1];
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
        f=C1.C.f;
        idx_t0 = t_cue>-200&t_cue<-50;
        for ilat=1:2
            idx0=~isnan(squeeze(mean(mean(coh,2),3)))&param>0.4&latency> lat_band(ilat,1)&latency<= lat_band(ilat,2);        
            if strcmp(Var,'ArrayOn')
                t=C1.C.t_array*1000;
                coh1=cat(1,C1.C.CueOn{iarea},C2.C.CueOn{iarea});
                latency1=[squeeze(V1.Vislatency.CueOn{iarea}(:,1)); squeeze(V2.Vislatency.CueOn{iarea}(:,1))];
                param1=[squeeze(V1.Vislatency.CueOn{iarea}(:,2)); squeeze(V2.Vislatency.CueOn{iarea}(:,2))];    

                idx1=~isnan(squeeze(mean(mean(coh1,2),3)))&param1>0.4&latency1>lat_band(ilat,1)&latency1<= lat_band(ilat,2);
                baseline =squeeze(mean(mean(coh1(idx1,idx_t0,:),2),1));             
            else
                t=t_cue;
                baseline=squeeze(mean(mean(coh(idx0,idx_t0,:),2),1));            
            end
            for iwin=1:2
                idx_t = t>window(iwin,1)&t<window(iwin,2);
                mm=squeeze(mean(coh(idx0,idx_t,:),2));
                mm=mm-repmat(baseline',size(mm,1),1)./repmat(baseline',size(mm,1),1);
                mm_t=mean(mm,1);
                ee=std(mm,0,1)/sqrt(size(mm,1));

                subplot(2,3,iarea+(iwin-1)*3);
                patchplot(f,mm_t,ee,clr2(ilat,:)); hold on;
            end
        end
    end
end
