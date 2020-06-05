clear all; close all;
% folder='Z:\Rujia\Results\';
folder='C:\Data\';
monkeyset={'Mikey','Vasco'};
cuetype={'endo','exo'};
Group={'VisResp','VisMotor','MotorResp','AllVisual'};
igroup=4;
groupName=Group{igroup};

params.tapers=[2,3];
params.pad=1;
params.Fs=1000;
params.fpass=[1 100];
params.err=[1 0.05];
params.trialave=1;

for imonkey=1%:2
    monkey=monkeyset{imonkey};
    for isexo=0:1
        load([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_cueResponsive_SU.mat']);
        Vars=fields(SpikeLFPInfo);
        for ivar=2:5
            var=Vars{ivar};
            for idate=1:size(SpikeLFPInfo.CueOnSpike,1)
                for iarea=1:3
                    ivar,idate,iarea
                    for icell=1:numel(SpikeLFPInfo.(var){idate,iarea})
                        num=0;
                        data1=transpose(SpikeLFPInfo.(Vars{ivar+4}){idate,iarea}{icell});
                        data1=rmlinesc(data1,params,[], 0, 60);
                        data1=data1-repmat(mean(data1,1),size(data1,1),1);
                        SpikeLFPInfo.(Vars{ivar+4}){idate,iarea}{icell}=transpose(data1);
                        for itrl=1:size(SpikeLFPInfo.(var){idate,iarea}{icell},1)
                            idx=find(SpikeLFPInfo.(var){idate,iarea}{icell}(itrl,101:end-100)>0);
                            if ~isempty(idx)
%                                 num=num+1;
%                                 subplot(6,6,num);
%                                 plot(data1(:,itrl),'b'); hold on;
                                for ii=idx
%                                     plot(ii+100,2,'b+'); hold on;
                                    lfp_raw=data1(ii:ii+200,itrl);
                                    xx=[1:98 104:200]; %+/-2,5-ms
                                    xx1=1:201;
                                    yy1=interp1(xx,lfp_raw(xx),xx1,'linear');
                                    SpikeLFPInfo.(Vars{ivar+4}){idate,iarea}{icell}(itrl,ii+97:ii+103)=yy1(98:104);
                                end
%                                 plot(SpikeLFPInfo.(Vars{ivar+4}){idate,iarea}{icell}(itrl,:),'r'); hold on;
%                                 xlim([600 800]);
                            end
                        end
                    end
                end
            end
        end
        save([folder monkey '_' cuetype{isexo+1} '_' groupName '_SpikeLFPInfo_AllChannelReference_DeSpike_linear_cueResponsive_SU.mat'],'SpikeLFPInfo','-v7.3');
    end
end
%%
figure;
subplot(1,3,1);
mm=mean(LFPs,1);
ee=std(LFPs,0,1)/sqrt(size(LFPs,1));
Time=-100:100;
% plot(Time,LFPs); hold on;
plot(Time,mm,'linewidth',2.0); hold on;
xx=[1:90 110:200];
yy=LFPs(:,xx);
plot(xx-101,mean(yy,1),'*','linewidth',1.5); hold on;
ylim([-0.5, 0]);
xlim([-100 100]);

xx1=1:201;
yy1=LFPs;
for ii=1:size(LFPs,1)
    yy1(ii,:)=interp1(xx,yy(ii,:),xx1,'spline');    
end
subplot(1,3,2);
mm=mean(yy1);
ee=std(yy1,0,1)/sqrt(size(yy1,1));
plot(xx-101,mean(yy,1),'*','linewidth',1.5); hold on;
plot(xx1-101,mm,'linewidth',2.0); hold on;
ylim([-0.5, 0]);
xlim([-100 100]);
title('Interpolation: SPLINE');

yy2=LFPs;
for ii=1:size(LFPs,1)
    yy2(ii,:)=interp1(xx,yy(ii,:),xx1,'linear');    
end
subplot(1,3,3);
mm=mean(yy2);
ee=std(yy2,0,1)/sqrt(size(yy2,1));
plot(xx-101,mean(yy,1),'*','linewidth',1.5); hold on;
plot(xx1-101,mm,'linewidth',2.0); hold on;
ylim([-0.5, 0]);
xlim([-100 100]);
title('Interpolation: LINEAR');






