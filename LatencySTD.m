function lat = LatencySTD(time, data1, bin, figOn)
% time is in units of ms;
% data is in form of cell/trial * Time
if nargin<2
    error('There are at least two inputs!');
elseif nargin<3
    figOn=0;
end
if length(time)~=size(data1,2)
    error('The length of the time and neural data must be the same!')
end
bin=ceil(bin);
data1=data1*1000;
idx0 = time<= 0&time>=-200;
baseR=mean(data1(:,idx0),2);
baseE = std(baseR);
lat=-1;
for itime = 1:500
    sig=[0 0 0];
    for ii=1:3
        idxT=time>=itime+bin*ii-bin*1.5&time<itime+bin*ii-bin*0.5;
        mresp=mean(data1(:,idxT),2);
        sig(ii) = ttest(mresp,baseR,'alpha',0.05,'tail','right');                   
    end
    if abs(sum(sig))==3
        lat=itime;
        break;
    end
end

if lat<30&&lat>0
    figOn=1;
end

if figOn ==1
    figure;
    for icell = 1:size(data1,1)
        data1(icell,:)=conv(data1(icell,:),ones(1,bin)/bin,'same');
    end
   mm=mean(data1(:,1:bin:end),1);            
   ee=std(data1(:,1:bin:end),0,1)/sqrt(size(data1,1));
   time_downsample=time(1:bin:end);
   idx = time_downsample>=-100&time_downsample<=300;
   TimePlt=time_downsample(idx);
   Vx=[TimePlt fliplr(TimePlt) TimePlt(1)];
   yy1=(mm(idx)-ee(idx));
   yy2=(mm(idx)+ee(idx));
   Vy=[yy1 fliplr(yy2) yy1(1)];
   patch(Vx,Vy,'r','edgecolor','r','facealpha',0.5,'edgealpha',0.5); hold on;
   plot(TimePlt, mm(idx),'color','r','linewidth', 1.5); hold on;   
   line([-100 300],[mean(baseR)+baseE*1.96, mean(baseR)+baseE*1.96],'color','k','linestyle','--'); hold on;
   line([-100 300],[mean(baseR)-baseE*1.96, mean(baseR)-baseE*1.96],'color','k','linestyle','--'); hold on;
   ylim([0 50]);
   xlim([-100 300]);
   set(gca,'linewidth',2,'fontsize',10);   
end
