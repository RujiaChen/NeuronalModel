AIdata=D.adjDirects;
AIdata=AIdata(~isnan(AIdata));
RefreshT=allEventTimes(EventCode==1)*1000;
%%
BlackScreen=AIdata>2000;
MaxTrl=5000;
OnsetT=zeros(1,numel(Correct.RewardT));
ChangeT=diff(BlackScreen);
CueOn=zeros(1,numel(Correct.RewardT));
CueOnsetT=zeros(1,numel(Correct.RewardT));
refreshIdx=zeros(1,numel(Correct.RewardT));
for itrl=1:numel(Correct.RewardT)
    TimeIdx=floor(Correct.RewardT(itrl));
%     idx=find(BlackScreen(TimeIdx-MaxTrl:TimeIdx)==1,1,'last');    
    idx=find(ChangeT(TimeIdx-MaxTrl:TimeIdx)==-1,1,'last');    
    OnsetT(itrl)=idx+TimeIdx-MaxTrl-1;
    idx00=find(ChangeT(OnsetT(itrl)-1400:OnsetT(itrl)-500)==-1,1);
    if ~isempty(idx00)
        CueOn(itrl)=1;
        CueOnsetT(itrl)=idx00+OnsetT(itrl)-1400-1;
    end   
    refreshIdx(itrl)=find(RefreshT<TimeIdx,1,'last');
end

%%
close all;
diffTime=round(OnsetT-Correct.ArrayT);
figure;
% subplot(2,1,1);
TrlOrder=1:length(diffTime);
plot(TrlOrder,diffTime);  hold on;
% plot(TrlOrder(CueOn==1),diffTime(CueOn==1),'r.'); hold on;

subplot(2,1,2);
% plot(Correct.ArrayT,diffTime,'.');  hold on;
% p=polyfit(Correct.ArrayT,diffTime,1);
% yy=polyval(p,Correct.ArrayT);
% plot(Correct.ArrayT,yy,'r'); hold on;

plot(refreshIdx,diffTime,'.');  hold on;
p=polyfit(refreshIdx(9:end),diffTime(9:end),1);
yy=polyval(p,refreshIdx);
plot(refreshIdx,yy,'r'); hold on;
save('Z:\RujiaChen\Results\p_arrayOn_trigger_correct.mat','p','-v7.3');

%%
TimeCost=floor(refreshIdx*p(1)*60/1000)*1000/60+30;
figure;
plot(refreshIdx,diffTime,'b','linewidth',1.5);  hold on;
plot(refreshIdx,TimeCost,'r','linewidth',1.5);  hold on;
legend('Raw','Estimated');

xlabel('Trial Order');
ylabel('Onset delay to trigger');
set(gca,'box','off');

%%
idx=find(CueOnsetT~=0);
cuety=unique(Correct.Cue(idx));
delayCue=CueOnsetT(idx)-Correct.CueT(idx);

figure;

plot(delayCue(Correct.Cue(idx)==8),'r.'); hold on;
plot(delayCue(Correct.Cue(idx)==32),'b.'); hold on;

legend('Exo','Endo');
xlabel('Trial order');
ylabel('Cue onset delay');
set(gca,'box','off')




