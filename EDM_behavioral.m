%%%%% load Behavioral data 'LOG_B_VascoXXXX.txt';
date='0619';
filepath=['Y:\Rujia\Logfiles\Training\061919\LOG_B_Vasco' date '.txt'];
fileID=fopen(filepath, 'r');
formatSpec='%f %f %f %f %f %f %f %f %f %f';
startRow=1;
TrlInfo=textscan(fileID, formatSpec,  'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

TrlID=TrlInfo{1};
CueType=TrlInfo{2}; %1 is exo; 2 is endo;
CueValid=TrlInfo{3};  % 1 is valid; 2 is invalid same object; 3 is invalid different object
CueDelay=TrlInfo{4};
isCorrect=TrlInfo{5};
ReactionT=TrlInfo{6};
isTargetOnestExo=TrlInfo{7};  % for exo trials
isTargetOnestEndo=TrlInfo{8};  % for exo trials
EndoCueLoc=TrlInfo{9};  % 0 is for exo trials
TrlDuration=TrlInfo{10};


%% load all infomation 'LOG_VascoXXXX.txt';
date='0918';
filepath=['Y:\Rujia\Logfiles\Training\' date '19\LOG_Vasco' date '.txt'];
fileID=fopen(filepath, 'r');
formatSpec='%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
startRow=1;
TrlInfo=textscan(fileID, formatSpec,  'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

TrlID=TrlInfo{1};
isEndo=TrlInfo{2}; % 1 is exo; 2 is endo
CueValid=TrlInfo{3};  % 1 is valid; 2 is invalid same object; 3 is invalid different object
CueDelay=TrlInfo{6};
TargetDelay=TrlInfo{7};
isCorrect=TrlInfo{8};  %1 is eye fix error; 2 is early release; 3 is correct; 4/5 is late release 5  
ReactionT=TrlInfo{9};

isTargetOnestExo=TrlInfo{10};  % for exo trials, whether target has been onset
isTargetOnestEndo=TrlInfo{11};  % for exo trials

BarOrientation=TrlInfo{12}; % 1 is vertical, 2 is horizontal
isCatchTrl=TrlInfo{13};  % 1 is catch trial; >=2 is normal trial
CuePosInCatchExo=TrlInfo{14};
TargetPosExo=TrlInfo{15};

CuePosInCatchEndo=TrlInfo{17};
TargetPosEndo=TrlInfo{18};

TrlDuration=TrlInfo{19};


%% plot behavioral results: hit rate in different locations
close all;
StartTrl=find(TrlID==1);
TrlNum=diff(StartTrl);
EndTrl=[TrlID(StartTrl(2:end)-1); TrlID(end)];
SessionId=find(EndTrl>100);
TrlOrder=1:numel(TrlID);
HitRate=zeros(numel(SessionId),4);

ValidTrl.TargetDelay=[];
ValidTrl.TargetPosExo=[];
ValidTrl.isCorrect=[];
ValidTrl.CueValid=[];

icue=1;
for isession=1:numel(SessionId)
    figure;
    idxSession=TrlOrder>=StartTrl(SessionId(isession))&TrlOrder<=StartTrl(SessionId(isession))+EndTrl(SessionId(isession))-1;
    for iloc=1:4
        idx=isEndo==1 &isCatchTrl<2 &isTargetOnestExo==1&CuePosInCatchExo==iloc & idxSession'==1;  % for catch trials
        idxHit=isEndo==1 &isCatchTrl<2 &isTargetOnestExo==1 &isCorrect==3&CuePosInCatchExo==iloc & idxSession'==1;    
        
%         idx00=isEndo==1 &isCatchTrl>=2 &isTargetOnestExo==1 ; % for detection trials
%         idx=idx00&CueValid==icue&TargetPosExo==iloc&idxSession'==1;  
%         idxHit=idx00&CueValid==icue&isCorrect==3&TargetPosExo==iloc&idxSession'==1;           
        HitRate(isession,iloc)=sum(idxHit)/sum(idx);
    end
     
    bar(1:4, HitRate(isession,:));
    xlabel('Cue position');
    ylabel('Hit rate');
    set(gca,'box','off');
    
    ValidTrl.TargetDelay=[ValidTrl.TargetDelay; TargetDelay(idxSession'&idx00)];
    ValidTrl.TargetPosExo=[ValidTrl.TargetPosExo; TargetPosExo(idxSession'&idx00)];
    ValidTrl.isCorrect=[ValidTrl.isCorrect; isCorrect(idxSession'&idx00)];
    ValidTrl.CueValid=[ValidTrl.CueValid; CueValid(idxSession'&idx00)];    
end
% save(['Y:\Rujia\Logfiles\Training\ValidTrlInfo\ValidTrl_' date '.mat'] ,'ValidTrl','-v7.3');

%% plot behavioral results: hit rate in valid and invalid cue conditions
close all;
for isession=1:numel(SessionId)
    figure;
    HitRate=zeros(4,3);
    idxSession=TrlOrder>=StartTrl(SessionId(isession))&TrlOrder<=StartTrl(SessionId(isession))+EndTrl(SessionId(isession))-1;
    
    for iloc=1:4
        for icue=1:3
            idx=isEndo==1 &isCatchTrl>=2 &isTargetOnestExo==1&CueValid==icue&idxSession'==1&TargetPosExo==iloc; % &BarOrientation==2;    
            idxHit=isEndo==1 &isCatchTrl>=2 &isTargetOnestExo==1 &isCorrect==3&CueValid==icue&idxSession'==1&TargetPosExo==iloc;  % &BarOrientation==2;     
            HitRate(iloc,icue)=sum(idxHit)/sum(idx);
        end
        
        subplot(2,2,iloc);
        bar(1:3, HitRate(iloc,:));
        xticks(1:3);
        xticklabels({'Cued','SO','DO'});
        ylabel('Hit rate');
        set(gca,'box','off');
    end
    
end

%% get the rhythmicity of behavioral oscillation for individual session
close all;
Tcenter=425:10:1775;

params.tapers=[2,3];
params.pad=1;
params.Fs=50;
params.fpass=[1 30];
params.err=[1 0.05];
params.trialave=1;

isession=1;
idxSession=TrlOrder>=StartTrl(SessionId(isession))&TrlOrder<=StartTrl(SessionId(isession))+EndTrl(SessionId(isession))-1;
HitRate=zeros(4,length(Tcenter));
icue=1;
for iloc=1:4
    num=0;
    for tdelay=Tcenter
        num=num+1;
        idx=isEndo==1&isCatchTrl>=2&isTargetOnestExo==1&TargetDelay>=tdelay-15&TargetDelay<tdelay+15&CueValid==icue&TargetPosExo==iloc&idxSession'==1;        
        idxHit=isEndo==1&isCatchTrl>=2&isTargetOnestExo==1&TargetDelay>=tdelay-15&TargetDelay<tdelay+15&CueValid==icue&isCorrect==3&TargetPosExo==iloc&idxSession'==1;
        if sum(idxHit)>0
            HitRate(iloc, num)=sum(idxHit)/sum(idx);
        end
        
    end
    
    ss=HitRate(iloc,:)-mean(HitRate(iloc,:));
    [S,f,~]=mtspectrumc(ss', params);
    if iloc==2||iloc==3
        ifig=5-iloc;
    else
        ifig=iloc;
    end
    subplot(2,2,ifig);
%     plot(Tcenter, HitRate(iloc,:));
%     xlabel('Cue-target delay (ms)');
%     ylabel('Hit rate');
    plot(f,S, 'k', 'linewidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['Cue Loc = ' num2str(iloc)]);
    
end
%%%% converge trials in all 4 locations
num=0;
icue=1;  % 1 is valid; 2 is invalid same object; 3 is invalid different object
HitRateAll=zeros(1, length(Tcenter));
for tdelay=Tcenter
    num=num+1;
    idx=isEndo==1&isCatchTrl>=2&isTargetOnestExo==1&TargetDelay>=tdelay-15&TargetDelay<tdelay+15&CueValid==icue;   %&idxSession'==1;
    idxHit=isEndo==1&isCatchTrl>=2&isTargetOnestExo==1&TargetDelay>=tdelay-15&TargetDelay<tdelay+15&CueValid==icue&isCorrect==3;  %&idxSession'==1;
    if sum(idxHit)>0
        HitRateAll(num)=sum(idxHit)/sum(idx);
    end
    
end
figure;
subplot(2,1,1);
plot(Tcenter, HitRateAll, 'linewidth', 1.5); hold on;
line([400 2000],[0 0],'linestyle','--'); hold on;
set(gca,'box','off');
xlim([400 2000])
xlabel('Cue-target delay (ms)');
ylabel('Hit rate');

ss=HitRateAll-mean(HitRateAll);
[S,f,~]=mtspectrumc(ss', params);

subplot(2,1,2);
plot(f,S, 'k', 'linewidth', 1.5);
set(gca,'box', 'off');
xlabel('Frequency (Hz)');
ylabel('Power');

%% converge all trials in different sessions
clear all; close all;
folder='Y:\Rujia\Logfiles\Training\ValidTrlInfo\';
filepath=dir([folder 'ValidTrl*.mat']);
load([folder filepath(1).name]) ;
Var=fieldnames(ValidTrl);

for ivar=1:numel(Var)
    AllInfo.(Var{ivar})=[];
    for ifile=5:numel(filepath)
        load([folder filepath(ifile).name]) ;
        AllInfo.(Var{ivar})=[AllInfo.(Var{ivar}); ValidTrl.(Var{ivar})];
    end
end

%% plot the hit rate at different locations and on different cued conditions
close all;
icue=1;
HitLocRatio=zeros(1,4);
for iloc=1:4
    idx=AllInfo.CueValid==icue&AllInfo.TargetPosExo==iloc;
    idxHit=AllInfo.CueValid==icue&AllInfo.isCorrect==3&AllInfo.TargetPosExo==iloc;
    HitLocRatio(iloc)=sum(idxHit)/sum(idx);    
end
figure;
bar(1:4, HitLocRatio);
xlabel('Cued location')
ylabel('Hit rate');


HitRatio=zeros(1,3);
for icue=1:3
    idx=AllInfo.CueValid==icue;
    idxHit=AllInfo.CueValid==icue&AllInfo.isCorrect==3;
    HitRatio(icue)=sum(idxHit)/sum(idx);
end
figure;
plot(1:3, HitRatio,'-*','linewidth',2.0); hold on;
%     xlabel('Target location')
ylabel('Hit rate');
xticks(1:3);
xticklabels({'Cued','SO','DO'});
set(gca,'box','off');
 
%% get the rhythmicity of behavioral oscillation
close all;
Tcenter=425:10:1775;

params.tapers=[2,3];
params.pad=1;
params.Fs=50;
params.fpass=[1 40];
params.err=[1 0.05];
params.trialave=1;

HitRate=zeros(4,length(Tcenter));
icue=1;
figure;
clr4=colormap(hsv(4));
for iloc=1:4
    num=0;
    for tdelay=Tcenter
        num=num+1;
        idx=AllInfo.TargetDelay>=tdelay-25&AllInfo.TargetDelay<tdelay+25&AllInfo.CueValid==icue&AllInfo.TargetPosExo==iloc;        
        idxHit=AllInfo.TargetDelay>=tdelay-25&AllInfo.TargetDelay<tdelay+25&AllInfo.CueValid==icue&AllInfo.isCorrect==3&AllInfo.TargetPosExo==iloc;
        if sum(idxHit)>0
            HitRate(iloc, num)=sum(idxHit)/sum(idx);
        end  
    end
    
    ss=HitRate(iloc,:)-mean(HitRate(iloc,:));
    [S,f,~]=mtspectrumc(ss', params);
    if iloc==2||iloc==3
        ifig=5-iloc;
    else
        ifig=iloc;
    end
    subplot(2,2,ifig);
%     plot(Tcenter,ss  , 'color','b', 'linewidth', 1.5); hold on;  %   HitRate(iloc,:)  clr4(iloc,:)
%     xlabel('Cue-target delay (ms)');
%     ylabel('Hit rate');
%     xlim([400 1800]);
%     line([400 1800], [0 0],'linestyle','--');
    plot(f,S, 'k', 'linewidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    xlim([0 25]);
    title(['Cue Loc = ' num2str(iloc)]);
    set(gca,'box','off');
end

%% 
%%%% converge trials in all 4 locations
close all;
figure;
clr3=eye(3);
for icue=2:3
    HitRateAll=zeros(1, length(Tcenter));
    num=0;
    for tdelay=Tcenter
        num=num+1;
        idx=AllInfo.TargetDelay>=tdelay-25&AllInfo.TargetDelay<tdelay+25&AllInfo.CueValid==icue;  
        idxHit=AllInfo.TargetDelay>=tdelay-25&AllInfo.TargetDelay<tdelay+25&AllInfo.CueValid==icue&AllInfo.isCorrect==3;  
        if sum(idxHit)>0
            HitRateAll(num)=sum(idxHit)/sum(idx);
        end    
    end
%     figure;
%     subplot(2,3,icue);
    ss=HitRateAll-mean(HitRateAll);
    pp(icue-1)=plot(Tcenter, ss, 'linewidth', 1.5,'color', clr3(icue,:)); hold on;    
end
    line([400 1800],[0 0],'linestyle','--');
    legend(pp,{'SO', 'DO'});
    xlim([400 1800]);
    set(gca,'box','off');
    xlabel('Cue-target delay (ms)');
    ylabel('Hit rate');
    
%     subplot(2,3,3+icue);
%     [S,f,~]=mtspectrumc(ss', params);    
%     plot(f,S, 'k', 'linewidth', 1.5);
%     xlim([0 25]);
%     set(gca,'box', 'off');
%     xlabel('Frequency (Hz)');
%     ylabel('Power');    
% end
 
 