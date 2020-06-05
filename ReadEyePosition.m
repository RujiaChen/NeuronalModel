clear all;
close all;
TimeAlignment={'ToArray', 'ToCue'};   
CueType={ 'Exo', 'Endo'};
monkeyset={'Vasco','Mikey'};  

for imonkey=1:2
    monkey=monkeyset{imonkey};
    if strcmp(monkey, 'Mikey')
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        folder='Y:\anne\monkeys\Mikey\logfiles\';
    elseif strcmp(monkey, 'Vasco')
        dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        %     dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
        folder='Y:\anne\monkeys\presentation files\Vasco\logfiles\';
    end
    for ialign=2 % :2
        TimeZero = TimeAlignment{ialign};
        for iCueCnd=1:2
            CueCondition=CueType{iCueCnd};
            MSinterval=[];
            if strcmp(CueCondition, 'Exo')
                isexo=1;
            else
                isexo=0;
            end
            for idate=1 :numel(dateUsed)
                date=[ '20' dateUsed{idate}(5:6)  dateUsed{idate}(1:4)];
                % filepath='Y:\anne\monkeys\presentation files\Vasco\logfiles\20180510\flanker_eyepos051018145810.txt';
                subfolder=[folder date '\'];
                filename=dir([subfolder 'flanker_eyepos*']);
                bytes=zeros(1, numel(filename));
                for ilog=1:numel(filename)
                    bytes(ilog)=filename(ilog).bytes;
                end
                [~, largestFile]=max(bytes);
                filepath=[subfolder filename(largestFile).name];
                fileID=fopen(filepath, 'r');
                formatSpec='%f %f %f %f %f %f %f %f %f ';
                startRow=1;
                TrlInfo=textscan(fileID, formatSpec,  'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                fclose(fileID);
                
                % % %% plot the eye trace on one particular condition
                load(['Z:\RujiaChen\Results\flanker_TrlParam_' dateUsed{idate} '.mat']);
                load(['Z:\RujiaChen\Results\CorrectTrialInfo_' dateUsed{idate} '_new.mat']);
                load(['Z:\RujiaChen\Results\CorrectTrialParam_' dateUsed{idate} '.mat']);
                TargetPosY=[200,-200];
                cueLoc=[2,5];  %cue-on and off
                saccadeDir=cell(2);
                MS=cell(2);
                MStime=cell(2);
                SacDirection=cell(2);
                CueDelayDuration=cell(2);
                SaccadeTime=cell(2);
                ArrayDelay=cell(2);
                for ipos=1:2
                    %     figure;
                    for icue=1:2
                        idx= find(TrlParam.trial_response==-105&TrlParam.rfyi==TargetPosY(ipos)&TrlParam.cueloc==cueLoc(icue)&TrlParam.isexocue==isexo &TrlParam.is_mouse_trial==0);  %
                        idx2=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
                        if numel(idx)~=numel(idx2)
                            error('The numbers of trials in different parameter files do not match!');
                        end
                        Velocity=cell(1,numel(idx));
                        MoveDirection=cell(1, numel(idx));
                        CueDelayDuration{ipos, icue}=zeros(1, numel(idx));
                        SaccadeTime{ipos,icue}=zeros(1,numel(idx));
                        ArrayDelay{ipos,icue}=zeros(1,numel(idx));
                        SacDirection{ipos, icue}=cell(1,numel(idx));
                        MS{ipos, icue}=cell(1,numel(idx));
                        for itrl=1:numel(idx)
                            %                  figure;
                            iidx=find(TrlInfo{1}==idx(itrl));  %&TrlInfo{5}==3, 1, 'first'
                            RawTime=TrlInfo{6}(iidx);
                            CueT=TrlInfo{6}(find(TrlInfo{1}==idx(itrl)&TrlInfo{5}==8, 1, 'first'));
                            DelayT=TrlInfo{6}(find(TrlInfo{1}==idx(itrl)&TrlInfo{5}==3, 1, 'first'));
                            if (strcmp(monkey, 'Vasco')&& idate<=7) || strcmp(monkey, 'Mikey')  %%%%% this should be modified if 'dateUsed' is changed
                                ArrayT=DelayT+round(Correct.ArrayOnsetT(idx2(itrl)))-Correct.DelayT(idx2(itrl));
                                CueDelayDuration{ipos, icue}(itrl)=Correct.ArrayOnsetT(idx2(itrl))-Correct.CueT(idx2(itrl));
                                ArrayDelay{ipos,icue}(itrl)=TrlParam.holdshape_duration(idx(itrl))-(Correct.ArrayOnsetT(idx2(itrl))-Correct.ArrayT(idx2(itrl)));
                            else
                                ArrayT=DelayT+round(Correct.ArrayT(idx2(itrl)))-Correct.DelayT(idx2(itrl));
                                CueDelayDuration{ipos, icue}(itrl)=Correct.ArrayT(idx2(itrl))-Correct.CueT(idx2(itrl));
                                ArrayDelay{ipos,icue}(itrl)=TrlParam.holdshape_duration(idx(itrl));
                            end
                            
                            t2cue=RawTime-CueT;    % aligned by cue onset
                            t2array=RawTime-ArrayT;   % aligned by array onset
                            
                            
                            Xtrace=TrlInfo{7}(iidx)/17;   %+tt(1) : iidx+tt(end)
                            Xtrace=Xtrace-mean(Xtrace(t2cue>-400&t2cue<0));
                            Xtrace=conv(Xtrace,ones(1,10)/10, 'same');
                            XVel=(Xtrace(4:end-1)+Xtrace(5:end)-Xtrace(2:end-3)-Xtrace(1:end-4))/0.006;
                            %             XVel=(Xtrace(4:end-1)-Xtrace(2:end-3))/0.002;
                            
                            Ytrace=TrlInfo{8}(iidx)/17;   %+tt(1) : iidx+tt(end)
                            Ytrace=Ytrace-mean(Ytrace(t2cue>-400 &t2cue<0));
                            Ytrace=conv(Ytrace,ones(1,10)/10, 'same');
                            YVel=(Ytrace(4:end-1)+Ytrace(5:end)-Ytrace(2:end-3)-Ytrace(1:end-4))/0.006;  %
                            %             YVel=(Ytrace(4:end-1)-Ytrace(2:end-3))/0.002;
                            
                            if strcmp(TimeZero, 'ToCue')
                                TimeWin=[-600 1800];  %relative to cue onset
                                idxT=find(t2cue(3:end-2)>=TimeWin(1)&t2cue(3:end-2)<=TimeWin(2));  %set time window
                                realTime1=reshape(t2cue(idxT+2), 1, []);  %%%%% some eye position missed, try to interp a full trace
                            end
                            
                            if strcmp(TimeZero, 'ToArray')
                                TimeWin=[-900 900]; %[0 1300];  % for Saccade time  %[-900 900]; %for MS detection  %relative to array onset
                                idxT=find(t2array(3:end-2)>=TimeWin(1)&t2array(3:end-2)<=TimeWin(2));  %set time window
                                realTime1=reshape(t2array(idxT+2), 1, []);   %%%% relative to array onset
                            end
                            
                            realVelocity1=reshape(sqrt((XVel(idxT).^2+YVel(idxT).^2)/2), 1, []);
                            idxNonOverlap=diff([realTime1 0])~=0;
                            realTime=realTime1(idxNonOverlap) ; %%%%% remove extra points showing the same time
                            realVelocity=realVelocity1(idxNonOverlap);
                            
                            newTime=TimeWin(1):max(realTime);
                            Velocity{itrl}=interp1(realTime,realVelocity, newTime, 'linear');
                            Velocity{itrl}=conv(Velocity{itrl}, ones(1,10)/10, 'same');  %% smoothing again
                            newXVel=interp1(realTime, XVel(idxT(idxNonOverlap))', newTime, 'linear');
                            newYVel=interp1(realTime, YVel(idxT(idxNonOverlap))', newTime, 'linear');
                            if size(newXVel,1)~=1
                                newXVel=reshape(newXVel, 1, []);
                                newYVel=reshape(newYVel, 1, []);
                            end
                            MoveDirection{itrl}=[newXVel; newYVel];   % saccade direction
                            %                         muV_x=mean(newXVel(newTime<200&~isnan(newXVel)));
                            %                         eeV_x=std(newXVel(newTime<200&~isnan(newXVel)));
                            %                         idx_X=abs(newXVel-muV_x)>3*eeV_x;
                            %
                            %                         muV_y=mean(newYVel(newTime<200&~isnan(newYVel)));
                            %                         eeV_y=std(newYVel(newTime<200&~isnan(newYVel)));
                            %                         idx_Y=abs(newYVel-muV_y)>3*eeV_y;
                            
                            muV=mean(Velocity{itrl}(newTime<200&~isnan(Velocity{itrl})));   %(newTime>-200&newTime<400)
                            eeV=std(Velocity{itrl}(newTime<200&~isnan(Velocity{itrl})));    %(newTime>-200&newTime<400)
                            idxSig=find(abs(Velocity{itrl}-muV)>3*eeV);         % | idx_X |idx_Y
                            
                            idxSd=find(newTime>ArrayDelay{ipos,icue}(itrl)&Velocity{itrl}>max(Velocity{itrl})*3/4,1,'first');
                            SaccadeTime{ipos, icue}(itrl)=newTime(idxSd)-10;
                            
                            SigT=newTime(idxSig);
                            Duration=diff(SigT);
                            DiffPoint=[0 find(Duration>1) length(SigT)];
                            
                            iMS=0;
                            for ipoint=1:numel(DiffPoint)-1
                                if DiffPoint(ipoint+1)-DiffPoint(ipoint) > 20   % saccade duration >20 ms
                                    iMS=iMS+1;
                                    VelPiece=Velocity{ itrl}(idxSig(DiffPoint(ipoint)+1:DiffPoint(ipoint+1)));
                                    [~, maxIdx]=max(VelPiece);
                                    MS{ipos, icue}{itrl}(iMS)=SigT(maxIdx+DiffPoint(ipoint));
                                    SacDirection{ipos, icue}{itrl}=MoveDirection{itrl}(:, ismember(newTime, MS{ipos, icue}{itrl}));
                                end
                            end
                            if iMS>1
                                MSinterval=[MSinterval  diff(MS{ipos, icue}{ itrl})];
                            end
                                                        
                            %%%%%% plot the eye trace or the velocity of eye movement
                            %                     subplot(2,1,1);
                            %                     ValidIdx=abs(Xtrace)<1000&abs(Ytrace)<1000;
                            %                      if strcmp(TimeZero, 'ToCue')
                            %                          tt=t2cue;
                            %                      else
                            %                          tt=t2array;
                            %                      end
                            %                     plot(tt(ValidIdx), Xtrace(ValidIdx));  hold on;
                            %                     if iMS>0
                            %                         plot(MS{ipos, icue}{itrl}, ones(1, numel(MS{ipos, icue}{itrl})),'r.'); hold on;
                            %                     end
                            % %                     plot(tt(3:end-2), XVel);  hold on;
                            % %                     plot(tt(3:end-2), (muV_x+3*eeV_x)*ones(1, length(t2array)-4));  hold on;
                            % %                     plot(tt(3:end-2), (muV_x-3*eeV_x)*ones(1, length(t2array)-4));  hold on;
                            %
                            % %                     muX=mean(Velocity{itrl}(newTime>-200&newTime<400));
                            % %                     eeX=std(Velocity{itrl}(newTime>-200 &newTime<400));
                            % %                     plot(newTime, Velocity{itrl}); hold on;
                            % %                     line(TimeWin, [muX-2*eeX muX-2*eeX], 'linestyle', '--', 'color', 'k'); hold on;
                            % %                     line(TimeWin, [muX+2*eeX muX+2*eeX], 'linestyle', '--', 'color', 'k'); hold on;
                            %                     xlim(TimeWin);
                            % %                     ylim([-3 3 ]);
                            %                     xlabel('Time to array onset (ms)');
                            %                     ylabel('Horizontal position');
                            % %                     ylabel('Horizontal velocity');
                            %                     set(gca,'box', 'off');
                            %
                            %                     subplot(2,1,2);
                            %                     plot(tt(ValidIdx), Ytrace(ValidIdx));  hold on;
                            %                     if iMS>0
                            %                         plot(MS{ipos, icue}{itrl}, ones(1, numel(MS{ipos, icue}{itrl})),'r.'); hold on;
                            %                     end
                            % %                     plot(tt(3:end-2), YVel);  hold on;
                            % %                     plot(tt(3:end-2), (muV_y+3*eeV_y)*ones(1, length(t2array)-4));  hold on;
                            % %                     plot(tt(3:end-2), (muV_y-3*eeV_y)*ones(1, length(t2array)-4));  hold on;
                            %                     xlim(TimeWin);
                            % %                     xlim([-1500 600]);
                            % %                     ylim([-3 3]);
                            %                     xlabel('Time to array onset (ms)');
                            %                     ylabel('Vertical position');
                            % %                     ylabel('Vertical velocity');
                            %                     set(gca,'box', 'off');
                            
                        end
                    end
                end
                %             save(['Z:\RujiaChen\Results\SaccadeTime_'  TimeZero '_' CueCondition '_' dateUsed{idate} '.mat'], 'SaccadeTime', '-v7.3');
                save(['Z:\RujiaChen\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' dateUsed{idate} '.mat'], 'MS', '-v7.3');
                save(['Z:\RujiaChen\Results\SaccadeDirection_'  TimeZero '_' CueCondition '_' dateUsed{idate} '.mat'], 'SacDirection', '-v7.3');
                %             save(['Z:\RujiaChen\Results\CueDelayDuration_' CueCondition '_' dateUsed{idate} '.mat'], 'CueDelayDuration', '-v7.3');
            end
            
            %         %%%% plot the distribution of inter-microsaccade interval
            %         figure;
            %         [hh, cc]=hist(MSinterval,30);
            %         bar(cc, hh/numel(MSinterval));
            %         % xlim([0 1000]);
            %         set(gca, 'box', 'off');
            %         xlabel('Inter-microsaccade interval (ms)');
            %         ylabel('Probability');
            
        end
    end
end


%%
figure;
subplot(1,2,1);
hist(TrlParam.holdshape_duration); hold on;
subplot(1,2,2);
hist(TrlParam.delay_duration);
%% plot the distribution of saccade time (reaction time)
saccadeT=[];
CueCondition='Exo';
TimeZero='ToArray';
monkey='Mikey';
if strcmp(monkey, 'Mikey')
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
elseif strcmp(monkey, 'Vasco')
    dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
end

for idate=1:numel(dateUsed)
    load(['Z:\Rujia\Results\SaccadeTime_'  TimeZero '_' CueCondition '_' dateUsed{idate} '.mat']);
    ss=[SaccadeTime{1,1} SaccadeTime{1,2} SaccadeTime{2,1} SaccadeTime{2, 2}];
    saccadeT=[saccadeT ss];
end
figure;
hist(saccadeT);
xlabel('Saccade time relative to array onset (ms)');
ylabel('Occurences'); 
set(gca, 'box', 'off');
%% compare the delay duration in eye position log and plexon recorded log
% date='112018';
date='042418';  %'050118';   %
load(['Z:\RujiaChen\Results\CueDelayDuration_Exo_' date '.mat']);
load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat']);
load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);

offset=cell(2);
figure;
for ipos=1:2
    for icue=1:2
        subplot(2,2, 2*ipos+icue-2);
        idx=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0;
        delayTime=Correct.ArrayOnsetT(idx)-Correct.CueT(idx);
        offset{ipos, icue}= CueDelayDuration{ipos, icue}-delayTime; 
%         offset{ipos, icue}= CueDelayDuration{ipos, icue}-transpose(CorrectTrialParam.delay_duration(idx)); 
        hist(offset{ipos, icue});
%         xlim([-10 10]);
    end
end

%% scatter plot of MS and the averaged velocity trace
% dateUsed={'042418','051218','050118','050418','050518','042118','051718'};  % for Vasco
dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
% dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
close all;
TimeZero='ToCue';  %'ToArray';  %
figure;
Condition={'Exo', 'Endo'};
for icnd=1
    CueCondition=Condition{icnd};
    MSDirAll=cell(2);
    MSTimeAll=cell(2);
    FirstMS=[];
    interMS=[];
    for idate=1:numel(dateUsed)
        date=dateUsed{idate};
        if exist(['Z:\Rujia\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat'], 'file')
            load(['Z:\Rujia\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat']);
            load(['Z:\Rujia\Results\SaccadeDirection_'  TimeZero '_' CueCondition '_' date '.mat']);
            load(['Z:\Rujia\Results\CueDelayDuration_' CueCondition '_' date '.mat']);
            %         figure;
            for ipos=1:2
                for icue=1:2
%                     subplot(2, 2, 2*ipos+icue-2);
                    for ii=1: numel(MS{ipos, icue})
                        if ~isempty(MS{ipos, icue}{ii})
%                             plot(MS{ipos, icue}{ii}, ii*ones(1, numel(MS{ipos, icue}{ii})), 'k.'); hold on;
%                             set(gca, 'box', 'off') ; hold on;
                            MSTimeAll{ipos, icue}=[MSTimeAll{ipos, icue} MS{ipos, icue}{ii} ];
                            idxMS=find(MS{ipos, icue}{ii}>0&MS{ipos, icue}{ii}<CueDelayDuration{ipos, icue}(ii));   % find MS after cue onset                            
                            if numel(idxMS)>=2  %~isempty(idxMS)   %&&idxMS+1<=numel(MS{ipos, icue}{ii})
                                MSDirAll{ipos, icue}=[MSDirAll{ipos, icue}  SacDirection{ipos, icue}{ii}(:, idxMS(1)) ];
                                FirstMS=[ FirstMS MS{ipos, icue}{ii}(idxMS(1))];
                                interMS=[ interMS MS{ipos, icue}{ii}(idxMS(2))-MS{ipos, icue}{ii}(idxMS(1))];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % %% plot the time course of MS
    % close all;
    % figure;
    % for ipos=1:2
    %     for icue=1:2
    %         subplot(2,2, 2*ipos+icue-2);
    %         hist(MSTimeAll{ipos, icue},30);
    %         set(gca, 'box', 'off');
    %         xlabel('MS time (ms)');
    %         ylabel('Occurence');
    %     end
    % end
    %%%%% converage all data in 4 conditions together
    subplot(1,2,icnd);
    dataAll=[MSTimeAll{1,1} MSTimeAll{1,2} MSTimeAll{2,1} MSTimeAll{2,2}];
    if strcmp(TimeZero, 'ToCue')
        cc=-600:20:1800;
    else
        cc=-900:20:900;
    end
    [hh, cc]=hist(dataAll, cc);
    bar(cc, hh);
    xlim([min(cc)+100 max(cc)-100]);
    ylim([0 50])
    set(gca, 'box', 'off');
    xlabel('MS time (ms)');
    ylabel('Occurence');

end
%% plot the distribution of the Micro saccade directions 
close all;
figure;
xrange=400;
num=0;
for ipos=1:2
    for icue=1:2
        num=num+1;
        subplot(2,2,num);
        idx1=MSDirAll{ipos, icue}(1, :)>=0;
        theta1=mod(atan(MSDirAll{ipos, icue}(2, idx1)./MSDirAll{ipos, icue}(1, idx1)), 2*pi);
        theta2=atan(MSDirAll{ipos, icue}(2, idx1==0)./MSDirAll{ipos, icue}(1, idx1==0))+pi;
        theta=[theta1  theta2];

        cc=0: pi/16 : 2*pi;
        [hh, ~]= hist(theta, cc); 
        polar(cc , hh); hold on;
        polar(pi/4*(-1)^(ipos-1)+(icue-1)*pi, max(hh), 'ro'); hold on;
        
% %%%%% plot the distribution of the saccade velocities
%         xx=MSDirAll{ipos, icue}(1, :);
%         yy=MSDirAll{ipos, icue}(2, :);
%         plot(xx, yy, '.'); hold on;
% %         p=polyfit(xx, yy,1);
% %         x1=linspace(min(xx),max(xx));
% %         y1=polyval(p,x1);
% %         plot(x1, y1, 'r'); hold on;
%         line([-xrange xrange],[0 0],'color','k'); hold on;
%         line([0 0],[-xrange xrange],'color','k'); hold on;
% %         xlim([-xrange xrange]);
% %         ylim([-xrange xrange]);
%         xlabel('Horizontal velocity');
%         ylabel('Vertical velocity');
%         set(gca,'box', 'off');
    end
end

%%%%%% plot the distribution of the first saccade times
figure;
hist(FirstMS, 30);
% xlabel('Time of first MS (ms)');
xlabel('Time of second MS (ms)');
ylabel('Occurence');
xlim([0 1400])
figure;
hist( interMS,30);
xlabel('Inter-MS interval for first two MS');
ylabel('Occurence');


%% time course of the MS direction relative to cue direction
monkeyset={'Vasco','Mikey'};
imonkey=2;
monkey=monkeyset{imonkey};
if strcmp(monkey , 'Vasco')
    dateUsed={ '042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119' };   % for Vasco
elseif strcmp(monkey , 'Mikey')
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
end
close all;
TimeZero='ToArray';  %'ToCue';    %
CueDirection=[pi/4, 5*pi/4; 7*pi/4, 3*pi/4]; 
Condition={'Exo', 'Endo'};
figure;
delayAll=[];
for icnd=1:2
    CueCondition=Condition{icnd};    
    MSAngle=[];
    MSTime=[];
    for idate=1:numel(dateUsed)
        date=dateUsed{idate};
        if exist(['Z:\RujiaChen\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat'], 'file')
            load(['Z:\RujiaChen\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat']);
            load(['Z:\RujiaChen\Results\SaccadeDirection_'  TimeZero '_' CueCondition '_' date '.mat']);
            load(['Z:\RujiaChen\Results\CueDelayDuration_' CueCondition '_' date '.mat']);
            for ipos=1:2
                for icue=1:2
%                     delayAll=[delayAll CueDelayDuration{ipos, icue}];                    
                    for ii=1: numel(MS{ipos, icue})
                        if ~isempty(MS{ipos, icue}{ii})
                            MSTime=[MSTime MS{ipos, icue}{ii}];
                            
                            data=SacDirection{ipos, icue}{ii};
                            theta=zeros(1, size(data,2));
                            idx1=data(1, :)>=0;
                            theta(idx1)=mod(atan(data(2, idx1)./data(1, idx1)), 2*pi);
                            theta(idx1==0)=atan(data(2, idx1==0)./data(1, idx1==0))+pi;
                            Angle=mod(abs(theta-CueDirection(ipos, icue)), 2*pi);                            
                            MSAngle=[MSAngle Angle];
                        end
                    end
                end
            end
        end
    end
    %%%% plot the relative angle to cue direction
    idx00=MSAngle>pi;
    MSAngle(idx00)=2*pi-MSAngle(idx00);
    cc=pi/20:pi/10:pi;
    Nangle=numel(cc);
    Steps=30;
    if strcmp(TimeZero,'ToCue')
        TimeSlot=-600:Steps:1800;  % for cue
    elseif strcmp(TimeZero, 'ToArray')
        TimeSlot=-900:Steps:900;  % for array
    end
    AngleMatrix = zeros(Nangle, numel(TimeSlot));
    IsoDirection=zeros(1, numel(TimeSlot));
    AntiDirection=zeros(1, numel(TimeSlot));
    OrthDirection=zeros(1, numel(TimeSlot));
    Vertical=zeros(1, numel(TimeSlot));
    for itime=1:numel(TimeSlot)
        idx=MSTime>=TimeSlot(itime)&MSTime<TimeSlot(itime)+Steps;
        [hh, cc]= hist(MSAngle(idx), cc);
        AngleMatrix(:, itime)=hh;
        
        IsoDirection(itime)=sum(MSAngle(idx)<pi/4);  %pi/6
        AntiDirection(itime)=sum(MSAngle(idx)>pi*3/4);   %pi*5/6
        OrthDirection(itime)=sum(MSAngle(idx)>=pi/4&MSAngle(idx)<pi*3/4)/2;   %sum(MSAngle(idx)>=pi*3/8&MSAngle(idx)<pi*5/8);
        Vertical(itime)= (sum(MSAngle(idx)>=pi/6&MSAngle(idx)<pi/3)+sum(MSAngle(idx)>=pi*2/3&MSAngle(idx)<pi*5/6))/2;
    end
    
    subplot(1,2,icnd);
    imagesc(TimeSlot, cc*180/pi, AngleMatrix);
    xlabel(['Time to '  lower(TimeZero(3:end)) ' onset (ms)']);
    ylabel('|MS - cue direction|');
    xlim([min(TimeSlot)+100 max(TimeSlot)-200]);
    axis xy;
    caxis([0 15]);
    colorbar;
      
%     smWin=3;
%     plot(TimeSlot, conv(IsoDirection, ones(1,smWin)/smWin, 'same'), 'r', 'linewidth', 1.5); hold on;
%     plot(TimeSlot, conv(AntiDirection, ones(1,smWin)/smWin, 'same'), 'b', 'linewidth', 1.5); hold on;
% %     plot(TimeSlot, conv(OrthDirection, ones(1,smWin)/smWin, 'same'), 'color', [1 1 1]*0.3, 'linewidth', 1); hold on;
%     % plot(TimeSlot, conv(Vertical, ones(1,smWin)/smWin, 'same'),'k', 'linewidth', 1.5); hold on;
%     if strcmp(monkey , 'Vasco')
%         Vx=[1000 1300 1300 1000 1000];
%     elseif strcmp(monkey , 'Mikey')
%         Vx=[700 1300 1300 700 700];
%     end
%     Vy=[0 0 35 35 0];
%     pp=patch(Vx, Vy, [1 1 1]*0.8, 'FaceAlpha', 0.5, 'EdgeColor', 'w'); hold on;    
%     xlim([TimeSlot(1)+100  TimeSlot(end)-100]);
%     ylim([0 15]);
%     xlabel(['Time to '  lower(TimeZero(3:end)) ' onset (ms)']);
%     ylabel('Occurrence of MS');
%     ll=legend('Peri-cue', 'Anti-cue', 'location','best');  %, 'Orth-cue',  'Vertical'  
%     set(ll, 'box', 'off')
%     set(gca,'box', 'off');

end

%% get the neural response during cue delay period on saccade-to and saccade-away conditions
clear all; 
close all;
folder='Z:\RujiaChen\Results\';
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
clr2=[0 0 1; 1 0 0];
CueDirection=[pi/4, 5*pi/4; 7*pi/4, 3*pi/4]; 
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
TimeCue=-600:1200;
TimeArray=-505:1305; 
TimeZero='ToCue';  %'ToArray';  % 
if strcmp(TimeZero, 'ToCue')
    Time=TimeCue;
elseif strcmp(TimeZero, 'ToArray')
    Time=TimeArray;
end

FirstMS=[];
SecondMS=[];
CueCondition='Exo';
Selection={'','_5sigma','_SingleUnit'};
iselect=3;
Group={'VisResp','VisMotor','MotorResp','AllVisual'};  
igroup=4;
groupName=Group{igroup};
CellTypeName={'NarrowUnits','BroadUnits'};
icelltype=2;
cuetype={'endo','exo'};
isexo=1;
for imonkey=1:2
    monkey=monkeyset{imonkey};
    if strcmp(monkey, 'Vasco')
        dateUsed={ '042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119' };   % for Vasco
    elseif strcmp(monkey, 'Mikey')
        dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
    end
    Ncell=[0 0 0];
    for idate= 1:numel(dateUsed)
        date=dateUsed{idate};
        if iselect==1
            load([folder 'Unsorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
        elseif iselect==2
            load([folder 'Unsorted_SpikeTrain_allArea_' date '_5sigma_TriggerCorrected.mat']);
        elseif iselect==3
            load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
        end
        load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
        load([folder 'bVisResp_' date Selection{iselect} '.mat']);
        load([folder 'bMotorResp_' date Selection{iselect} '.mat']);
        load([folder 'bVisMotor_' date Selection{iselect} '.mat']);
        load([folder 'RFOn_' date Selection{iselect} '.mat']);
        load([folder 'SaccadeOn_' date Selection{iselect} '.mat']);
        load([folder 'CueOn_' date Selection{iselect} '.mat']);
        load([folder 'MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat']);
        load([folder 'SaccadeDirection_'  TimeZero '_' CueCondition '_' date '.mat']);
        load([folder 'CueDelayDuration_' CueCondition '_' date '.mat']);
        load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
        
        for iarea=1:3            
            if strcmp(monkey, 'Mikey')
                isession=find(strcmp(RawInfo(1,:), date)==1);
                channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];
                channelID=find(ismember(channel, RawInfo{iarea+1,isession}));
            elseif strcmp(monkey, 'Vasco')
                channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
            end
            
            for ich=channelID   %1:size(SpikeTrain.cueAlign,2)  %
                count=zeros(1,4);
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
                    Ncell(iarea)=Ncell(iarea)+1;
                    if SpikeShape.CellType{iarea}(Ncell(iarea))==icelltype
                        ipos=ceil(CueOn(ich)/2);
                        icue=2-mod(CueOn(ich),2);
                        idxOn=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
                        CueOnResp=squeeze(SpikeTrain.cueAlign(idxOn,ich,:));       %%% trial*time
                        ArrayOn=squeeze(SpikeTrain.arrayAlign(idxOn,ich,:));
                        mArray=conv(mean(ArrayOn,1),ones(1,50)/50,'same');
                        PeakResp=max(mArray(TimeArray>0&TimeArray<200));

                        for itrl=1:numel(idxOn)
                            idx=find(MS{ipos, icue}{itrl}>200 &MS{ipos, icue}{itrl}<min([1000,CueDelayDuration{ipos, icue}(itrl)-100]));  % >200 to remove the initial response
                            if ~isempty(idx)  % && numel(idx)>1   %%%%% distinguish different saccade directions
                                direction=SacDirection{ipos,icue}{itrl}(:, idx);
                                idx1=direction(1,:)>=0;
                                theta=zeros(1, numel(idx));
                                theta(idx1)=mod(atan(direction(2, idx1)./direction(1, idx1)), 2*pi);
                                theta(idx1==0)=mod(atan(direction(2, idx1==0)./direction(1, idx1==0))+pi, 2*pi);
                                Angle=mod(abs(theta-CueDirection(ipos, icue)), 2*pi);
                                Angle(Angle>pi)=2*pi-Angle(Angle>pi);

                                for ii=1:numel(idx)
                                    if Angle(ii)<=pi/4  %%%% MS to cue
                                        count(1)=count(1)+1;
                                        idxT=find(Time-MS{ipos, icue}{itrl}(idx(ii))>=-100&Time-MS{ipos, icue}{itrl}(idx(ii))<=100); % find the time window around MS
                                        CutResp1=(CueOnResp(itrl,idxT)-mean(CueOnResp(itrl,Time>-200&Time<0)));  %/PeakResp;
                                        Response2MS.CueOn_MStoCue{iarea}{Ncell(iarea)}(count(1),:)=CutResp1; % saccade to cued location
                                        %                             FirstMS=[FirstMS; MS{ipos, icue}{itrl}(idx(1))];
                                    elseif Angle(ii)>=pi*3/4
                                        count(2)=count(2)+1;
                                        idxT=find(Time-MS{ipos, icue}{itrl}(idx(ii))>=-100&Time-MS{ipos, icue}{itrl}(idx(ii))<=100);
                                        CutResp2=(CueOnResp(itrl,idxT)-mean(CueOnResp(itrl,Time>-200&Time<0)));  %/PeakResp;
                                        Response2MS.CueOn_MSawayCue{iarea}{Ncell(iarea)}(count(2),:)=CutResp2;  % saccade away cued location
                                        %                             SecondMS=[SecondMS; MS{ipos, icue}{itrl}(idx(2))];
                                    end
                                end
                            end
                        end                  
                    
                        idxOff=find(CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0);
                        CueOffResp=squeeze(SpikeTrain.cueAlign(idxOff,ich,:));       %%% trial*time
                        ArrayOff=squeeze(SpikeTrain.arrayAlign(idxOff,ich,:));
                        mArray=conv(mean(ArrayOff,1),ones(1,50)/50,'same');
                        PeakResp=max(mArray(TimeArray>0&TimeArray<200));

                        for itrl=1:numel(idxOff)
                            idx=find(MS{ipos, 3-icue}{itrl}>200 &MS{ipos, 3-icue}{itrl}<min([1000, CueDelayDuration{ipos, 3-icue}(itrl)-100]));  % >200 to remove the initial response
                            if ~isempty(idx)  % && numel(idx)>1   %%%%% distinguish different saccade directions
                                direction=SacDirection{ipos,3-icue}{itrl}(:, idx);
                                idx1=direction(1,:)>=0;
                                theta=zeros(1, numel(idx));
                                theta(idx1)=mod(atan(direction(2, idx1)./direction(1, idx1)), 2*pi);
                                theta(idx1==0)=mod(atan(direction(2, idx1==0)./direction(1, idx1==0))+pi, 2*pi);
                                Angle=mod(abs(theta-CueDirection(ipos, 3-icue)), 2*pi);
                                Angle(Angle>pi)=2*pi-Angle(Angle>pi);

                                for ii=1:numel(idx)
                                    if Angle(ii)<=pi/4  %%%% MS to cue
                                        count(3)=count(3)+1;
                                        idxT=find(Time-MS{ipos, 3-icue}{itrl}(idx(ii))>=-100&Time-MS{ipos, 3-icue}{itrl}(idx(ii))<=100); % find the time window around MS
                                        CutResp1=(CueOffResp(itrl,idxT)-mean(CueOffResp(itrl,Time>-200&Time<0)));  %/PeakResp;
                                        Response2MS.CueOff_MStoCue{iarea}{Ncell(iarea)}(count(3),:)=CutResp1; % saccade to cued location                                                     
                                    elseif Angle(ii)>=pi*3/4
                                        count(4)=count(4)+1;
                                        idxT=find(Time-MS{ipos, 3-icue}{itrl}(idx(ii))>=-100&Time-MS{ipos, 3-icue}{itrl}(idx(ii))<=100);
                                        CutResp2=(CueOffResp(itrl,idxT)-mean(CueOffResp(itrl,Time>-200&Time<0)));  %/PeakResp;
                                        Response2MS.CueOff_MSawayCue{iarea}{Ncell(iarea)}(count(4),:)=CutResp2;  % saccade away cued location                                                     
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    save([folder monkey '_' cuetype{isexo+1} '_' groupName '_Response2MS_' CellTypeName{icelltype} '_cueResponsive_SU.mat'],'Response2MS','-v7.3');
end

%% plot the neural activities relative to saccade to different directions
close all;
TT=-100:100;
SmoothWin=50;
folder='Z:\RujiaChen\Results\';
monkeyset={'Mikey','Vasco'};
AreaName={'LIP','PUL','FEF'};
groupName='AllVisual';
isexo=1;
icelltype=2;
for imonkey=1:2
    monkey=monkeyset{imonkey};
    load([folder monkey '_' cuetype{isexo+1} '_' groupName '_Response2MS_' CellTypeName{icelltype} '_cueResponsive_SU.mat']);   %  
    
    Vars=fields(Response2MS);
    figure; 
    for iarea=1:3
        count=zeros(1,4);        
        for ivar=1:numel(Vars)
            Ncell=numel(Response2MS.(Vars{ivar}){iarea});
            mCueResp=zeros(Ncell, length(TT));
            for icell=1:Ncell
                if ~isempty(Response2MS.(Vars{ivar}){iarea}{icell})
                    count(ivar)=count(ivar)+1;
                    mCueResp(count(ivar),:)= conv(mean(Response2MS.(Vars{ivar}){iarea}{icell},1),ones(1,SmoothWin)/SmoothWin,'same');
                end                
            end
            mCueResp_population.(Vars{ivar})=mean(mCueResp(1:count(ivar),:),1)*1000;
            ee.(Vars{ivar})=std(mCueResp(1:count(ivar),:),0,1)/sqrt(count(ivar))*1000;
        end
               
        subplot(2,3,iarea);
        patchplot(TT, mCueResp_population.CueOn_MStoCue, ee.CueOn_MStoCue, 'r'); hold on;
        patchplot(TT, mCueResp_population.CueOff_MSawayCue, ee.CueOff_MSawayCue, 'b'); hold on;
        line([0 0],[min(mCueResp_population.CueOff_MSawayCue) max(mCueResp_population.CueOn_MStoCue)+1],'color','k','linestyle','--'); hold on;
        % legend('CueOn','CueAway');
        set(gca,'box','off');
        xlabel('Time to MS onset (ms)');
        ylabel('Firing rate (Hz)');
        title(AreaName{iarea});
        axis tight;
        
        subplot(2,3,iarea+3);
        pp(1)=patchplot(TT, mCueResp_population.CueOn_MSawayCue, ee.CueOn_MSawayCue, 'r'); hold on;
        pp(2)=patchplot(TT, mCueResp_population.CueOff_MStoCue, ee.CueOff_MStoCue, 'b'); hold on;
        line([0 0],[min(mCueResp_population.CueOff_MStoCue) max(mCueResp_population.CueOn_MSawayCue)+1],'color','k','linestyle','--'); hold on;
        ll=legend(pp,'CueOn','CueAway');
        set(ll,'box','off');
        set(gca,'box','off');
        xlabel('Time to MS onset (ms)');
        ylabel('Firing rate (Hz)');
        axis tight;
    end
end
%% plot the neural response (MUA) when saccading to different directions
close all;
figure;
clear fig1; clear fig2;
if strcmp(monkey, 'Vasco')
    clr2=[0 0 1; 1 0 0];
elseif strcmp(monkey, 'Mikey')
    clr2=[ 1 0 0; 0 0 1];
end
TT=-199:199;
% TT=-400:600; 
yrange=[-10 20];
for icue=1:2
    RespTo=[];
    RespAway=[];
    for ipos=1:2
        for idate=1:numel(Resp_MStoCue{ipos, icue})
            if ~isempty(Resp_MStoCue{ipos, icue}{idate})&&size(Resp_MStoCue{ipos, icue}{idate},2)>1
                RespTo=[ RespTo; Resp_MStoCue{ipos, icue}{idate}];
            elseif  ~isempty(Resp_MStoCue{ipos, icue}{idate})&&size(Resp_MStoCue{ipos, icue}{idate},2)==1
                RespTo=[ RespTo; reshape(Resp_MStoCue{ipos, icue}{idate}, 1, []) ];
            end
            
            if ~isempty(Resp_MSawayCue{ipos, icue}{idate})&&size(Resp_MSawayCue{ipos, icue}{idate},2)>1
                RespAway=[ RespAway; Resp_MSawayCue{ipos, icue}{idate}];
            elseif  ~isempty(Resp_MSawayCue{ipos, icue}{idate})&&size(Resp_MSawayCue{ipos, icue}{idate},2)==1
                RespAway=[ RespAway; reshape(Resp_MSawayCue{ipos, icue}{idate}, 1, []) ];
            end
        end
    end
        
        subplot(1,2,1);
%         for icell=1:size(RespTo)
%             RespTo(icell,:) = conv(RespTo(icell,:) , ones(1,20)/20, 'same');
%         end
        idx=find(mean(RespTo,2)~=0&~isnan(RespTo(:,1))&~isinf(RespTo(:,1)));
        mm1=mean(RespTo(idx,:),1);
        ee=std(RespTo(idx,:), [], 1)/sqrt(size(RespTo(idx,:), 1));
        patch([TT fliplr(TT) TT(1)], [mm1-ee fliplr(mm1+ee) mm1(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
        fig1(icue)=plot(TT, mm1, 'color', clr2(icue,:), 'linewidth', 1.5); hold on;    
        
        if strcmp(monkey, 'Vasco')&&icue==2
%             line([0 0],yrange); 
%             ylim(yrange);
            ll=legend(fig1, 'Attend away', 'Attend to');            
            set(ll,'box', 'off', 'location', 'best');
        elseif strcmp(monkey, 'Mikey')&&icue==2
%             line([0 0], yrange);
%             ylim(yrange);
            ll=legend(fig1, 'Attend to', 'Attend away');
            set(ll,'box', 'off', 'location', 'best');
        end        
        xlabel('Time to MS (ms)');
        ylabel('Firing rate (Hz)');   
        set(gca, 'box' ,'off');
        
        subplot(1,2,2);
%         for icell=1:size(RespAway)
%             RespAway(icell,:) = conv(RespAway(icell,:) , ones(1,20)/20, 'same');
%         end
        idx=find(mean(RespAway,2)~=0&~isnan(RespAway(:,1))&~isinf(RespAway(:,1)));
        mm2=mean(RespAway(idx,:),1);
        ee=std(RespAway(idx,:), [], 1)/sqrt(size(RespAway(idx,:), 1));
        patch([TT fliplr(TT) TT(1)], [mm2-ee fliplr(mm2+ee) mm2(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
        fig2(icue)=plot(TT, mm2, 'color',clr2(icue,:), 'linewidth', 1.5,'linestyle', '-'); hold on;  %
        if strcmp(monkey, 'Vasco')&&icue==2
%             line([0 0], yrange);
%             ylim(yrange);
            ll=legend(fig2, 'Attend away', 'Attend to','');
            set(ll,'box', 'off', 'location', 'best');
        elseif strcmp(monkey, 'Mikey')&&icue==2
%             line([0 0], yrange);
%             ylim(yrange);
            ll=legend(fig2, 'Attend to', 'Attend away');
            set(ll,'box', 'off', 'location', 'best');
        end
%         legend([fig1(icue) fig2(icue)],'MS to cue', 'MS away cue');
        xlabel('Time to MS (ms)');
        ylabel('Firing rate (Hz)');
        set(gca, 'box' ,'off');
end
    

%% get the LFP response on saccade-to and saccade-away conditions
clear all;
% close all;
folder='Z:\RujiaChen\Results\';
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');
clr2=[0 0 1; 1 0 0];
CueDirection=[pi/4, 5*pi/4; 7*pi/4, 3*pi/4]; 
monkeyset={'Vasco', 'Mikey'};
areaName={'LIP', 'PUL', 'FEF'};
TimeZero='ToCue';  %'ToArray';  % 
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
imonkey=2;
monkey=monkeyset{imonkey};  
if strcmp(monkey, 'Vasco')
    dateUsed={ '042418','051218','050118','050418','050518','042118','051718', '112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119' };   % for Vasco   
elseif strcmp(monkey, 'Mikey')
    dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
end
TimeCue=-600:1200;
TimeArray=-905:905;
if strcmp(TimeZero, 'ToCue')
    Time=TimeCue;
else
    Time=TimeArray;
end
CueCondition='Exo';
for iarea=1:3
    Area=areaName{iarea};
    Resp_MStoCue=cell(2);
    Resp_MSawayCue=cell(2);
    for idate=1:numel(dateUsed)
        date=dateUsed{idate};
        load([folder 'LFP_allArea_' date '_TriggerCorrected.mat']);
        load([folder 'bArrayResp_' date '.mat']);
        load([folder 'bVisMotor_' date '.mat']);
        load([folder 'CueOn_' date]);
        load(['Z:\RujiaChen\Results\CorrectTrialParam_' date '.mat']);
        if strcmp(monkey, 'Mikey')
            isession=find(strcmp(RawInfo(1,:), date)==1);
            channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];
            channelID=find(ismember(channel, RawInfo{iarea+1,isession}));
        else
             channelID=(iarea-1)*32+1:iarea*32;  % for Vasco
        end
        CueResp=cell(2);
        ArrayResp=cell(2);
        for ipos=1:2
            num=0;
            for ich=channelID
                if  bVisMotor(ich, ipos)>0 % bArrayResp(ich, ipos)>0   %CueOn(ich,ipos)>0    
                    num=num+1;  
                    for icue=1:2
                        idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==1&CorrectTrialParam.is_mouse_trial==0;
                        CueResp{ipos, icue}(num,:,:)=squeeze(LFPinfo.cueAlign(idx1,ich,:));       %%% channel*trial*time
                        ArrayResp{ipos, icue}(num,:,:)=squeeze(LFPinfo.ArrayAlign(idx1,ich,:));
                    end
                end
            end
        end
        
        load(['Z:\RujiaChen\Results\MicroSaccade_time_'  TimeZero '_' CueCondition '_' date '.mat']);
        load(['Z:\RujiaChen\Results\SaccadeDirection_'  TimeZero '_' CueCondition '_' date '.mat']);
        load(['Z:\RujiaChen\Results\CueDelayDuration_' CueCondition '_' date '.mat']);
        
        for ipos=1:2
            for icue=1:2
                sResp1=[];
                sResp2=[];
                count=[0 0];
                if ~isempty(CueResp{ipos, icue})
                    for itrl=1:numel(MS{ipos, icue})
                        idx=find(MS{ipos, icue}{itrl}>200 &MS{ipos, icue}{itrl}<min([700 CueDelayDuration{ipos, icue}(itrl)]));  % >200 to remove the initial response
                        Resp0=squeeze(CueResp{ipos, icue}(:, itrl, :));
                        Resp1=squeeze(ArrayResp{ipos, icue}(:, itrl, :));
                        if ~isempty(idx)  % && numel(idx)>1   %%%%% to be continue, distinguish different saccade directions
                            direction=SacDirection{ipos,icue}{itrl}(:, idx);
                            idx1=direction(1, :)>=0;
                            theta=zeros(1, numel(idx));
                            theta(idx1)=mod(atan(direction(2, idx1)./direction(1, idx1)), 2*pi);
                            theta(idx1==0)=mod(atan(direction(2, idx1==0)./direction(1, idx1==0))+pi, 2*pi);
                            Angle=mod(abs(theta-CueDirection(ipos, icue)), 2*pi);
                            Angle(Angle>pi)=2*pi-Angle(Angle>pi);
                            
                            for ii=1:numel(idx)
                                if Angle(ii)<pi/4
                                    count(1)=count(1)+1;
                                    idxT=find(Time-MS{ipos, icue}{itrl}(idx(ii))>-200&Time-MS{ipos, icue}{itrl}(idx(ii))<200);
%                                     idxT=find(Time>=-400&Time<=600);
                                    if size(Resp0,2)==1
                                        CutResp=reshape(Resp0(idxT)-mean(Resp0(TimeCue>-200&TimeCue<0)), 1, []);
                                    else
                                        CutResp=Resp0(:, idxT)-repmat(mean(Resp0(:, TimeCue>-200&TimeCue<0),2), 1, length(idxT));
                                    end
                                    sResp1(count(1),:,:)=CutResp;
                                end
                                
                                if Angle(ii)>pi*3/4
                                    count(2)=count(2)+1;
                                    idxT=find(Time-MS{ipos, icue}{itrl}(idx(ii))>-200&Time-MS{ipos, icue}{itrl}(idx(ii))<200);
%                                     idxT=find(Time>=-400&Time<=600);
                                    if size(Resp0,2)==1
                                        CutResp2=reshape(Resp0(idxT)-mean(Resp0(TimeCue>-200&TimeCue<0)), 1, []);
                                    else
                                        CutResp2=Resp0(:, idxT)-repmat(mean(Resp0(:, TimeCue>-200&TimeCue<0),2), 1, length(idxT));
                                    end
                                    sResp2(count(2),:,:)=CutResp2;
                                end
                            end
                        end
                    end
                end
                Resp_MStoCue{ipos, icue}{idate}=squeeze(mean(sResp1,1));
                Resp_MSawayCue{ipos, icue}{idate}=squeeze(mean(sResp2,1));
            end
        end
    end
    
    save(['Z:\RujiaChen\Results\LFP_MStoCue_CueAligned_' monkey '_' Area '.mat'], 'Resp_MStoCue', '-v7.3');
    save(['Z:\RujiaChen\Results\LFP_MSawayCue_CueAligned_' monkey '_' Area '.mat'], 'Resp_MSawayCue', '-v7.3');
    
end

%% plot the neural response (LFP) when saccading to different directions
figure;
clear fig1; clear fig2;
if strcmp(monkey, 'Vasco')
    clr2=[0 0 1; 1 0 0];
elseif strcmp(monkey, 'Mikey')
    clr2=[ 1 0 0; 0 0 1];
end
TT=-199:199;
% TT=-400:600; 
NN=zeros(1, 2);
MM=zeros(1, 2);
yrange=[-2 1];
for icue=1:2
    RespTo=[];
    RespAway=[];
    for ipos=1:2
        for idate=1:numel(Resp_MStoCue{ipos, icue})
            if ~isempty(Resp_MStoCue{ipos, icue}{idate})&&size(Resp_MStoCue{ipos, icue}{idate},2)>1
                RespTo=[ RespTo; Resp_MStoCue{ipos, icue}{idate}];
            elseif  ~isempty(Resp_MStoCue{ipos, icue}{idate})&&size(Resp_MStoCue{ipos, icue}{idate},2)==1
                RespTo=[ RespTo; reshape(Resp_MStoCue{ipos, icue}{idate}, 1, []) ];
            end
            
            if ~isempty(Resp_MSawayCue{ipos, icue}{idate})&&size(Resp_MSawayCue{ipos, icue}{idate},2)>1
                RespAway=[ RespAway; Resp_MSawayCue{ipos, icue}{idate}];
            elseif  ~isempty(Resp_MSawayCue{ipos, icue}{idate})&&size(Resp_MSawayCue{ipos, icue}{idate},2)==1
                RespAway=[ RespAway; reshape(Resp_MSawayCue{ipos, icue}{idate}, 1, []) ];
            end
        end
    end
    NN(icue)=size(RespTo,1);
    MM(icue)=size(RespAway,1);
    
%     %%%%% plot the mean and ste of the LFP signals
%     subplot(1,2,1);
%     mm1=mean(RespTo,1);
%     ee=std(RespTo, [], 1)/sqrt(size(RespTo, 1));
%     patch([TT fliplr(TT) TT(1)], [mm1-ee fliplr(mm1+ee) mm1(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     fig1(icue)=plot(TT, mm1, 'color', clr2(icue,:), 'linewidth', 1.5); hold on;    
%     if strcmp(monkey, 'Vasco')
%         ll=legend(fig1, 'Attend away', 'Attend to');
%         line([0 0], yrange);
%         ylim(yrange);
%     else
%         ll=legend(fig1, 'Attend to', 'Attend away');
%         line([0 0], yrange);
%         ylim(yrange);
%     end
%     set(ll,'box', 'off', 'location', 'best');
%     xlabel('Time to MS (ms)');
%     ylabel('Firing rate (Hz)');
%     set(gca, 'box' ,'off');
%     
%     subplot(1,2,2);
%     mm2=mean(RespAway,1);
%     ee=std(RespAway, [], 1)/sqrt(size(RespAway, 1));
%     patch([TT fliplr(TT) TT(1)], [mm2-ee fliplr(mm2+ee) mm2(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     fig2(icue)=plot(TT, mm2, 'color', clr2(icue,:), 'linewidth', 1.5,'linestyle', '-'); hold on;
%     if strcmp(monkey, 'Vasco')
%         ll=legend(fig2, 'Attend away', 'Attend to');
%         line([0 0], yrange);
%         ylim(yrange);
%     else
%         ll=legend(fig2, 'Attend to', 'Attend away');
%         line([0 0], yrange);
%         ylim(yrange);
%     end
%     set(ll,'box', 'off', 'location', 'best');
%     xlabel('Time to MS (ms)');
%     ylabel('Firing rate (Hz)');
%     set(gca, 'box' ,'off');
    
end
    
%% check the numbers of different events
Variables={'nTrl', 'cueType', 'cueX', 'cueY', 'EventCode', 'time', 'posX', 'posY','response'};
load(['Z:\RujiaChen\Results\flanker_TrlParam_' dateUsed{1} '.mat']);
for ii=1:numel(TrlInfo)
    AllInfo{ii}=unique(TrlInfo{ii});
end
nEvt=zeros(1, numel(AllInfo{5}));
for ievent=1:numel(AllInfo{5})
    nEvt(ievent)=numel(unique(TrlInfo{1}(TrlInfo{5}==AllInfo{5}(ievent))));    
end
TrlID=unique(TrlInfo{1}(TrlInfo{5}==AllInfo{5}(1)));
EventTime=cell(1, numel(AllInfo{5}));
SaccadeTrl=zeros(1, nEvt(1));
for ievent=1:numel(AllInfo{5})
    numTrl=0;
    for itrl=1:numel(TrlID)    
        idx=find(TrlInfo{1}==itrl&TrlInfo{5}==AllInfo{5}(ievent),1, 'first');
        if ~isempty(idx)
            numTrl=numTrl+1;
            EventTime{ievent}(numTrl)=TrlInfo{6}(idx);
            if ievent==1
                SaccadeTrl(numTrl)=itrl;
            end
        end
    end
end

%% check the consistency of time recorded in eye position log file and presentation files
load(['Z:\RujiaChen\Results\TrialTiming_' dateUsed{1} '.mat']);
load(['Z:\RujiaChen\Results\CorrectTrialInfo_' dateUsed{1} '_new.mat']);
num=0;
delayWin1=[]; delayWin2=[];
delayWin3=[]; 
for itrl=1:numel(TrlID)
    idx1=find(TrlInfo{1}==itrl&TrlInfo{5}==3,1, 'first');  %event=3, delay;
    idx2=find(TrlInfo{1}==itrl&TrlInfo{5}==4,1, 'first');  %event=4, array;
    idx3=find(TrlInfo{1}==itrl&TrlInfo{5}==8,1, 'first');  %event=8, cue onset;
    
    if ~isempty(idx1)&&~isempty(idx2)
        num=num+1;
        delayWin1(num)=TrlInfo{6}(idx2)-TrlInfo{6}(idx1); % time in eye positions log file
        cueOnset(num)=TrlInfo{6}(idx1)-TrlInfo{6}(idx3);
        delayWin2(num)=TrlParam.delay_duration(itrl);  % parameter set for trials                
%         delayWin3(num)=double(dataRaw{3}(dataRaw{1}==itrl&dataRaw{2}==4)-dataRaw{3}(dataRaw{1}==itrl&dataRaw{2}==3)); % time recorded in presentation        
    end
end
delayWin4=Correct.ArrayT-Correct.DelayT;
figure;
% TT=delayWin2-delayWin1;
% plot(TT);
% hist(delayWin1)
% hist(delayWin2);
plot(delayWin1, delayWin2, '.'); hold on;  % delayWin2 is as same as delayWin3


%% read the flanker_fix.log file
logfolder='Y:\anne\monkeys\presentation files\Vasco\logfiles\20180424\';
filename=[logfolder 'flanker_fix_1.log'];
fileID=fopen(filename, 'r');
formatSpec='%d %s %d %s %d %s %s %s %s %s %s %s %s ';
for ii=1: 5
    tline=fgets(fileID);
    if ii==4
        headInfo=textscan(tline,'%s','delimiter','\t');
    end
end

TrlOnlineInfo=textscan(fileID,formatSpec, 'EndOfLine', '\n');  % 'delimiter',{' ', '\t'},
fclose(fileID);

TrlID=unique(TrlOnlineInfo{3});
CueT=zeros(1, numel(TrlID));
CueDelayT=zeros(1, numel(TrlID));
ArrayDelayT=zeros(1, numel(TrlID));
num=0;
for ii=1:numel(TrlID)        
    if ~isempty(find(strcmp(TrlOnlineInfo{4},'array')&TrlOnlineInfo{3}==ii, 1))  %holdshapedim
        num=num+1;
        CueT(num)=TrlOnlineInfo{5}((strcmp(TrlOnlineInfo{4},'cue') | strcmp(TrlOnlineInfo{4},'endocue'))&TrlOnlineInfo{3}==ii);      % or endocue  
        DelayOnset=TrlOnlineInfo{5}(strcmp(TrlOnlineInfo{4},'delay')&TrlOnlineInfo{3}==ii);    
        ArrayT=TrlOnlineInfo{5}(strcmp(TrlOnlineInfo{4},'array')&TrlOnlineInfo{3}==ii);
        CueDelayT(num)=(ArrayT- DelayOnset)/10;
%         DimT=TrlOnlineInfo{5}(strcmp(TrlOnlineInfo{4},'holdshapedim')&TrlOnlineInfo{3}==ii);
%         ArrayDelayT(num)=(DimT-ArrayT)/10;
    end
end
CueDelayT(CueDelayT==0)=[];
figure;
hist(CueDelayT);



%%

TimeNew=-300:500;
figure;

for ipos=1:2
    subplot(1,2,ipos);
    for icue =1:2
        mm=mean(Resp_MS{ipos, icue},1);
        mm=conv(mm,ones(1,30)/30,'same'); hold on;
         plot(TimeNew, mm, 'color', clr2(icue,:)); hold on;  % for MS aligned response
         
    end
    xlim([-200 400]);
    xlabel('Time to MS (ms)');
     ylabel('Firing rate (Hz)');
     
end




%%
% dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};  % for Vasco
dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};
for idate=1:numel(dateUsed)
    date=dateUsed{idate};
    load(['Z:\RujiaChen\Results\flanker_TrlParam_' date '.mat']);
    NN(idate)=sum(TrlParam.trial_response==-105);
    
    load(['Z:\RujiaChen\Results\CorrectTrialInfo_' date '_new.mat']); 
    MM(idate)=numel(Correct.Cue);    
end


%%
clear all;
close all;
filepath='Y:\anne\monkeys\presentation files\Vasco\logfiles\20180510\flanker_fix_1.log';
fileID=fopen(filepath, 'r');
formatSpec='%d %s %d %s %f %f %f %f %f %f %s %s %f';
num=0;
trlCount=0;
while num<7000
     tline=fgetl(fileID);
     num=num+1;
     if ~isempty(strfind(tline, 'correct'))
         dataRaw=textscan(tline,formatSpec);
         trlCount=trlCount+1;
         CorrectTime(trlCount)=dataRaw{5};
         
         tline=fgetl(fileID);
         num=num+1;
         tline=fgetl(fileID);
         num=num+1;
         if ~isempty(strfind(tline, 'juice'))
             dataRaw=textscan(tline,formatSpec);
             JuiceTime(trlCount)=dataRaw{5};
         end
     end
%     disp(tline); 
end
