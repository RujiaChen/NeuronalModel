%% read plexon raw data
clear D;
folder = 'Y:\Rujia\Vasco_data\';
date = '022820';
% folder= 'Y:\anne\monkeys\plexondata\';
% date = '012119';

ifile = 1;
filepath=dir([folder date '\' date '*_merged_spl_noWB_*4sigma.pl2']);

%%
fileName=[folder date '\' filepath(ifile).name];
% fileName=[folder date '\bar_3.pl2'];
muaDataDirRoot=['Z:\RujiaChen\Vasco\' date '\'];   %042118
sessionName='merged';
areaName='SC';   %'PUL';  %
isLoadSpikes=1;
isLoadMua=0;
isLoadLfp=0;
isLoadSpkc=0;
isLoadDirect=0;
spikeChannelPrefix = 'SPK_SPKC';
% spikeChannelPrefix = 'SPK';

channel = 1:32;
spikeChannelsToLoad=channel;
muaChannelsToLoad=channel;
lfpChannelsToLoad = channel;
spkcChannelsToLoad = channel;
directChannelsToLoad = channel;
D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);

%% identify the bar session ID and trim the recording to that session
All_files = dir([folder date '\*.pl2']);
All_files_cell = transpose(struct2cell(All_files));
idx = cellfun(@(x) ~contains(x, 'merge'), All_files_cell(:,1));
Files_no_merge = All_files_cell(idx,:);
Sorted_files = sortrows(Files_no_merge, 3);
Nfile = numel(D.blockStartTimes);
if ifile == 1
    Files_merged = Sorted_files(1:Nfile,:);
elseif ifile == 2
    Files_merged = Sorted_files(end-Nfile+1: end,:);
else
    print('Please check the merged file ID');
end
Bar_ID =  cellfun(@(x) contains(x, 'bar'), Files_merged(:,1));

RFtask = find(Bar_ID == 1);
isession = 2;
if numel(RFtask) == 1
    isession =1;
end
D_session = trimSpikeTimesAndEvents(D, RFtask(isession));
Plexon_session_time =  Files_merged{RFtask(isession), 3};

%%%%%% read the events times in plexon
allEventTimes = floor(cell2mat(D_session.events([1:8])')*10000);
allEventTimes = unique(allEventTimes);
BinaryCode=cell(size(allEventTimes));
EventCode=zeros(size(allEventTimes));
for i = 1:numel(allEventTimes)
    for j = 1:8
        if any(abs(D_session.events{j}*10000 - allEventTimes(i)) < 1)
            BinaryCode{i}=[BinaryCode{i};j];     
        end 
    end
    idigit=1;
    port=0;
    while idigit<=numel(BinaryCode{i})
        if BinaryCode{i}(idigit)<9
            EventCode(i)=2^(BinaryCode{i}(idigit)-1)+ EventCode(i);
        else
            EventCode(i)=2^(BinaryCode{i}(idigit)-1-8)+ EventCode(i);
            port=1;
        end
        idigit=idigit+1;
    end
end

%% check the trial numbers in plexon and log files , get the stimulus onset time
clear FlashOnT;
date_log = [date(5:6) date(1:4)];
log_files = dir([folder date '\' date_log '_test_bars*.txt']);
log_files_cell = transpose(struct2cell(log_files));
idx = cellfun(@(x) ~contains(x, 'EYE'), log_files_cell(:,1));
Param_log = log_files_cell(idx,:);

log_time_offset = cellfun(@(x) sum(abs(datenum(x) - datenum(Plexon_session_time))), Param_log(:,3));
[~, log_file_ID] = min(log_time_offset);

logPresentation = [folder date '\' Param_log{log_file_ID,1}];
fileID=fopen(logPresentation,'r');
tline=fgets(fileID);
headInfo=textscan(tline,'%s','delimiter','\t');
headInfo_new = headInfo{1}(~cellfun(@isempty, headInfo{1}));

formatSpec='%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
startRow=1;

% dataRaw=textscan(fileID,formatSpec);
TrlInfo=textscan(fileID, formatSpec,  'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
idx = ~isnan(TrlInfo{2});
TrlInfo_new= cell(1, numel(TrlInfo));
for icell = 1: numel(TrlInfo)
    TrlInfo_new{icell} = TrlInfo{icell}(idx);
end

for ii = 1:10 
    Param.(headInfo_new{ii})=TrlInfo_new{ii};
end
for jj = 11:numel(headInfo_new)
    Param.(headInfo_new{jj})=TrlInfo_new{jj+1};
end


Trial_plexon = sum(EventCode==1); % event=1 is trial onset
Success_trials = sum(Param.success_trial==1);
assert(Trial_plexon == max(Param.trial_attempt));
assert(sum(EventCode==64) == Success_trials);
    
Valid_trial = unique(Param.trial_attempt);
FlashOnT=[];
TrialOnsetTime=allEventTimes(EventCode==1);
for itrl = 1:numel(Valid_trial)
    if Valid_trial(itrl) < numel(TrialOnsetTime)
        onset_idx = find((EventCode==2) & allEventTimes > TrialOnsetTime(Valid_trial(itrl)) & allEventTimes < TrialOnsetTime(Valid_trial(itrl)+1));
    else
        onset_idx = find((EventCode==2) & allEventTimes > TrialOnsetTime(Valid_trial(itrl)));
    end
    n_repet = sum(Param.trial_attempt == Valid_trial(itrl));
    if ~isempty(onset_idx) && n_repet ~=0
        FlashOnT= [FlashOnT; allEventTimes(onset_idx(1:n_repet))/10];
    end
end

%%%%
clear SpikeTrain;
for itrl=1:numel(FlashOnT)
    Duration=[FlashOnT(itrl)-100 FlashOnT(itrl)+300];  %% the param is recorded by trials
    for ich=1:numel(D.allSpikeStructs)
        idx00=find(D.allSpikeStructs{ich}.ts*1000>=Duration(1)& D.allSpikeStructs{ich}.ts*1000<=Duration(2));
        if ~isempty(idx00)
            spkTime = D.allSpikeStructs{ich}.ts(idx00)*1000-FlashOnT(itrl);
            BinT=-100:10:300; %
            SpikeTrain(itrl,ich,:)=hist(spkTime,BinT);
%             Resp(itrl,ich)=numel(D.allSpikeStructs{ich}.ts(idx00)*1000);
        end
    end 
end

%%
close all;
stim_x= Param.stim_x + Param.x_center;
stim_y= Param.stim_y + Param.y_center;

xx = unique(stim_x);
yy = unique(stim_y);
mResp_x = zeros(length(xx), size(SpikeTrain,2), size(SpikeTrain,3));
Smooth_kernal = ones(1)/1;
for ii = 1: numel(xx)
    id_trl=find(stim_x == xx(ii) & stim_y == Param.y_center);
    mm = squeeze(mean(SpikeTrain(id_trl,:,:),1));
    mResp_x(ii,:,:) = conv2(mm, Smooth_kernal, 'same');
end
mResp_y = zeros(length(yy), size(SpikeTrain,2), size(SpikeTrain,3));
for ii = 1: numel(yy)
    id_trl=find(stim_y == yy(ii) & stim_x == Param.x_center);
    mm = squeeze(mean(SpikeTrain(id_trl,:,:),1));
    mResp_y(ii,:,:) = conv2(mm, Smooth_kernal, 'same');
end

%%%%% plot the timecourse for different conditions
figure;
ch = [1:19 21: 32];
for ich = 1:size(mResp_x ,2)
    subplot(6,6 , ich);
    imagesc(BinT, xx, squeeze(mResp_x(:, ich,:)));
    title(['ch = ' num2str(ch(ich))]);
    axis xy;    
end
figure;    
for ich = 1:size(mResp_y ,2)
    subplot(6,6, ich);
    imagesc(BinT, yy, squeeze(mResp_y(:, ich,:)));
    title(['ch = ' num2str(ch(ich))]);
    axis xy;    
end

%%
%%%%% plot the tuning curves
close all;
figure;
X_tune = zeros(31, 2);
Y_tune = zeros(31, 2);
for ich = 1:size(mResp_x ,2)
    subplot(6,6 , ich);
    mm = mean(mResp_x(:, ich,BinT>=30 & BinT <=100),3)-mean(mResp_x(:, ich,BinT>=-100 & BinT <=0),3);
    mm = conv(mm, ones(1,3)/3, 'same');
    plot( xx, mm, 'ro'); hold on;
    if sum(abs(mm))~=0
        [fit_param1, yfit]= gaussian_fit(xx, mm);
        X_tune(ich, :) = [fit_param1.mu, fit_param1.Fitness];
        plot(xx, yfit, 'r', 'linewidth', 2.0); hold on;
    end
    
    nn = mean(mResp_y(:, ich,BinT>=30 & BinT <=100),3)-mean(mResp_y(:, ich,BinT>=-100 & BinT <=0),3);
    nn = conv(nn, ones(1,3)/3, 'same');
    plot(yy, nn, 'bo'); hold on;
    if sum(abs(nn))~=0
        [fit_param2, yfit]= gaussian_fit(yy, nn);
        Y_tune(ich, :) = [fit_param2.mu, fit_param2.Fitness];
        plot(yy, yfit, 'b', 'linewidth', 2.0); hold on;
        title(['ch = ' num2str(ch(ich)) ', ' num2str(min(fit_param1.Fitness, 100), '%5.1f') ', ' num2str(min(fit_param2.Fitness, 100), '%5.1f')]);
    end
end

%%
figure;
idx1 = X_tune(:,2)<=0.4;
chID = [1:19 21:32]; 
plot(X_tune(idx1, 1), chID(idx1), 'ro','markerfacecolor','r'); hold on;
idx2 = Y_tune(:,2)<=0.4;
plot(Y_tune(idx2, 1), chID(idx2), 'bo','markerfacecolor','b'); hold on;
xlim([-30 30])
ylim([0 32])
set(gca,'box','off');



















