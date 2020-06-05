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
Bar_ID =  cellfun(@(x) contains(x, '_sac'), Files_merged(:,1));

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
log_files = dir([folder date '\' date_log '*delsacc*.txt']);
log_files_cell = transpose(struct2cell(log_files));
idx = cellfun(@(x) ~contains(x, 'EYE'), log_files_cell(:,1));
Param_log = log_files_cell(idx,:);

logfile_time_offset = cellfun(@(x) sum(abs(datenum(x) - datenum(Plexon_session_time))), Param_log(:,3));
[~, log_file_ID] = min(logfile_time_offset);

logPresentation = [folder date '\' Param_log{log_file_ID,1}];
fileID=fopen(logPresentation,'r');
tline=fgets(fileID);
headInfo=textscan(tline,'%s','delimiter','\t');
headInfo_new = headInfo{1}(~cellfun(@isempty, headInfo{1}));

formatSpec='%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s';
startRow=1;

% dataRaw=textscan(fileID,formatSpec);
TrlInfo=textscan(fileID, formatSpec,  'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
idx = ~isnan(TrlInfo{2});
TrlInfo_new= cell(1, numel(TrlInfo));
for icell = 1: numel(TrlInfo)
    TrlInfo_new{icell} = TrlInfo{icell}(idx);
end

for ii = 1:numel(headInfo_new) 
    Param.(headInfo_new{ii})=TrlInfo_new{ii};
end

%%
clear SpikeTrain;
Trial_plexon = sum(EventCode==1); % event=1 is trial onset
Success_trials = sum(Param.pass_sacc_resp>0);
% assert(Trial_plexon == max(Param.trial_attempt));
assert(sum(EventCode==64) == Success_trials);
correct_idx = find(Param.pass_sacc_resp>0);
for ii = 1:numel(headInfo_new) 
    Correct_param.(headInfo_new{ii})=Param.(headInfo_new{ii})(correct_idx);
end

Reward_time = allEventTimes(EventCode==64);
N = length(Reward_time);
CueOn=zeros(1, N);
SaccOn =zeros(1, N); 
for itrl = 1:N
    idx1 = find(EventCode==2 & allEventTimes < Reward_time(itrl),1,'last' );
    CueOn(itrl) = allEventTimes(idx1)/10;
    idx2 = find(EventCode==16 & allEventTimes < Reward_time(itrl),1,'last' );
    SaccOn(itrl) = allEventTimes(idx2)/10;
    
    Duration=[CueOn(itrl)-100 CueOn(itrl)+300];  %% the param is recorded by trials
    Dura_2 = [SaccOn(itrl)-300 SaccOn(itrl)+600];
    for ich=1:numel(D.allSpikeStructs)
        idx00=find(D.allSpikeStructs{ich}.ts*1000>=Duration(1)& D.allSpikeStructs{ich}.ts*1000<=Duration(2));
        if ~isempty(idx00)
            spkTime = D.allSpikeStructs{ich}.ts(idx00)*1000-CueOn(itrl);
            BinT=-100:10:300; %
            SpikeTrain.CueOn(itrl,ich,:)=hist(spkTime,BinT);
        end
        
        idx11=find(D.allSpikeStructs{ich}.ts*1000 >= Dura_2(1)& D.allSpikeStructs{ich}.ts*1000 <= Dura_2(2));
        if ~isempty(idx11)
            spkTime = D.allSpikeStructs{ich}.ts(idx11)*1000- SaccOn(itrl);
            BinT=-300:10:600; 
            SpikeTrain.SaccOn(itrl,ich,:)=hist(spkTime,BinT);
        end
    end 
end

%%
chs= [1:19 21:32];
t1= -100:10:300;
t2= -300:10:600;
mm = squeeze(mean(SpikeTrain.CueOn,1));
nn = squeeze(mean(SpikeTrain.SaccOn,1));
figure;
subplot(1,2,1);
imagesc(t1, chs, mm);
subplot(1,2,2)
imagesc(t2, chs, nn)




%%

theta = atan(Correct_param.cue_x./Correct_param.cue_y);
idx1 = Correct_param.cue_x<0;
theta(idx1) = theta(idx1)+pi;
theta_all = unique(theta);
t1= -100:10:300;
t2= -300:10:600;
num=0;
CueResp = zeros(length(theta_all), 2);
SaccResp = zeros(length(theta_all), 2);
for ii =1:numel(theta_all)    
    idx = find(theta == theta_all(ii));
    CueResp(ii,:) = squeeze(mean(mean(SpikeTrain.CueOn(idx,:,t1>30&t1<=100),1),3))-squeeze(mean(mean(SpikeTrain.CueOn(idx,:,t1>-100&t1<=0),1),3));
    SaccResp(ii,:) = squeeze(mean(mean(SpikeTrain.SaccOn(idx,:,t1>30&t1<=100),1),3))-squeeze(mean(mean(SpikeTrain.SaccOn(idx,:,t1>-100&t1<=0),1),3));
end
























