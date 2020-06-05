close all;

dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey 
load('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat');

% channel=1:96;  % allArea for Vasco
areaName='allArea';
% folder='Y:\anne\monkeys\plexondata\';  % for vasco
folder='Y:\anne\monkeys\pigz_plexondata\mikey_left_hemisphere\';  % for Mikey
for idate=1% : numel(dateUsed)
date=dateUsed{idate};
% if ~exist(['Z:\RujiaChen\Results\Unsorted_SpikeTrain_' areaName  '_' date '_TriggerCorrected.mat'],'file')
fileName=[folder date '\merged_spl_noWB_001-01.pl2'];
muaDataDirRoot=folder;   %['Z:\RujiaChen\Vasco\' date '\'];   %042118
sessionName='merged_spl_noWB_001';
isession=find(strcmp(RawInfo(1,:), date)==1);
channel=[RawInfo{2,isession} RawInfo{3,isession} RawInfo{4,isession}];  % 1:96   
% channel=1:96;
isLoadSpikes=1;
isLoadMua=0;
isLoadLfp=0;
isLoadSpkc=0;
isLoadDirect=0;
spikeChannelPrefix = 'SPK_SPKC';
% spikeChannelPrefix = 'SPK';
spikeChannelsToLoad=channel;
muaChannelsToLoad=channel;
lfpChannelsToLoad = channel; 
spkcChannelsToLoad = channel;
directChannelsToLoad =3;
D = loadPL2(fileName, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad);
    
    
%%
for icell=1:numel(D.allSpikeStructs)
    timepoint=D.allSpikeStructs{icell}.ts;
    wave=D.allSpikeStructs{icell}.wf;

    
end