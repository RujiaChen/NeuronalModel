folder='Y:\anne\monkeys\Mikey\logfiles\';
dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey 
delayCue=zeros(numel(dateUsed), 2);
delayArray=zeros(numel(dateUsed), 2);
for idate=1:numel(dateUsed)
    date=dateUsed{idate};
    logDate=['20' date(5:6) date(1:4)];
    % logFolder=['Y:\anne\monkeys\presentation files\Vasco\logfiles\' logDate '\'];
    logFolder=['Y:\anne\monkeys\Mikey\logfiles\' logDate '\'];
    filename=dir([logFolder 'flanker_params*.txt']);
    
    if ~isempty(filename)
        logPresentation=[logFolder filename(1).name];
        fileID=fopen(logPresentation,'r');
        tline=fgetl(fileID);
        tline=fgetl(fileID);
        Info=textscan(tline, '%d');
        
        delayCue(idate, :)=[Info{1}(6)  Info{1}(7)];
        delayArray(idate, :)=[Info{1}(8)  Info{1}(9)];
        
    end
end