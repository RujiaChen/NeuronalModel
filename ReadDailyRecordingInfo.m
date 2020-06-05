%%% read the recording information in each session for Mikey
clear all;
folder='Y:\anne\monkeys\monkey forms\';
filename='daily recording sheets.xlsx';
[num, text, raw]=xlsread([folder filename]);
isession=0;
for jcell=3:size(text,2)
    if ~isempty(strfind(text{1, jcell}, '/201'))
        isession=isession+1;
        tt=text{1, jcell};
        idx11=strfind(tt,'/');
        if idx11(1)==2
            tt=['0' tt];
        end
        idx00=strfind(tt,'/');
        if idx00(2)<6
            tt=[tt(1:idx00(1)) '0' tt(idx00(1)+1:end)]; 
        end
        RawInfo{1,isession}=strrep(tt, '/', '');
        RawInfo{1,isession}=strrep(RawInfo{1,isession}, '2017', '17');
        RawInfo{1,isession}=strrep(RawInfo{1,isession}, '2018', '18');

        Depth{1, isession}=RawInfo{1,isession};
                
        for kk=1:3
            pp=text{kk*3+2, jcell};
            idx=strfind(pp, '-');
            if isempty(idx)
                RawInfo{kk+1,isession}=str2double(pp);
            else
                elec1=str2double(pp(1:idx-1));
                elec2=str2double(pp(idx+1:end));
                RawInfo{kk+1,isession}=elec1:elec2;
            end
            
            Depth{kk+1, isession}=str2double(text{kk*3+3, jcell});
        end
    end
end
save('Z:\RujiaChen\Results\Mikey_RecordingDepth.mat', 'Depth', '-v7.3');
save('Z:\RujiaChen\Results\Mikey_RecordingInfo.mat', 'RawInfo', '-v7.3');

%%
date='122317';
xx=cellfun(@(x) strcmp(x, date), RawInfo(1,:));

yy=find(strcmp(RawInfo(1,:), date)==1);