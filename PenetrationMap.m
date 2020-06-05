close all;
folder='C:\Rujia\Vasco\';
fileID=fopen([folder 'RecordingDepth_Location.csv'],'r');
tline=fgets(fileID);
headInfo=textscan(tline,'%s','delimiter',',');
formatSpec='%s %s %f %f %f';
dataRaw=textscan(fileID,formatSpec,'delimiter',',');
fclose(fileID);
AreaName={'LIP','PUL','FEF'};
dateUsed={'051018','042418','050718','051218','050118','050418','050518','042118','051718','041718'};
figure;
stats=zeros(3, 4);
params=zeros(3);
for iarea=1:3
    xx=dataRaw{4}(iarea:3:end);
    yy=dataRaw{5}(iarea:3:end);
    depth=dataRaw{3}(iarea:3:end);
    xloc=-3:0.5:2.5;
    yloc=-3:0.5:2;
    zz=zeros(length(yloc), length(xloc));
    
    SessionName=dataRaw{1}([iarea:3:end]);
    for ifile=1:numel(SessionName)
        ll=length(SessionName{ifile});
        if ll==4
            tt=['0' SessionName{ifile}(1) SessionName{ifile}(3:4) '18'];
        elseif ll==3
            tt=['0' SessionName{ifile}(1) '0' SessionName{ifile}(3) '18'];
        end
        ss{ifile}=tt;
    end
    
    for ii=1:numel(yloc)
        for jj=1:numel(xloc)
            idx=find(yy==yloc(ii)&xx==xloc(jj));
            if ~isempty(idx)&&ismember(ss{idx}, dateUsed)
                zz(ii, jj)=depth(idx);
            end
        end
    end
    subplot(1,3,iarea);
    imagesc(xloc,yloc,zz);
    caxis([min(depth)-1, max(depth)]);
    axis xy;
    colorbar;
    title(AreaName{iarea});

% % % % % regress the relationship between location and depth;
X=[ones(size(xx)) xx yy];
[b,bint,r,rint,stats(iarea,:)] = regress(depth,X);
scatter3(xx, yy, depth, 'filled');  hold on;
x1fit = min(xx):0.5:max(xx);
x2fit = min(yy):0.5:max(yy);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT ;
mesh(X1FIT,X2FIT,YFIT)
title(AreaName{iarea});
params(iarea,:)=b;
end

save('Z:\RujiaChen\Results\Regress_location_depth_params.mat', 'params', '-v7.3');

%%

for ii=1:numel(dataRaw{1})
    ll=length(dataRaw{1}{ii});
    if ll==4
        tt=['0' dataRaw{1}{ii}(1) dataRaw{1}{ii}(3:4) '18'];
    elseif ll==3
        tt=['0' dataRaw{1}{ii}(1) '0' dataRaw{1}{ii}(3) '18'];
    end    
    dataRaw{1}{ii}=tt;
end

save('Z:\RujiaChen\Results\RecordingDepth.mat','dataRaw','-v7.3');

