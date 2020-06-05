clear all;
close all;
folder='Z:\Rujia\Results\';
monkeyset={'Mikey','Vasco'};
load('Z:\Rujia\Results\Mikey_RecordingInfo.mat');
Group={'VisResp','VisMotor','MotorResp','AllVisual','AllUnit'}; 
border=0.15;
for igroup=4
    groupName=Group{igroup};
    for imonkey=2
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end
        for iarea=1:3
            num=0;           
            for ifile=1:numel(dateUsed) %1 :numel(filename)
                date=dateUsed{ifile};
%                 load([folder 'UnSorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                load([folder 'bVisResp_' date '_SingleUnit.mat']);
                load([folder 'bMotorResp_' date '_SingleUnit.mat']);
                load([folder 'bVisMotor_' date '_SingleUnit.mat']);
                load([folder 'CueOn_' date '_SingleUnit.mat']);
                load([folder 'RFOn_' date '_SingleUnit.mat']);
                load([folder 'SaccadeOn_' date '_SingleUnit.mat']);
                if imonkey==1
                    isession=find(strcmp(RawInfo(1,:), date)==1);
                    channelID=find(ismember(SpikeTrain.channel, RawInfo{iarea+1,isession}));  % for Mikey
                else
                    channelID=find(ismember(SpikeTrain.channel,(iarea-1)*32+1:iarea*32));  % for Vasco
                end
                for ich=channelID   %1:numel(SpikeTrain.channel)  %
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
                        if size(SpikeTrain.waveform,1)==size(SpikeTrain.cueAlign,1)
                            waveshape=squeeze(mean(SpikeTrain.waveform,1));
                        else
                            waveshape=SpikeTrain.waveform;
                        end
                        
                        num=num+1;
                        [SpikeShape.ValueLeft{iarea}(num),SpikeShape.PeakLeft{iarea}(num)]=max(waveshape(ich,1:17));  % 17 is the index of valley;
                        [SpikeShape.ValueRight{iarea}(num),SpikeShape.PeakRight{iarea}(num)]=max(waveshape(ich,17:end));  % 17 is the index of valley;
                        
                        SpikeShape.WaveShape{iarea}(num,:)=waveshape(ich,:);  %SpikeTrain.waveform(ich,:);
                        SpikeShape.ElecID{iarea}(num)=SpikeTrain.channel(ich);
                        idx00=SpikeShape.ValueRight{iarea}(num)>=0.1;
                        if SpikeShape.PeakRight{iarea}(num)>=5&&SpikeShape.PeakRight{iarea}(num)<=10 &&idx00                  
                            SpikeShape.CellType{iarea}(num) = 1;  % Narrow-spike, inhibitory
                        elseif SpikeShape.PeakRight{iarea}(num)>=12&&SpikeShape.PeakRight{iarea}(num)<=25 && idx00 
                            SpikeShape.CellType{iarea}(num) = 2;  % broad-spike
                        else
                            SpikeShape.CellType{iarea}(num)=0;                            
                        end                        
                        Valley = abs(min(SpikeShape.WaveShape{iarea}(num,:)));                        
                        NormLeft=SpikeShape.ValueLeft{iarea}(num)/Valley;
                        NormRight=SpikeShape.ValueRight{iarea}(num)/Valley;
                        if NormLeft<=border && NormRight<=border
                            SpikeShape.WaveformClass{iarea}(num)=1;
                        elseif NormLeft<=border && NormRight>border
                            SpikeShape.WaveformClass{iarea}(num)=2;
                        elseif NormLeft>border && NormRight<=border
                            SpikeShape.WaveformClass{iarea}(num)=3;
                        elseif NormLeft>border && NormRight>border
                            SpikeShape.WaveformClass{iarea}(num)=4;
                        end                                                 
                    end
                end
            end
        end
        save([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat'],'SpikeShape','-v7.3');  %' groupName '_4sigma_MUA
        clear SpikeShape;
    end
end
%% plot the distribution of spike widths and amplitudes for each monkey
close all;
folder='Z:\Rujia\Results\';
Group={'VisResp','VisMotor','MotorResp','AllVisual','AllUnit'}; 
igroup=4;
groupName=Group{igroup};
monkeyset={'Mikey','Vasco'};
Ncell = cell(1,2);
for imonkey=1:2
    monkey=monkeyset{imonkey};
    load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);  %cueResponsive_
    figure;
    for iarea=1:3          
        subplot(2,3,iarea);        
%         Amplitude = transpose(SpikeShape.ValueRight{iarea})- min(SpikeShape.WaveShape{iarea},[],2);
        idx1=SpikeShape.ValueRight{iarea}>=0&SpikeShape.ValueRight{iarea}<=3;
        idx2=SpikeShape.ValueRight{iarea}<SpikeShape.ValueLeft{iarea};
        idx0=SpikeShape.PeakRight{iarea}<=30&idx2;         %&idx1
        hist(SpikeShape.PeakRight{iarea}(idx0)/40*1000,20); hold on;
        xticks(200:200:1000)
        axis tight;
        Ncell{imonkey}(iarea, 1)=sum(SpikeShape.CellType{iarea}==1&idx2);
        Ncell{imonkey}(iarea, 2)=sum(SpikeShape.CellType{iarea}==2&idx2);
%         
%         subplot(2,3,3+iarea);        
%         plot(SpikeShape.PeakRight{iarea}(idx0)/40*1000,SpikeShape.ValueRight{iarea}(idx0),'.');  hold on;  %Amplitude
% %         line([100 1000],[0.15 0.15],'color','r'); hold on;
%         axis tight;
%         xlim([100 750]);
%         ylim([0 2])
        
        subplot(2,3,3+iarea);
        yy=transpose(SpikeShape.WaveShape{iarea}(SpikeShape.PeakRight{iarea}<=10&idx1&idx2,:));
        plot(yy./repmat(abs(min(yy,[],1)),size(yy,1),1),'color','b'); hold on; 
%         plot(yy,'color','b'); hold on; 
        zz=transpose(SpikeShape.WaveShape{iarea}(SpikeShape.PeakRight{iarea}>=12&SpikeShape.PeakRight{iarea}<=30&idx1&idx2,:));
        plot(zz./repmat(abs(min(zz,[],1)),size(zz,1),1),'color','r'); hold on;
%         plot(zz,'color','r'); hold on;
        axis tight
    end
end

%% plot the distribution of spike widths and amplitudes for two monkeys
close all;
folder='Z:\Rujia\Results\';
Ncell = zeros(3,2);
Group={'VisResp','VisMotor','MotorResp','AllVisual','AllUnit'}; 
igroup=3;
groupName='AllUnit';  %Group{igroup};
S1= load([folder 'Mikey_' groupName '_SpikeShape_SU.mat']);  %cueResponsive_
S2= load([folder 'Vasco_' groupName '_SpikeShape_SU.mat']);  %cueResponsive_
figure;
for iarea=1:3
    subplot(2,3,iarea);
    AmpRight = [S1.SpikeShape.ValueRight{iarea}, S2.SpikeShape.ValueRight{iarea}];
    PeakRight = [S1.SpikeShape.PeakRight{iarea}, S2.SpikeShape.PeakRight{iarea}];
    AmpLeft = [S1.SpikeShape.ValueLeft{iarea}, S2.SpikeShape.ValueLeft{iarea}];
    PeakLeft = [S1.SpikeShape.PeakLeft{iarea}, S2.SpikeShape.PeakLeft{iarea}];
    idx1=AmpRight>=0.1&AmpRight<2;    
    idx2=AmpRight>AmpLeft;
    idx0=PeakRight<=30;  %&idx2;  %&idx1
    hist(PeakRight(idx0)/40*1000,20); hold on;
    xticks(200:200:1000)
    axis tight;
    Ncell(iarea, 1)=sum(PeakRight>=5&PeakRight<=10&idx2);
    Ncell(iarea, 2)=sum(PeakRight>=12&PeakRight<=25&idx2);
    
    subplot(2,3,3+iarea);
    plot(PeakRight(idx0)/40*1000,AmpRight(idx0),'.');  hold on;  %Amplitude
    line([100 1000],[0.3 0.3],'color','r'); hold on;
    axis tight;
    xlim([100 800]);
    ylim([0 2])
    
    %         figure;
    %         yy=transpose(SpikeShape.WaveShape{iarea}(SpikeShape.PeakRight{iarea}<=10&idx1&idx2,:));
    %         plot(yy./repmat(abs(min(yy,[],1)),size(yy,1),1),'color','b'); hold on;
    %         % plot(yy,'color','b'); hold on;
    %         zz=transpose(SpikeShape.WaveShape{iarea}(SpikeShape.PeakRight{iarea}>=12&SpikeShape.PeakRight{iarea}<=25&idx1&idx2,:));
    %         plot(zz./repmat(abs(min(zz,[],1)),size(zz,1),1),'color','r'); hold on;
    %         % plot(zz,'color','r'); hold on;
end


%% plot the distribution of narrow or broad spikes across different depths of probe
load([folder monkey '_' groupName '_SpikeShape_cueResponsive_SU.mat']);
% close all;
figure;
for iarea=1:3
idx1=SpikeShape.ValueRight{iarea}>=0;
yy=transpose(SpikeShape.ElecID{iarea}(SpikeShape.PeakRight{iarea}<=10&idx1));
subplot(3,2,(iarea-1)*2+1);
hist(yy);
xlim([iarea*32-31 iarea*32]);

zz=transpose(SpikeShape.ElecID{iarea}(SpikeShape.PeakRight{iarea}>=12&SpikeShape.PeakRight{iarea}<=25&idx1));
subplot(3,2,(iarea-1)*2+2);
hist(zz);
xlim([iarea*32-31 iarea*32]);
end

%% scatter plot the right amplitude to the left amplitude
close all;
folder='Z:\Rujia\Results\';
Ncell = zeros(3,4);
Group={'VisResp','VisMotor','MotorResp','AllVisual','AllUnit'}; 
AreaName = {'LIP', 'PUL', 'FEF'};
igroup=4;
groupName=Group{igroup};
S1= load([folder 'Mikey_' groupName '_SpikeShape_cueResponsive_SU.mat']);  %cueResponsive_
S2= load([folder 'Vasco_' groupName '_SpikeShape_cueResponsive_SU.mat']);  %cueResponsive_
figure;
clr4 = colormap(hsv(4));
border = 0.15;
for iarea=1:3    
    AmpRight = [S1.SpikeShape.ValueRight{iarea}, S2.SpikeShape.ValueRight{iarea}];
    PeakRight = [S1.SpikeShape.PeakRight{iarea}, S2.SpikeShape.PeakRight{iarea}];
    AmpLeft = [S1.SpikeShape.ValueLeft{iarea}, S2.SpikeShape.ValueLeft{iarea}];
    PeakLeft = [S1.SpikeShape.PeakLeft{iarea}, S2.SpikeShape.PeakLeft{iarea}];
    SpikeWave = [S1.SpikeShape.WaveShape{iarea}; S2.SpikeShape.WaveShape{iarea}];
    Valley = abs(min(SpikeWave, [],2));
    NormWave = SpikeWave./repmat(Valley, 1, size(SpikeWave,2));
    NormLeft=AmpLeft./Valley';
    NormRight=AmpRight./Valley';    
    idx0=PeakRight<=25;
    idx1 = NormLeft<=border&NormRight<=border&idx0;
    idx2 = NormLeft<=border&NormRight>border&idx0;
    idx3 = NormLeft>border&NormRight<=border&idx0;
    idx4 = NormLeft>border&NormRight>border&idx0;    
%     subplot(2,2,1)
%     plot(transpose(NormWave(idx1,:)), 'color', clr4(1,:)); hold on;    
%     subplot(2,2,2)
%     plot(transpose(NormWave(idx2,:)), 'color', clr4(2,:)); hold on;    
%     subplot(2,2,3)
%     plot(transpose(NormWave(idx3,:)), 'color', clr4(3,:)); hold on;    
%     subplot(2,2,4)
%     plot(transpose(NormWave(idx4,:)), 'color', clr4(4,:)); hold on;
    
%     subplot(3,4,iarea*4-3);
%     Ncell(iarea,1)=sum(idx1);
% %     plot(transpose(NormWave(idx1,:)), 'color', clr4(1,:)); hold on;
%     hist(PeakRight(idx1)/40*1000,20); hold on;
%     h = findobj(gca,'Type','patch');
%     h.FaceColor = clr4(1,:);
%     title([AreaName{iarea} ':' num2str(Ncell(iarea,1))]);
%     axis tight;
% %     ylim([-1 1.2]);   
%     subplot(3,4,iarea*4-2);
%     Ncell(iarea,2)=sum(idx2);
% %     plot(transpose(NormWave(idx2,:)), 'color', clr4(2,:)); hold on;
%     hist(PeakRight(idx2)/40*1000,20); hold on;
%     h = findobj(gca,'Type','patch');
%     h.FaceColor = clr4(2,:);
%     title([AreaName{iarea} ':' num2str(Ncell(iarea,2))]);axis tight;
% %     ylim([-1 1.2]);    
%     subplot(3,4,iarea*4-1);
%     Ncell(iarea,3)=sum(idx3);
% %     plot(transpose(NormWave(idx3,:)), 'color', clr4(3,:)); hold on;
%     hist(PeakRight(idx3)/40*1000,20); hold on;
%     h = findobj(gca,'Type','patch');
%     h.FaceColor = clr4(3,:);
%     title([AreaName{iarea} ':' num2str(Ncell(iarea,3))]);axis tight;
% %     ylim([-1 1.2]);    
%     subplot(3,4,iarea*4);
%     Ncell(iarea,4)=sum(idx4);
% %     plot(transpose(NormWave(idx4,:)), 'color', clr4(4,:)); hold on;
%     hist(PeakRight(idx4)/40*1000,20); hold on;
%     h = findobj(gca,'Type','patch');
%     h.FaceColor = clr4(4,:);
%     title([AreaName{iarea} ':' num2str(Ncell(iarea,4))]);axis tight;
% %     ylim([-1 1.2]);
    
    %%%%% plot the right peak amplitude against left peak
    subplot(2,3,iarea);
    plot(AmpRight(idx0), AmpLeft(idx0), 'r.'); hold on;
    line([border, border], [-.1, 1.3], 'color', 'k'); hold on;
    line([-.1, 1.3],[border, border], 'color', 'k'); hold on;
    axis tight
    xlabel('Right Peak');
    ylabel('Left Peak');
    title(AreaName{iarea});     
    subplot(2,3,3+iarea);
%     subplot(3,1,iarea);
    plot(NormRight(idx1), NormLeft(idx1), '.', 'color', clr4(1,:));  hold on; 
    plot(NormRight(idx2), NormLeft(idx2), '.', 'color', clr4(2,:));  hold on; 
    plot(NormRight(idx3), NormLeft(idx3), '.', 'color', clr4(3,:));  hold on; 
    plot(NormRight(idx4), NormLeft(idx4), '.', 'color', clr4(4,:));  hold on; 
    line([border, border], [-0.1 1.3], 'color', 'k'); hold on;
    line([-0.1 1.3],[border, border], 'color', 'k'); hold on;
    xlabel('Normalized Right Peak');
    ylabel('Normalized Left Peak');
    axis tight;    
end

