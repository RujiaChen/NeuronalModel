close all;
clear all;
monkeyset={'Vasco', 'Mikey'};
areaName={'LIP', 'PUL', 'FEF'};
imonkey=1;
monkey=monkeyset{imonkey}; 

iarea=3;
Area=areaName{iarea};
% load(['Z:\RujiaChen\Results\Resp_MStoCue_' monkey '_' Area '.mat']);
% load(['Z:\RujiaChen\Results\Resp_MSawayCue_' monkey '_' Area '.mat']);
% TT=-199:199;
load(['Z:\RujiaChen\Results\LFP_MStoCue_CueAligned_' monkey '_' Area '.mat']);
load(['Z:\RujiaChen\Results\LFP_MSawayCue_CueAligned_' monkey '_' Area '.mat']);
TT=-400:600;

if strcmp(monkey, 'Vasco')
    clr2=[0 0 1; 1 0 0];
elseif strcmp(monkey, 'Mikey')
    clr2=[ 1 0 0; 0 0 1];
end

NN=zeros(1, 2);
MM=zeros(1, 2);
yrange=[-2 1];
legendName={'Attend to','Attend away'};
figure;
for icue=1:2  %:2
%     figure;
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
    
    %%%%% plot the mean and ste of the LFP signals
    subplot(1,2,1);
    mm1=mean(RespTo,1);
    ee=std(RespTo, [], 1)/sqrt(size(RespTo, 1));
    patch([TT fliplr(TT) TT(1)], [mm1-ee fliplr(mm1+ee) mm1(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
    fig1(icue)=plot(TT, mm1, 'color', clr2(icue,:), 'linewidth', 1.5); hold on;
    if strcmp(monkey, 'Vasco')
        line([0 0], yrange);
        ylim(yrange);
        ll=legend(fig1, 'Attend away', 'Attend to');
        set(ll,'box', 'off', 'location', 'best');
    else
        line([0 0], yrange);
        ylim(yrange);
        ll=legend(fig1, 'Attend to', 'Attend away');
        set(ll,'box', 'off', 'location', 'best');
    end
    xlim([TT(1) TT(end)]);
    xlabel('Time to cue onset (ms)');  %MS
    ylabel('Evoked potential');
    set(gca, 'box' ,'off');
    
    subplot(1,2,2);
    mm2=mean(RespAway,1);
    ee=std(RespAway, [], 1)/sqrt(size(RespAway, 1));
    patch([TT fliplr(TT) TT(1)], [mm2-ee fliplr(mm2+ee) mm2(1)-ee(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
    plot(TT, mm2, 'color', clr2(icue,:), 'linewidth', 1.5); hold on;
    if strcmp(monkey, 'Vasco')       
        line([0 0], yrange);
        ylim(yrange);
        ll=legend(fig1, 'Attend away', 'Attend to');
        set(ll,'box', 'off', 'location', 'best');
    else
        line([0 0], yrange);
        ylim(yrange);
        ll=legend(fig1, 'Attend to', 'Attend away');
        set(ll,'box', 'off', 'location', 'best');
    end
    xlim([TT(1) TT(end)]);
    xlabel('Time to cue onset (ms)');  %MS
    ylabel('Evoked potential');
    set(gca, 'box' ,'off');
    
%     %%%%%% plot the power spectrum of LFP 
%     params.tapers=[2,3];
%     params.pad=0;
%     params.Fs=1000;
%     params.fpass=[0 100];
%     params.err=[1 0.05];
%     params.trialave=1;
%     movingwin=[0.2, 0.02];   %for multi-taper functions
%     
%     %%%% plot the power spectrum
%     subplot(1,2,1);
%     freq=3:2:50;
%     Amp=zeros(1, length(freq));
%     steAmp=zeros(1, length(freq));
%     num=0;
%     for ifreq=3:2:50 
%         num=num+1;
%         signal_filtered=bandpass(RespTo', [ifreq-2 ifreq+2], params.Fs);
%         signal_hilbert=hilbert(signal_filtered);
%         aa=max(abs(signal_hilbert),[],1);
%         Amp(num)=mean(aa);
%         steAmp(num)=std(aa)/sqrt(numel(aa));        
%     end 
%     xx=freq;
%     yy1=Amp-steAmp;
%     yy2=Amp+steAmp;
%     patch([xx fliplr(xx) xx(1)], [yy1 fliplr(yy2) yy1(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     plot(xx, Amp,'color', clr2(icue,:)); hold on;
%     xlabel('Frequency (Hz)');
%     ylabel('Power spectrum');
%     
%     subplot(1,2,2);
%     freq=3:2:50;
%     Amp=zeros(1, length(freq));
%     steAmp=zeros(1, length(freq));
%     num=0;
%     for ifreq=3:2:50 
%         num=num+1;
%         signal_filtered=bandpass(RespAway', [ifreq-2 ifreq+2], params.Fs);
%         signal_hilbert=hilbert(signal_filtered);
%         aa=mean(abs(signal_hilbert),1);
%         Amp(num)=mean(aa);
%         steAmp(num)=std(aa)/sqrt(numel(aa));        
%     end 
%     xx=freq;
%     yy1=Amp-steAmp;
%     yy2=Amp+steAmp;
%     patch([xx fliplr(xx) xx(1)], [yy1 fliplr(yy2) yy1(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     plot(xx, Amp,'color', clr2(icue,:)); hold on;
%     xlabel('Frequency (Hz)');
%     ylabel('Power spectrum');    
    
%     %%%%%% plot the power and phase of the peak frequency
%     if strcmp(monkey, 'Vasco')
%         ifreq=15;
%     else
%         ifreq=19;
%     end
%     xx=TT;
%     signal_filtered=bandpass(RespTo', [ifreq-2 ifreq+2], params.Fs);
%     signal_hilbert=hilbert(signal_filtered);
%     
%     %%%%% plot an example trials
% %     subplot(1,3,1);
% %     plot(TT, real(signal_hilbert)); hold on;
% %     line([0 0],[-0.5 0.5]); hold on;
% %     title('EP');
% %     subplot(1,3,2);
% %     plot(TT,angle(signal_hilbert)); hold on;
% %     line([0 0],[-4 4]); hold on;
% %     title('Phase');
% %     subplot(1,3,3);
% %     plot(TT, abs(signal_hilbert)); hold on;
% %     line([0 0],[0 0.5]); hold on;
% %     title('Amplitude');
% %%%%%%%%%end of example plot
% 
%     subplot(1,2,1);
%     data=abs(signal_hilbert);   % signal_filtered
%     Amp=mean(data,2);
%     steAmp=std(data, [],2)/sqrt(size(data,2));
%     Ang=angle(signal_hilbert);
% %     cc=-pi:pi/10:pi;
% %     histogram(Ang(xx==0,:),cc);
% %     phase_cos=mean(cos(Ang),2);
% %     plot(TT, phase_cos,'color', 'g'); hold on; %clr2(icue,:)
% 
% 
%     yy1=transpose(Amp-steAmp);
%     yy2=transpose(Amp+steAmp);
%     patch([xx fliplr(xx) xx(1)], [yy1 fliplr(yy2) yy1(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     fig1(icue)=plot(xx, Amp,'color', clr2(icue,:)); hold on;
%     if strcmp(monkey, 'Vasco')
% %         ll=legend(fig1, 'Attend away', 'Attend to');
%         line([0 0], yrange);
%         ylim(yrange);
%     else
% %         ll=legend(fig1, 'Attend to', 'Attend away');
%         line([0 0], yrange);
%         ylim(yrange);
%     end
% %     set(ll,'box', 'off', 'location', 'best');
%     xlabel('Time to MS (ms)');
% %     ylabel('Evoked potential');
%     ylabel('Amplitude');
%     
%     subplot(1,2,2);
%     signal_filtered=bandpass(RespAway', [ifreq-2 ifreq+2], params.Fs);
%     signal_hilbert=hilbert(signal_filtered);
%     data=abs(signal_hilbert);   % signal_filtered
%     Amp=mean(data,2);
%     steAmp=std(data, [],2)/sqrt(size(data,2));
%     Ang=angle(signal_hilbert);
% %     histogram(Ang(xx==0,:),cc);
% %     phase_cos=mean(cos(Ang),2);
% %     plot(TT, phase_cos,'color', 'k'); hold on;  %clr2(icue,:)
% 
%      %%%%%% plot the averaged LFP amplitude
%     yy1=transpose(Amp-steAmp);
%     yy2=transpose(Amp+steAmp);
%     patch([xx fliplr(xx) xx(1)], [yy1 fliplr(yy2) yy1(1)],  clr2(icue,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w'); hold on;
%     fig2(icue)=plot(xx, Amp,'color', clr2(icue,:)); hold on;
%     if strcmp(monkey, 'Vasco')
% %         ll=legend(fig2, 'Attend away', 'Attend to');
%         line([0 0], yrange);
%         ylim(yrange);
%     else
% %         ll=legend(fig2, 'Attend to', 'Attend away');
%         line([0 0], yrange);
%         ylim(yrange);
%     end
% %     set(ll,'box', 'off', 'location', 'best');
%     xlabel('Time to MS (ms)');
% %     ylabel('Evoked potential');
%     ylabel('Amplitude');
end