%% get the d-prime for different groups 
clear all; 
close all;
folder='Z:\Rujia\Results\';
% dateUsed={'092718','092818','100418','100518','100718', '100818','101018','101118','101218','101518','101618', '101718', '102218','102318','102418','102518', '103118', '110118', '111618'};  % for Vasco
monkeyset={'Mikey','Vasco'};
TargetPosY=[200,-200];
cueLoc=[2,5];  %cue-on and off
RFtype=[1,2]; % different target shape
cueTyp=[0 1];
cuetype={'endo','exo'};
load('Z:\Rujia\Results\Mikey_RecordingInfo.mat');
Group={'VisResp','VisMotor','MotorResp','AllVisual'};    
Time=-600:1200; 
Time1=-505:1305;
window=[0.05 0.025];
Fs=1000;
for igroup=4 %1:3
    groupName=Group{igroup};
    for imonkey=1 
        monkey=monkeyset{imonkey};
        if strcmp(monkey, 'Mikey')
            dateUsed={ '111817','112017','112217', '112417', '112617', '112817','112917', '120117', '120217', '120317', '121817','121917', '122217', '122317'};  % for Mikey
        elseif strcmp(monkey, 'Vasco')
            dateUsed={'042418','051218','050118','050418','050518','042118','051718','112018',  '010719', '010919', '011119', '011419', '011519', '011719', '012119'};  % for Vasco  '112118',  '011819',
        end        
        for isexo=0:1
            DprimeCue=cell(3,1);
            DprimeArray=cell(3,1);            
%             ElecInfo=cell(3,2);
            DprimeCueDelay=cell(3,1);
            DprimeArrayDelay=cell(3,1);
            pValueCueDelay=cell(3,1);
            pValueArrayDelay=cell(3,1);
            for iarea=1:3
                num=0;
                for ifile = 1:numel(dateUsed) %1 :numel(filename)
                    date=dateUsed{ifile};
                    load([folder 'Sorted_SpikeTrain_allArea_' date '_TriggerCorrected.mat']);
                    load([folder 'CorrectTrialParam_' date '.mat']);
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
                    for ich=channelID
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
                            num=num+1;
                            %%%%%% prime position
                            if RFOn(ich)>0
                                ipos=ceil(RFOn(ich)/2);
                                icue=2-mod(RFOn(ich),2);
                            elseif SaccadeOn(ich)>0
                                ipos=ceil(SaccadeOn(ich)/2);
                                icue=2-mod(SaccadeOn(ich),2);
                            else
                                ipos=1;
                                icue=imonkey;
                            end
                            
                            idx1=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            idxCueT=Time>=-300&Time<=1000;
                            idxArrayT=Time1>=-300&Time1<=1000;
                            CueOnResp=squeeze(SpikeTrain.cueAlign(idx1,ich,idxCueT)); 
                            ArrayOnResp=squeeze(SpikeTrain.arrayAlign(idx1,ich,idxArrayT));  
                            BaseOn=repmat(squeeze(mean(SpikeTrain.cueAlign(idx1,ich,Time>=-300&Time<=0),3)),1,size(CueOnResp,2));
                            CueOnResp=CueOnResp-BaseOn;
                            ArrayOnResp=ArrayOnResp-BaseOn;
                             
                            idx2=CorrectTrialParam.rfyi==TargetPosY(ipos)&CorrectTrialParam.cueloc==cueLoc(3-icue)&CorrectTrialParam.isexocue==isexo &CorrectTrialParam.is_mouse_trial==0;
                            CueOffResp=squeeze(SpikeTrain.cueAlign(idx2,ich,idxCueT));
                            ArrayOffResp=squeeze(SpikeTrain.arrayAlign(idx2,ich,idxArrayT));
                            BaseOff=repmat(squeeze(mean(SpikeTrain.cueAlign(idx2,ich,Time>=-300&Time<=0),3)),1,size(CueOffResp,2));
                            CueOffResp=CueOffResp-BaseOff;
                            ArrayOffResp=ArrayOffResp-BaseOff;
                            
%                             [DprimeCue{iarea}(num,:),~,~]=getDPrime(CueOnResp,CueOffResp,window,Fs,0);
%                             [DprimeArray{iarea}(num,:),~,T]=getDPrime(ArrayOnResp,ArrayOffResp,window,Fs,0);

%                             ElecInfo{iarea,1}(num,:)=[SpikeTrain.channel(ich), ifile, RFOn(ich)];
                            
                            [DprimeCueDelay{iarea}(num),pValueCueDelay{iarea}(num)]=getDPrime(mean(CueOnResp(:,550:800),2),mean(CueOffResp(:,550:800),2));
                            [DprimeArrayDelay{iarea}(num),pValueArrayDelay{iarea}(num)]=getDPrime(mean(ArrayOnResp(:,550:800),2),mean(ArrayOffResp(:,550:800),2));
                        end
                    end
                end
            end
            
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCue_cueResponsive_SU.mat'],'DprimeCue','-v7.3');
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArray_cueResponsive_SU.mat'],'DprimeArray','-v7.3'); 
%             save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeT.mat'],'T','-v7.3');
            
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeCueDelay_250-500_cueResponsive_SU.mat'],'DprimeCueDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_DprimeArrayDelay_250-500_cueResponsive_SU.mat'],'DprimeArrayDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueCueDelay_250-500_cueResponsive_SU.mat'],'pValueCueDelay','-v7.3');
            save([folder monkey '_' cuetype{isexo+1} '_' groupName '_pValueArrayDelay_250-500_cueResponsive_SU.mat'],'pValueArrayDelay','-v7.3');
        end
    end
end
