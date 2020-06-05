folder='Z:\anne\monkeys\pigz_plexondata\mikey_left_hemisphere\';
subfolder=dir(folder);
TargetFolder='C:\Rujia\Mikey\Plexon\RawData\';
for ifolder=5:numel(subfolder)
    if exist([TargetFolder subfolder(ifolder).name],'dir')==0
        mkdir([TargetFolder subfolder(ifolder).name]);       
    end
    fileName=dir([folder subfolder(ifolder).name '\merged*.pl2']);
    if ~isempty(fileName)
        
        for ii=1:numel(fileName)
            if exist([TargetFolder subfolder(ifolder).name '\' fileName(ii).name],'file')==0
                copyfile([fileName(ii).folder '\' fileName(ii).name],[TargetFolder subfolder(ifolder).name '\' fileName(ii).name]);
            end
        end 
    end   
end