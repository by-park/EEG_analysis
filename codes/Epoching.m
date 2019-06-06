%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG Epoching
%
%   -> Epoching (baseline)
%
% created at 2019.05.12 PBY
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% parameters
% Name of the analysis
AnalyName = 'ep';

% EEG file location
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\뇌파데이터\pp_ica_removed\'; 
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex2_n2pc_dist\뇌파데이터\pp_ica_removed\';
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\뇌파데이터\pp_ica_removed\';

cd(EEGloc);

% behavior file location
% behavfile = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\행동데이터\raw16_matlab_cueloca.xlsx';
% behavfile = 'D:\EEG 파일 모음\내 실험\Ex2_n2pc_dist\행동데이터\raw16_matlab_cueloca.xlsx';
behavfile = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\행동데이터\raw16_matlab_cueloca.xlsx';


% find Files
runStr = '*.set';
fileNames = dir(runStr);
foldersize = size(fileNames,1);

% number of trials
trialNum = 576; % exp7 %1024; % exp1, exp2

% epoching range
grange1 = -100; % -0.1
grange1 = grange1 * 0.001;
grange2 = 500; % 0.5
grange2 = grange2 * 0.001;

% outlier limit
outlier_limit = 200;

%% load behav file (matlab "load data menu = 데이터 가져오기 메뉴")
[~, ~, raw] = xlsread(behavfile,'Sheet1');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기
% behavior file
behavFile = reshape([raw{:}],size(raw));
clearvars raw R;

behavTrial = 1;

%% eeglab
eeglab;

%% preprocessing
for file = 1:foldersize
        
    fprintf('\n\nprocessing file %d \n', file)  
    
    % 1. load EEG file
%     EEG = pop_loadcnt(fileNames(file).name , 'dataformat', 'auto', 'memmapfile', '');
%     EEG = eeg_checkset( EEG );
    EEG = pop_loadset('filename',fileNames(file).name,'filepath',EEGloc);
    EEG = eeg_checkset( EEG );
    
    % 2-1. Epoch (-100ms baseline 500ms epoch)
%    EEG = pop_epoch( EEG, {  '11'  '12'  '13'  '14'  '21'  '22'  '23'  '24'  '91'  '92'  '93'  '94'  }, [grange1 grange2], 'newname', 'CNT file resampled pruned with ICA epochs', 'epochinfo', 'yes');
    EEG = pop_epoch( EEG, {  '11'  '12'  '13'  '21'  '22'  '23'  '31'  '32'  '33'  '41'  '42'  '43'  }, [grange1 grange2], 'newname', 'CNT file resampled pruned with ICA epochs', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    
    % 2-2. Remove baseline (-100 ms)
    EEG = pop_rmbase( EEG, [grange1*1000 0]);
    EEG = eeg_checkset( EEG );    
    
    % 3. exclude practice trials (only 1024 trials are MAIN EXPERIMENT!!)
    EEG.event = EEG.event(size(EEG.event,2)-trialNum*2+1:size(EEG.event,2)); % 2 triggers in each trial
    EEG.data = EEG.data(:,:,size(EEG.data,3)-trialNum+1:size(EEG.data,3)); 
    EEG.epoch = EEG.epoch(size(EEG.epoch,2)-trialNum+1:size(EEG.epoch,2)); 
        
    % 4. check outlier (do not use eeglab because of behavior files)
    for trialEEG = 1 : 2 : size(EEG.event,2) % 2048 (cue trigger 1024 + target trigger 1024)
        if sum(sum(EEG.data(:,:,(trialEEG+1)/2) > outlier_limit))>0 || sum(sum(EEG.data(:,:,(trialEEG+1)/2) < -outlier_limit)) > 0
            EEG.event(trialEEG).type = nan; % triger1
            EEG.event(trialEEG+1).type = nan; % triger2
            EEG.data(:,:,(trialEEG+1)/2) = nan;
            behavFile(behavTrial,:) = [];
        else
            behavTrial = behavTrial + 1; % pass a behavior trial
        end 
    end
    
    % 5. remove checked list
    ind = 1;
    while ind ~= size(EEG.event, 2) + 1
        if isnan(EEG.event(ind).type)
            EEG.event(ind) = []; % triger1 
            EEG.event(ind) = []; % triger2
            EEG.data(:,:,(ind+1)/2) = [];
            EEG.epoch((ind+1)/2) = [];
        else
           ind = ind + 2; % 2 trigers
        end     
    end

    for nevent = 1:2:size(EEG.data,3)*2
        EEG.event(nevent).epoch = (nevent+1)/2;
        EEG.event(nevent+1).epoch = (nevent+1)/2;
        EEG.epoch((nevent+1)/2).event = [nevent nevent+1];
    end
    
    EEG = eeg_checkset( EEG );
    
    % 6. save
    if ~exist(AnalyName, 'dir')
        mkdir(AnalyName);
    end
    EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(fileNames(file).name, "."), '_', AnalyName, '.set'),'filepath',strcat(EEGloc, AnalyName));
    EEG = eeg_checkset( EEG );
    
end

% save edited behavior file
save([EEGloc,'\',AnalyName,'\','epoch'],'behavFile'); 
