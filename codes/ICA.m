%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG ICA
%
%   -> run ICA 
%
% created at 2019.04.25 PBY
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% parameters
% Name of the analysis
AnalyName = 'ica';

% EEG file location
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\뇌파데이터\pp\'; 
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex2_n2pc_dist\뇌파데이터\pp\';
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\뇌파데이터\pp\';

cd(EEGloc);

% find Files
runStr = '*.set';
fileNames = dir(runStr);
foldersize = size(fileNames,1);

% sampling rate
resamplerate = 500;

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
    
    % 2. ICA
    EEG = pop_runica(EEG, 'extended',1,'interupt','off');
    EEG = eeg_checkset( EEG );
    
    % 3. save
    if ~exist(AnalyName, 'dir')
        mkdir(AnalyName);
    end
    EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(fileNames(file).name, "."), '_', AnalyName, '.set'),'filepath',strcat(EEGloc, AnalyName));
    EEG = eeg_checkset( EEG );
end
