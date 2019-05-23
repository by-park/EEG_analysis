%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG Preprocessing
%
%   -> re-referencing, filtering, resampling, ICA 
%
% created at 2019.04.24 PBY
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% parameters
% Name of the analysis
AnalyName = 'pp_ica';

% EEG file location
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\뇌파데이터\'; 
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

    % 2. load channel location file
    % EEG=pop_chanedit(EEG, 'lookup','C:\\Users\\BY\\Documents\MATLAB\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
    EEG=pop_chanedit(EEG, 'lookup','D:\\백업\\개인자료\\프로그래밍_대학원\\뇌파\\64chan.ced');

    % 3. reference - offline (M1 M2)
    EEG = pop_reref( EEG, [33 43] );
    
    % 4. band pass filter (1, 30 Hz)
    EEG = pop_eegfiltnew(EEG, 0.1,30,[],0,[],0); % final argument 0 for not plotting & filtorder is empty
    
    % 5. resample (down sample) : 1000Hz -> 500 Hz
    EEG = pop_resample( EEG, resamplerate);
    
    % 6. ICA
    EEG = pop_runica(EEG, 'extended',1,'interupt','off', 'icatype', 'binica');
    EEG = eeg_checkset( EEG );
    
    % baseline
    if ~exist(AnalyName, 'dir')
        mkdir(AnalyName);
    end
    EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(fileNames(file).name, "."), '_', AnalyName, '.set'),'filepath',strcat(EEGloc, AnalyName));
    EEG = eeg_checkset( EEG );
end
