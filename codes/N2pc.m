%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERP analysis
%
%   -> averging N2pc components and visualization
%
% created at 2019.05.12 PBY
% updated at 2019.06.26 PBY
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% parameters
% Name of the analysis
AnalyName = 'N2pc';

% EEG file location
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\뇌파데이터\pp_ica_removed\ep\'; 
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex2_n2pc_dist\뇌파데이터\pp_ica_removed\ep\';
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\뇌파데이터\pp_ica_removed\ep\';

cd(EEGloc);

% find Files
runStr = '*.set';
fileNames = dir(runStr);
foldersize = size(fileNames,1);

% electrodes for N2pc
LeftElec = 51; %PO7
RightElec = 57 ; %PO8

% sampling rate
resamplerate = 500;

% baseline
baserange = 100;
baserange = baserange * 0.001;

% range for N2pc
lrange = 180; % lower range
lrange = lrange * 0.001;
urange = 240; % upper range
urange = urange * 0.001;

% range for epoch (for graph)
grange1 = -100; % -0.1
grange1 = grange1 * 0.001;
grange2 = 500; % 0.5
grange2 = grange2 * 0.001;

% manipulated conditions
cond1 = 3; % cue type (4: ex1, ex2 / 3: ex7)
cond2 = 2; % ipsil/cont

% for averaging
total_subplot = zeros(cond1*cond2, round(grange2*resamplerate) + round(-1*grange1*resamplerate));
AnovaResult = zeros(foldersize*cond1*cond2, 4); % sub, cue type, ipsil/cont, eeg
Average_grand = zeros(2,round(grange2*resamplerate) + round(-1*grange1*resamplerate)); % 2 = PO7, PO8
nGrandAverage = 0;
nAnova = 1;

% for topography
fortopography = [];

% graph names for visualization
% graphNames = {'Distractor(blue)','Target(green)', 'Target(yellow)','Neutral(red)'};
% graphNames = {'Target(blue)','Distractor(green)', 'Distractor(yellow)','Neutral(red)'};
% graphNames = {'Distractor(red)','Target(blue)', 'Target(yellow)','Neutral(green)'};
% graphNames = {'Distractor(red)','Target(blue)', 'Target(green)','Neutral(yellow)'};
graphNames = {'Target(red)','Distractor(green)', 'Neutral(blue)'};


%% load behav file
behavFile = load([EEGloc, 'epoch.mat']);
behavFile = behavFile.behavFile;
trialBehav = 0;

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
    
    % 2. avearging
    plotResult_LeftElec = zeros(cond1, round(grange2*resamplerate)-round(grange1*resamplerate), cond2);
    plotResult_RightElec = zeros(cond1, round(grange2*resamplerate)-round(grange1*resamplerate), cond2);
    plotResult_AllElec = zeros(cond1, round(grange2*resamplerate)-round(grange1*resamplerate));
    checkNum = zeros(cond1,cond2); % 1234 bgyr (dist, tar, tar, neu) & ipsil/cont
    
    for removedTrialEEG = 1: 2 : size(EEG.event,2)
        
        % trigger at the trial
        trigger = EEG.event(removedTrialEEG).type;

        % for topography
        trialBehav = trialBehav + 1;
        temp = EEG.data(:,round(baserange*resamplerate)+round(100*0.001*resamplerate)+1:round(baserange*resamplerate)+round(400*0.001*resamplerate), (removedTrialEEG+1)/2);
        fortopography = [fortopography; temp(:)', behavFile(trialBehav,6), rem(trigger,10)]; % time point, location, cue type
        
        % exp7 trigger
        if (floor(trigger/10) == 1) || (floor(trigger/10) == 3)
            continue;
        elseif floor(trigger/10) == 4
            trigger = trigger-30;
        end
        
        if trigger < 90 % only for left and right electrdes        
            % plotResult_LeftElec(trigger 뒷자리수 = cue type, eeg, 앞자리수 = location)
            plotResult_LeftElec(rem(trigger,10),...
            :,floor(trigger/10))...
            ...
            = plotResult_LeftElec(rem(trigger,10),...
            :,floor(trigger/10)) + ...
            mean(EEG.data(LeftElec,:,(removedTrialEEG+1)/2),1);

            plotResult_RightElec(rem(trigger,10),...
            :,floor(trigger/10))...
            ...
            = plotResult_RightElec(rem(trigger,10),...
            :,floor(trigger/10)) + ...
            mean(EEG.data(RightElec,:,(removedTrialEEG+1)/2),1);

            plotResult_AllElec(rem(trigger,10),:)...
            = plotResult_AllElec(rem(trigger,10),:)...
            + mean(EEG.data([LeftElec RightElec], :, (removedTrialEEG+1)/2),1);

            % counting for averaging
            checkNum(rem(trigger,10),floor(trigger/10)) ...
             = checkNum(rem(trigger,10),floor(trigger/10)) + 1;

        
            % Grand average (앞 -100ms 뒤 500ms해서 전체 154 데이터 포인트)
            Average_grand = Average_grand + mean(EEG.data([LeftElec RightElec],:,:),3);
            % counting grand averaging
            nGrandAverage = nGrandAverage + 1;
        
        end       
    
    end
    
    % 3. plot for each sub
    h = figure;
    
    subplot_matrix = zeros(cond1*cond2, round(grange2*resamplerate) + round(-1*grange1*resamplerate));
    AnovaResult_sub = zeros(cond1*cond2, 3); % cond1, cond2, eeg = 3    
        
    for nPlot = 1:cond1
        % matrix for plot
        subplot_matrix(nPlot*2-1,:) = (plotResult_LeftElec(nPlot,:,1)./checkNum(nPlot,1) + plotResult_RightElec(nPlot,:,2)./checkNum(nPlot,2))./2;
        subplot_matrix(nPlot*2,:) = (plotResult_LeftElec(nPlot,:,2)./checkNum(nPlot,2) + plotResult_RightElec(nPlot,:,1)./checkNum(nPlot,1))./2;
        
        % plot 2 lines
        subplot(cond1, 1, nPlot);
        % ipsilateral
        plot(grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000,subplot_matrix(nPlot*2-1, :),'--');
        hold on;
        % contralateral
        plot(grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000,subplot_matrix(nPlot*2,:));
                
        % parameters for plot
        set(gca,'YDir','Reverse')
        legend('ipsilateral','contralateral')
        title(graphNames(nPlot))

        % 전체 피험자용 그래프 grand average
        total_subplot(nPlot*2-1,:) = total_subplot(nPlot*2-1,:) +  subplot_matrix(nPlot*2-1,:);
        total_subplot(nPlot*2,:) = total_subplot(nPlot*2,:) +  subplot_matrix(nPlot*2,:);

        % for anova
        AnovaResult_sub(nPlot*2-1:nPlot*2, 1) = nPlot; 
        AnovaResult_sub(nPlot*2-1,2) = 1; 
        AnovaResult_sub(nPlot*2,2) = 2;
        AnovaResult_sub(nPlot*2-1, 3) = mean(subplot_matrix(nPlot*2-1,round(baserange*resamplerate)+round(lrange*resamplerate) ...
        :round(baserange*resamplerate)+round(urange*resamplerate)));
        AnovaResult_sub(nPlot*2, 3) = mean(subplot_matrix(nPlot*2,round(baserange*resamplerate)+round(lrange*resamplerate) ...
        :round(baserange*resamplerate)+round(urange*resamplerate)));
    end   
        
    % Fill Anova matrix
    AnovaResult(nAnova:nAnova+cond1*cond2-1, 1) = file; % sub
    AnovaResult(nAnova:nAnova+cond1*cond2-1, 2:4) = AnovaResult_sub;
    
    nAnova = nAnova + cond1*cond2; % cond1*cond2 = 8
    
    % 6. save
    if ~exist(AnalyName, 'dir')
        mkdir(AnalyName);
    end
    
    % eeg files are not changed
    % EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(fileNames(file).name, "."), '_', AnalyName, '.set'),'filepath',strcat(EEGloc, AnalyName));
    % EEG = eeg_checkset( EEG );
    
    % save graphs
    saveas(h,[AnalyName,'\','sub', num2str(file), AnalyName],'jpg')
    
end

% Grand average plot (total subject!)
hold off;
g = figure;
plot(grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000,mean(Average_grand)/nGrandAverage);
title('Grand Average');
set(gca,'YDir','Reverse') % reverse axis!!

% save
saveas(g,[AnalyName,'\','Grand Average',AnalyName],'jpg')

% for excel graph (19.06.26 수정)
N2pc_matrix = [grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000];

% Grand average for each cue type
g2 = figure;
    for nGrandPlot = 1:cond1
        subplot(cond1,1,nGrandPlot);
        % ipsil
        plot(grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000,total_subplot(nGrandPlot*2-1,:)./foldersize);
        hold on;
        % cont
        plot(grange1*1000:(grange2*1000-grange1*1000)/(round(grange2*resamplerate)-round(grange1*resamplerate)-1):grange2*1000,total_subplot(nGrandPlot*2,:)./foldersize);
        set(gca,'YDir','Reverse')
        legend('ipsilateral','contralateral','Location','NorthWest')
        ylim([-10 15])
        title(graphNames(nGrandPlot))  
        
        N2pc_matrix = [N2pc_matrix;total_subplot(nGrandPlot*2-1,:)./foldersize; total_subplot(nGrandPlot*2,:)./foldersize];

    end
    
% save figure
saveas(g2,[AnalyName,'\','Grand Average_cue color',AnalyName],'jpg');   

% save excel
newfilename = [AnalyName,'\','AnovaResult',AnalyName,'.xlsx'];
xlswrite(newfilename,AnovaResult);

plots = [AnalyName,'\','N2pc_matrix',AnalyName,'.xlsx'];
xlswrite(plots,N2pc_matrix);

save([AnalyName,'\','fortopography_',AnalyName],'fortopography'); 