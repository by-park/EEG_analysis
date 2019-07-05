%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topography
%
%   -> EEG visualization
%
% created at 2019.07.06 PBY
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% Topography
% load data
eeglab;
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\뇌파데이터\pp_ica_removed\ep\';
eloc = readlocs('D:\EEG 파일 모음\topography\64chan.ced');
load([EEGloc, 'N2pc\fortopography_N2pc.mat']);

% plot
for numfig = 1 % target
    % 150 timepoints (100 ~ 400 ms) => linspace(100, 400, 150)
    f1 = figure;
    for numcond = 1:4 % 4 locations
        subplot(1,4,numcond)
        figdata = zeros(1,64);
        for timepoint = 41:70 % 180 ms ~ 240 ms
            figdata = figdata + mean(fortopography((fortopography(:,9602,1)==numfig)&(fortopography(:,9601,1)==numcond),(timepoint-1)*64+1-2:timepoint*64-2),1);
        end
        plotData = figdata./30;
        topoplot(plotData, eloc,'electrodes','off','headrad', 0.44,'plotrad',0.45,'maplimits',[-6 6]); % option: 'maplimits', 'maxmin'
        axis([-0.6 0.6 -0.6 0.6]);
        cbar('vert',0,[-6 6]); % [-1 1]*max(abs(plotData))
    end
end

% save figure
saveas(f1,[EEGloc, 'N2pc\topography'],'jpg');   
