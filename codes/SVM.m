%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERP analysis
%
%   -> averging EEG singal for SVM and permutation test
%
% created at 2019.05.19 PBY (현재 iteration 과정이 없음!)
% github: https://github.com/BY1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
% clear all
% close all
% clc

%% parameters
% Name of the analysis
AnalyName = 'SVM';

% EEG file location
EEGloc = 'D:\EEG 파일 모음\내 실험\Ex1_n2pc_tar\뇌파데이터\pp_ica_removed\ep\'; 
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex2_n2pc_dist\뇌파데이터\pp_ica_removed\ep\';
% EEGloc = 'D:\EEG 파일 모음\내 실험\Ex7_n2pc_tar_single\뇌파데이터\pp_ica_removed\ep\';

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
cond1 = 4; % cue type (4: ex1, ex2 / 3: ex7)
cond2 = 2; % ipsil/cont

% for averaging
total_subplot = zeros(cond1*cond2, round(grange2*resamplerate) + round(-1*grange1*resamplerate));
AnovaResult = zeros(foldersize*cond1*cond2, 4); % sub, cue type, ipsil/cont, eeg
Average_grand = zeros(2,round(grange2*resamplerate) + round(-1*grange1*resamplerate)); % 2 = PO7, PO8
nGrandAverage = 0;
nAnova = 1;

% for machine learning
fornnstart = [];

% graph names for visualization
graphNames = {'Distractor(blue)','Target(green)', 'Target(yellow)','Neutral(red)'};
% graphNames = {'Target(blue)','Distractor(green)', 'Distractor(yellow)','Neutral(red)'};
% graphNames = {'Distractor(red)','Target(blue)', 'Target(yellow)','Neutral(green)'};
% graphNames = {'Distractor(red)','Target(blue)', 'Target(green)','Neutral(yellow)'};

%% permutation matrix
permutationmatrix = [];


%% load behav file
behavFile = load([EEGloc, 'epoch.mat']);
behavFile = behavFile.behavFile;
trialBehav = 1;

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
    
    % 2. averaging
    temp_raw = EEG.data(:,round(baserange*resamplerate)+round(0*0.001*resamplerate):1:round(baserange*resamplerate)+round(300*0.001*resamplerate), 1);
    fornnstart_comb = zeros(size(1:2:size(EEG.event,2),2), 64+2, size(3:4:size(temp_raw,2)-2,2));
    for removedTrialEEG = 1: 2 : size(EEG.event,2) 
        % trigger at the trial
        trigger = EEG.event(removedTrialEEG).type;
        
        % SVM time points 0 ms ~ 300 ms
        temp_raw = EEG.data(:,round(baserange*resamplerate)+round(0*0.001*resamplerate):1:round(baserange*resamplerate)+round(300*0.001*resamplerate), (removedTrialEEG+1)/2);
        temp = zeros(64,size(3:4:size(temp_raw,2)-2,2));
        
        timepoint = 1;
        for temppoint = 3 : 4 : size(temp_raw,2)-2
            temp(:,timepoint) = mean(temp_raw(:,temppoint-2:temppoint+2),2);
            tempt = temp(:,timepoint)' ;    
            fornnstart_comb((removedTrialEEG+1)/2, :, timepoint) = [tempt, rem(trigger,10), behavFile(trialBehav,6)];
            
            timepoint = timepoint + 1;
        end
        trialBehav = trialBehav + 1;            
    end
      
    % 3. SVM training
    ttestmatrix = zeros(10,4,size(fornnstart_comb,3)); % 30 samples with 3 cue color at 5 time points

    for timepoint = 1:size(fornnstart_comb,3)
        fprintf('** Timepoint %d \n',timepoint)
        randind_target = randperm(size(fornnstart_comb,1));
        fornnstart2_target = [fornnstart_comb(:,:,timepoint) randind_target'];
        fornnstart2_target2 = [fornnstart_comb(:,:,timepoint) randind_target'];
        
        randind_dist = randperm(size(fornnstart_comb,1));
        fornnstart2_dist = [fornnstart_comb(:,:,timepoint) randind_dist'];        

        randind_neutral = randperm(size(fornnstart_comb,1));
        fornnstart2_neutral = [fornnstart_comb(:,:,timepoint) randind_neutral'];        
        
        % target list: remove other cues except target cues
%         fornnstart2_target((fornnstart2_target(:,end-2)==2) | (fornnstart2_target(:,end-2)==3) | (fornnstart2_target(:,end-2)==4),:) = [];

        fornnstart2_target((fornnstart2_target(:,end-2)==1) | (fornnstart2_target(:,end-2)==4)| (fornnstart2_target(:,end-2)==3),:) = [];
        fornnstart2_target2((fornnstart2_target2(:,end-2)==1) | (fornnstart2_target2(:,end-2)==4)| (fornnstart2_target2(:,end-2)==2),:) = [];

        % distractor list (23)
%         fornnstart2_dist((fornnstart2_dist(:,end-2)==1) | (fornnstart2_dist(:,end-2)==4),:) = [];
        fornnstart2_dist((fornnstart2_dist(:,end-2)==2) | (fornnstart2_dist(:,end-2)==3)  | (fornnstart2_dist(:,end-2)==4),:) = [];

        % neutral list (4)
        fornnstart2_neutral((fornnstart2_neutral(:,end-2)==1) | (fornnstart2_neutral(:,end-2)==2) | (fornnstart2_neutral(:,end-2)==3),:) = [];
%         fornnstart2_neutral((fornnstart2_neutral(:,end-2)==1) | (fornnstart2_neutral(:,end-2)==2),:) = [];

        % target sorting
        [~,sortOrder] = sort(fornnstart2_target(:,end),'descend');        
        out = fornnstart2_target(sortOrder,:);
        newout = [];
        for i = 1:4 % four loca
            outind = (out(:,end-1)==i);
            outplus = out(outind,1:end-3);
            newout = [newout; outplus, repmat(i,sum(outind),1)];
%             for j = 1:10:size(out,1)-9 % combine 10 trials
%                 outind = (out(:,end-1)==i);
%                 outplus = mean(out(outind(j:j+9),1:end-3),1);
%                 newout = [newout; outplus, i];
%             end
        end
        randind2 = randperm(size(newout,1));
        newout = [newout randind2'];
        [~,sortOrder] = sort(newout(:,end),'descend');
        newout2_target = newout(sortOrder,:);
        nanindex = isnan(newout2_target(:,1)); % delete NAN
        newout2_target(nanindex,:) = [];
        
        % target sorting2
        [~,sortOrder] = sort(fornnstart2_target2(:,end),'descend');        
        out = fornnstart2_target2(sortOrder,:);
        newout = [];
        for i = 1:4 % four loca
            outind = (out(:,end-1)==i);
            outplus = out(outind,1:end-3);
            newout = [newout; outplus, repmat(i,sum(outind),1)];
%             for j = 1:10:size(out,1)-9 % combine 10 trials
%                 outind = (out(:,end-1)==i);
%                 outplus = mean(out(outind(j:j+9),1:end-3),1);
%                 newout = [newout; outplus, i];
%             end
        end
        randind2 = randperm(size(newout,1));
        newout = [newout randind2'];
        [~,sortOrder] = sort(newout(:,end),'descend');
        newout2_target2 = newout(sortOrder,:);
        nanindex = isnan(newout2_target2(:,1)); % delete NAN
        newout2_target2(nanindex,:) = [];        
        
        % distractor sorting
        [~,sortOrder] = sort(fornnstart2_dist(:,end),'descend');        
        out = fornnstart2_dist(sortOrder,:);
        newout = [];
        for i = 1:4 % four loca
            outind = (out(:,end-1)==i);
            outplus = out(outind,1:end-3);
            newout = [newout; outplus, repmat(i,sum(outind),1)];
%             for j = 1:10:size(out,1)-9 % combine 10 trials
%                 outind = (out(:,end-1)==i);
%                 outplus = mean(out(outind(j:j+9),1:end-3),1);
%                 newout = [newout; outplus, i];
%             end
        end
        randind2 = randperm(size(newout,1));
        newout = [newout randind2'];
        [~,sortOrder] = sort(newout(:,end),'descend');
        newout2_dist = newout(sortOrder,:);
        nanindex = isnan(newout2_dist(:,1)); % delete NAN
        newout2_dist(nanindex,:) = [];
        
        % neutral sorting
        [~,sortOrder] = sort(fornnstart2_neutral(:,end),'descend');        
        out = fornnstart2_neutral(sortOrder,:);
        newout = [];
        for i = 1:4 % four loca
            outind = (out(:,end-1)==i);
            outplus = out(outind,1:end-3);
            newout = [newout; outplus, repmat(i,sum(outind),1)];
%             for j = 1:10:size(out,1)-9 % combine 10 trials
%                 outind = (out(:,end-1)==i);
%                 outplus = mean(out(outind(j:j+9),1:end-3),1);
%                 newout = [newout; outplus, i];
%             end
        end
        randind2 = randperm(size(newout,1));
        newout = [newout randind2'];
        [~,sortOrder] = sort(newout(:,end),'descend');
        newout2_neutral = newout(sortOrder,:);
        nanindex = isnan(newout2_neutral(:,1)); % delete NAN
        newout2_neutral(nanindex,:) = [];        
        
        % 1/10 list for each cue
        target_ten = 1:floor(size(newout2_target,1)*1/10):size(newout2_target,1)-floor(size(newout2_target,1)*1/10)+1;
        target_ten2 = 1:floor(size(newout2_target2,1)*1/10):size(newout2_target2,1)-floor(size(newout2_target2,1)*1/10)+1;
        dist_ten = 1:floor(size(newout2_dist,1)*1/10):size(newout2_dist,1)-floor(size(newout2_dist,1)*1/10)+1;
        neutral_ten = 1:floor(size(newout2_neutral,1)*1/10):size(newout2_neutral,1)-floor(size(newout2_neutral,1)*1/10)+1;
        
        % cross validation 10 fold
        crossorder=1;
        for cross = 1:10
        fprintf('** crossvalid %d \n',crossorder)
        % target training
        cuecol = 2;
        % target1
        crossind = target_ten(cross):target_ten(cross)+floor(size(newout2_target,1)*1/10)-1;        
        TestSet = newout2_target(crossind,1:end-2);
        TestAnswer = newout2_target(crossind,end-1);
        
        % target2
        crossind2 = target_ten2(cross):target_ten2(cross)+floor(size(newout2_target2,1)*1/10)-1;        
        TestSet2 = newout2_target2(crossind2,1:end-2);
        TestAnswer2 = newout2_target2(crossind2,end-1);
        
        % target 1 training set
        TrainingSet = newout2_target(:,1:end-2);
        GroupTrain = newout2_target(:,end-1);        
        TrainingSet(crossind,:)=[];
        GroupTrain(crossind,:)=[];
        % target 2 training set
        TrainingSet2 = newout2_target2(:,1:end-2);
        GroupTrain2 = newout2_target2(:,end-1);        
        TrainingSet2(crossind2,:)=[];
        GroupTrain2(crossind2,:)=[];
        % combine 2 training sets
        TrainingSet = [TrainingSet; TrainingSet2];
        GroupTrain = [GroupTrain; GroupTrain2];
        
        u=unique(GroupTrain);
        numClasses=length(u);
        result = zeros(length(TestSet(:,1)),1);        
        Mdl = fitcecoc(TrainingSet, GroupTrain);
        
        % target test
        for j=1:size(TestSet,1)
            result(j,1) = predict(Mdl,TestSet(j,:));
        end
        % save result
        ttestmatrix(crossorder, cuecol, timepoint) = sum(result == TestAnswer)/size(result,1);
        
        % target2 test
        cuecol = 4;
        crossind = target_ten2(cross):target_ten2(cross)+floor(size(newout2_target2,1)*1/10)-1;        
        TestSet = newout2_target2(crossind,1:end-2);
        TestAnswer = newout2_target2(crossind,end-1);
        result = zeros(length(TestSet(:,1)),1);   
        for j=1:size(TestSet,1)
            result(j,1) = predict(Mdl,TestSet(j,:));
        end
        ttestmatrix(crossorder, cuecol, timepoint) = sum(result == TestAnswer)/size(result,1);
        
        % distractor test
        cuecol = 1;
        crossind = dist_ten(cross):dist_ten(cross)+floor(size(newout2_dist,1)*1/10)-1;        
        TestSet = newout2_dist(crossind,1:end-2);
        TestAnswer = newout2_dist(crossind,end-1);
        result = zeros(length(TestSet(:,1)),1);   
        for j=1:size(TestSet,1)
            result(j,1) = predict(Mdl,TestSet(j,:));
        end
        ttestmatrix(crossorder, cuecol, timepoint) = sum(result == TestAnswer)/size(result,1);
        
        % neutral test
        cuecol = 3;
        crossind = neutral_ten(cross):neutral_ten(cross)+floor(size(newout2_neutral,1)*1/10)-1;        
        TestSet = newout2_neutral(crossind,1:end-2);
        TestAnswer = newout2_neutral(crossind,end-1);
        result = zeros(length(TestSet(:,1)),1);   
        for j=1:size(TestSet,1)
            result(j,1) = predict(Mdl,TestSet(j,:));
        end
        ttestmatrix(crossorder, cuecol, timepoint) = sum(result == TestAnswer)/size(result,1);                
        
        crossorder=crossorder+1;
        end
    end
    
    permutationmatrix = [permutationmatrix;mean(ttestmatrix, 1)];
        
end
    

%% permutation test
cd('C:\Users\BY\Documents\Github\EEG_Decoding\codes');

perm_pvalue = zeros(4,size(fornnstart_comb,3));
for timepoint = 1:size(fornnstart_comb,3)
    for cuecol = 1:4
        fprintf('one sample t-test timepoint %d cue color %d \n', timepoint, cuecol)
        % ttest
        % [h,p,ci,stats] = ttest(ttestmatrix(:,cuecol,timepoint),0.25,'Alpha',0.05,'Tail','right');
        % permutation test
        chance = repmat([0.25], [size(permutationmatrix(:, cuecol, timepoint), 1), 1]);
        [p, ~, ~] = permutationTest(permutationmatrix(:, cuecol, timepoint), chance, 1000)
        perm_pvalue(cuecol,timepoint) =p; 
    end
end
cd(EEGloc);
if ~exist(AnalyName, 'dir')
    mkdir(AnalyName);
end
save([AnalyName,'\permutationmatrix'],'permutationmatrix'); 

%% graphs of accuracy
% cue col 4
cd('C:\Users\BY\Documents\Github\EEG_Decoding\codes');

f1 = figure;
cuecol = 2; % target 1
aexample = zeros(foldersize,15);
for timepoint=1:15
    aexample(:,timepoint) = permutationmatrix(:,cuecol,timepoint);
end
%stdshade(aexample, 0.5,'g');
plot(mean(aexample, 1), 'g');
hold on;
plot([1 15],[0.25 0.25],'k');
axis([1 15 0 0.35])

f2 = figure;
cuecol = 1; % distractor
aexample = zeros(foldersize,15);
for timepoint=1:15
    aexample(:,timepoint) = permutationmatrix(:,cuecol,timepoint);
end
%stdshade(aexample, 0.5,'b');
plot(mean(aexample, 1), 'b');
hold on;
plot([1 15],[0.25 0.25],'k');
axis([1 15 0 0.35])

f3 = figure;
cuecol = 3; % neutral
aexample = zeros(foldersize,15);
for timepoint=1:15
    aexample(:,timepoint) = permutationmatrix(:,cuecol,timepoint);
end
%stdshade(aexample, 0.5,'r');
plot(mean(aexample, 1), 'r');
hold on;
plot([1 15],[0.25 0.25],'k');
axis([1 15 0 0.35])

f4 = figure;
cuecol = 4; % target2
aexample = zeros(foldersize,15);
for timepoint=1:15
    aexample(:,timepoint) = permutationmatrix(:,cuecol,timepoint);
end
%stdshade(aexample, 0.5,'y');
plot(mean(aexample, 1), 'y');
hold on;
plot([1 15],[0.25 0.25],'k');
axis([1 15 0 0.35])

% save graphs
cd(EEGloc);
saveas(f1,[AnalyName,'\SVM_target1'],'jpg');   
saveas(f2,[AnalyName,'\SVM_distractor'],'jpg');   
saveas(f3,[AnalyName,'\SVM_neutral'],'jpg');   
saveas(f4,[AnalyName,'\SVM_target2'],'jpg');   


