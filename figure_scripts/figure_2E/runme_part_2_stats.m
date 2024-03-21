% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure 2E.
%
% Note that the figures were post-processed in Adobe Illustator.
close all
clear
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

LOCAL_OUTPUT_FOLDER = './local_output'; % for Figure S5B
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    error('The local folder should exist, but does not!');
end

% We will process each excel file
availableFiles = dir(fullfile('./local_output', '*_accuracies.xlsx'));
nFiles = length(availableFiles);
for iFile = 1:nFiles
    ifn = fullfile(availableFiles(iFile).folder, availableFiles(iFile).name);
    
    T = readtable(ifn);
    
    R = compute_descriptive_statistics(T);
    ofn = strrep(ifn, '_accuracies.xlsx', '_descriptive_stats.xlsx');
    writetable(R, ofn);
    fprintf('Results saved to %s\n', ofn);

    R = compute_ttests(T);
    ofn = strrep(ifn, '_accuracies.xlsx', '_ttests.xlsx');
    writetable(R, ofn);
    fprintf('Results saved to %s\n', ofn);
end % iFile

function [R] = compute_ttests(T)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);

    M = zeros(nGroups,1);
    E = zeros(nGroups,1);
    R = [];
    k = 1;
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        x = T.accuracy(ismember(T.groupLabel, groupLabel),:);

        [h,p, ci, stats] = ttest(x, 0.5);

        R(iGroup).groupLabel = groupLabel;
        R(iGroup).h = h;
        R(iGroup).p = p;
        R(iGroup).ci = ci;
        R(iGroup).tstat = stats.tstat;
        R(iGroup).df = stats.df;
        R(iGroup).sd = stats.sd;
    end
    R = struct2table(R);
end

function [DS] = compute_descriptive_statistics(T)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);

    M = zeros(nGroups,1);
    E = zeros(nGroups,1);
    DS = [];
    k = 1;
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        S = T(ismember(T.groupLabel, groupLabel),:);
        x = S.accuracy * 100;
        meanX = mean(x);
        errX = std(x) ./ sqrt(length(x));
        M(iGroup) = meanX;
        E(iGroup) = errX;

        DS(k).groupLabel = groupLabel;
        DS(k).mean = meanX;
        DS(k).stderr = errX;
        DS(k).nSamples = length(x);
        k = k + 1;
    end
    DS = struct2table(DS);
end

function [signRankTests, R, S] = compute_signrank_stats(R0, MIN_CELLS)
    R = R0;
    R(R.nCellsUsed < MIN_CELLS,:) = [];
    
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);

    signRankTests = [];
    tk = 1;
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
    
        x = R.accuracy(ismember(R.groupLabel, groupLabel));
    
        N = length(x);
        [p, h, ~] = signrank(x, 0.5, 'tail', 'right');
        signRankTests(tk).groupLabel = groupLabel;
        signRankTests(tk).N = N;
        signRankTests(tk).h = h;
        signRankTests(tk).p = p;
        tk = tk + 1;

        fprintf('%s, N=%d, h = %d, p = %0.3f\n', groupLabel, N, h, p);
    end % iGroup
    signRankTests = struct2table(signRankTests);

    S = sortrows(R, {'groupLabel', 'isCalcium'});
end % function
