% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This computes stats for Figure 6C.

close all
clear all
clc

LOCAL_OUTPUT_FOLDER = './local_output';

availableFiles = dir(fullfile(LOCAL_OUTPUT_FOLDER, 'train_predict_different_days_cells_stable_*.xlsx'));
nFiles = length(availableFiles);

% Process each file that has predictions
for iFile = 1:nFiles
    ifn = fullfile( availableFiles(iFile).folder, availableFiles(iFile).name );

    % Read in the prediction results
    T = readtable(ifn);
    fprintf('Read in the prediction results from %s\n', ifn);

    R = compute_descriptive_stats(T);
    ofn = strrep(ifn, '.xlsx', '_descriptive_stats.xlsx');
    writetable(R, ofn);
    fprintf('Stats saved to %s\n', ofn);

    R = compute_ttests(T);
    ofn = strrep(ifn, '.xlsx', '_ttests.xlsx');
    writetable(R, ofn);
    fprintf('Stats saved to %s\n', ofn);
end

function [R] = compute_descriptive_stats(T)
    % Get the day pairs (instead of assuming)
    DayPairs = unique(T(:, ismember(T.Properties.VariableNames, {'trainGroupLabel', 'predictGroupLabel'})), 'rows');
    nPairs = size(DayPairs,1);

    % ttests
    R = [];
    for iPair = 1:nPairs
        dayPair = table2struct(DayPairs(iPair,:));

        x = T.accuracy(ismember(T.trainGroupLabel, dayPair.trainGroupLabel) & ismember(T.predictGroupLabel, dayPair.predictGroupLabel));

        meanX = mean(x);
        stdX = std(x);
        nX = length(x);
        errX = stdX ./ sqrt(nX);

        R(iPair).trainGroupLabel = dayPair.trainGroupLabel;
        R(iPair).predictGroupLabel = dayPair.predictGroupLabel;
        R(iPair).meanAccuracy = meanX;
        R(iPair).stdAccuracy = stdX;
        R(iPair).nSamples = nX;
        R(iPair).stderrorAccuracy = errX;
    end
    R = struct2table(R)
end % function

function [R] = compute_ttests(T)
    % Get the day pairs (instead of assuming)
    DayPairs = unique(T(:, ismember(T.Properties.VariableNames, {'trainGroupLabel', 'predictGroupLabel'})), 'rows');
    nPairs = size(DayPairs,1);

    % ttests
    R = [];
    for iPair = 1:nPairs
        dayPair = table2struct(DayPairs(iPair,:));

        x = T.accuracy(ismember(T.trainGroupLabel, dayPair.trainGroupLabel) & ismember(T.predictGroupLabel, dayPair.predictGroupLabel))

        % We just do a ttest
        tail = 'right';
        [h, p, ci, stats] = ttest(x, 0.5, 'tail', tail);
        
        R(iPair).trainGroupLabel = dayPair.trainGroupLabel;
        R(iPair).predictGroupLabel = dayPair.predictGroupLabel;
        R(iPair).tail = tail;
        R(iPair).h = h;
        R(iPair).p = p;
        R(iPair).ci_lo = ci(1);
        R(iPair).ci_hi = ci(2);
        R(iPair).tstat = stats.tstat;
        R(iPair).df = stats.df;
        R(iPair).sd = stats.sd;
    end
    R = struct2table(R);
end % function
