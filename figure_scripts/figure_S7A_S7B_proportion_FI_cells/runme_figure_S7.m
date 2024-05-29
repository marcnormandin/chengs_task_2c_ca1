% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes the figure of Proportions of cell types by day and dataset
% (Figure S7).
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

% Calcium Results
calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calData = load(fullfile(calciumAnalysisResultsFolder, 'analysis_results.mat'));

tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');
tetData = load(fullfile(tetrodesAnalysisResultsFolder, 'analysis_results.mat'));

%%
[TTet, MTet, SETet, groupLabelsTet] =  process_dataset(tetData);
[TCal, MCal, SECal, groupLabelsCal] =  process_dataset(calData);

% TTet and TCal have the percentStable per animal and day
% Save the results to Excel files
writetable(TTet, fullfile(LOCAL_OUTPUT_FOLDER, 'tetrodes_percent_stable_per_animal_per_day.xlsx'));
writetable(TCal, fullfile(LOCAL_OUTPUT_FOLDER, 'calcium_percent_stable_per_animal_per_day.xlsx'));

% Save the descriptive statistics
T = array2table([MTet', SETet'], 'VariableNames', {'meanPercentStable', 'stderrorPercentStable'});
T.groupLabel = groupLabelsTet';
writetable(T, fullfile(LOCAL_OUTPUT_FOLDER, 'tetrodes_percent_stable_per_animal_per_day_descriptive_stats.xlsx'));

% Save the descriptive statistics
T = array2table([MCal', SECal'], 'VariableNames', {'meanPercentStable', 'stderrorPercentStable'});
T.groupLabel = groupLabelsCal';
writetable(T, fullfile(LOCAL_OUTPUT_FOLDER, 'calcium_percent_stable_per_animal_per_day_descriptive_stats.xlsx'));



close all

hFig = figure('position', get(0,'screensize'));

ax(1) = subplot(1,2,1);
bar(1:length(groupLabelsTet), MTet)
hold on
errorbar(1:length(groupLabelsTet), MTet, SETet, 'color', 'k', 'linestyle', 'none', 'linewidth', 4);
xticklabels(groupLabelsTet);
set(gca, 'fontsize', 18, 'fontweight', 'bold');
a = axis;
axis([a(1), a(2), 0, 100]);
grid on
title('Electrophysiology')
ylabel('Percent FI')

ax(2) = subplot(1,2,2);
bar(1:length(groupLabelsCal), MCal)
hold on
errorbar(1:length(groupLabelsCal), MCal, SECal, 'color', 'k', 'linestyle', 'none', 'linewidth', 4);
xticklabels(groupLabelsCal);
set(gca, 'fontsize', 18, 'fontweight', 'bold');
a = axis;
axis([a(1), a(2), 0, 100]);
grid on
title('Calcium-Imaging')
ylabel('Percent FI')

fnPrefix = 'figure_S7_percents_FI_vs_day';
ml_savefig(hFig, OUTPUT_FOLDER, fnPrefix, {'png', 'svg', 'fig'});

%% Added to export to excel for natcomms
% Just read in the data and format it for natcomms
XT = readtable(fullfile(pwd, 'local_output', 'tetrodes_percent_stable_per_animal_per_day.xlsx'));
XT.datasetName = repmat(string('Tetrodes'), height(XT), 1);
XC = readtable(fullfile(pwd, 'local_output', 'calcium_percent_stable_per_animal_per_day.xlsx'));
XC.datasetName = repmat(string('Calcium'), height(XC), 1);

XA = [XT; XC];
XA = renamevars(XA,["groupLabel"],["dayName"]);
writetable(XA, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_S7.xlsx'), 'Sheet', 'figure_S7');

%% Functions
function [T, M,SE,groupLabels] =  process_dataset(dataset)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    BA = dataset.analysisResults.BestAligned;
    BA = helper_add_group_labels_to_table(BA, load_sessions_to_groups_table(dataset.analysisSettings), true);
    BA = sortrows(BA, 'groupId');
    BA(~ismember(BA.groupLabel, groupLabels),:) = [];
    
    
    animalNames = unique(BA.animalName);
    nAnimals = length(animalNames);
    nGroupLabels = length(groupLabels);
    
    % Make the table of results per animal
    T = [];
    k = 1;
    for iAnimal = 1:nAnimals
        animalName = animalNames{iAnimal};
        for iGroup = 1:nGroupLabels
            groupLabel = groupLabels{iGroup};
    
            BAMatch = BA(ismember(BA.animalName, animalName) & ismember(BA.groupLabel, groupLabel),:);
            BAMatch(isnan(BAMatch.bestCorrelation),:) = [];

            nCells = size(BAMatch,1);
            nStable = sum(BAMatch.bestCorrelation > 0.3);
    
            T(k).animalName = animalName;
            T(k).groupLabel = groupLabel;
            T(k).percentStable = nStable / nCells * 100;
            k = k + 1;
        end % iGroup
    end % iAnimal
    T = struct2table(T);
    
    % Arrange the data for the bar plot (take the averages per day)
    M = nan(1, nGroupLabels);
    SE = nan(1, nGroupLabels);
    for iGroup = 1:nGroupLabels
        groupLabel = groupLabels{iGroup};
    
        percents = T.percentStable(ismember(T.groupLabel, groupLabel));
        percents(isnan(percents)) = [];
        
        n = length(percents);
        M(iGroup) = mean(percents);
        SE(iGroup) = std(percents) ./ sqrt(n);
    end % iGroup
end % function
