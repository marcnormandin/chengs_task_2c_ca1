% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure S14 (first run the other script). It was previously
% made in Excel and this version can be made in Matlab.
close all
clear all
clc

run('../load_figure_config.m')

% Load the Excel file that the previous code made
INPUT_FOLDER = '.\local_output'; % where the Excel file is located.
fn = fullfile(INPUT_FOLDER, 'calcium_rate_matrices_within_across_per_animal_session.xlsx');

OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

T = readtable(fn);

groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroupLabels = length(groupLabels);
comparisonTypes = {'within', 'across'};
nComparisonTypes = length(comparisonTypes);

M = nan(nGroupLabels, 2);
SE = nan(nGroupLabels, 2);
for iGroup = 1:nGroupLabels
    groupLabel = groupLabels{iGroup};

    for iComparison = 1:nComparisonTypes
        comparisonType = comparisonTypes{iComparison};
        varName = sprintf('%sMean', comparisonType);
        x = T.(varName)(ismember(T.groupLabel, groupLabel));
        nX = length(x);

        M(iGroup, iComparison) = mean(x);
        SE(iGroup, iComparison) = std(x) / sqrt(nX);
    end % iComparison
end % iGroup

Y = M;
E = SE;
colours2 = colour_matrix_yellow_green();

hFig = figure('position', get(0, 'screensize'));
b = barplot_with_errors(Y', E', colours2);
ylabel('Absolute Rate Difference (Hz)')
xticks(1:nGroupLabels);
xticklabels(groupLabels)
legend({'Within Context', 'Across Context'}, 'location', 'southoutside', 'orientation', 'horizontal');
set(gca, 'fontsize', 18, 'fontweight', 'bold')

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S14', {'png', 'svg'});

%% Added to export data to excel for natcomms
NC = T;
NC(:, ismember(NC.Properties.VariableNames, {'sessionName', 'groupId', 'withinError', 'acrossError'})) = [];
NC = renamevars(NC, {'groupLabel', 'acrossMean', 'withinMean'}, {'dayName', 'absolute_rate_difference_within', 'absolute_rate_difference_across'});
NC
writetable(NC, fullfile(OUTPUT_FOLDER, "natcomms_excel_figure_S14.xlsx"), 'Sheet', 'figure_S14');
