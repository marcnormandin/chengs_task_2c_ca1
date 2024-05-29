% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes figure S12B
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

% Read in the table made by the other script
dataFn = fullfile(LOCAL_OUTPUT_FOLDER, 'predict_contexts_using_withheld_map.xlsx');
allV = readtable(dataFn);
fprintf('Loaded results from %s\n', dataFn);


%% Make the combined figure
hFig = figure('position', get(0, 'screensize'));
ax = [];

ax(1) = subplot(1,1,1);
plot_prediction_accuracy(allV)
set(gca, 'fontweight', 'bold', 'fontsize', 18);

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S12B', {'png', 'svg', 'fig'});

%% Added to export data to excel for natcomms
T = allV;
T = renamevars(T,["bestCorrelation_recomputed","groupLabel","recording_type"],["context_similarity","dayName","dataset"]);
T(:,ismember(T.Properties.VariableNames, {'animalSessionCellName', 'cellId', 'isStable', 'isValid', 'dayNum', 'sessionName', 'bestCorrelation', 'groupId'})) = [];
writetable(T, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_S12B.xlsx'), 'Sheet', 'figure_S12B');

function plot_prediction_accuracy(V)
    % Plot Accuracy 
    xp = 1:1:100;
    x = prctile(V.bestCorrelation_recomputed, xp);
    y = zeros(size(x));
    e = zeros(size(x));
    for k = 1:length(x)
        d = V.prediction_accuracy(V.bestCorrelation_recomputed <= x(k));
        y(k) = mean(d);
        e(k) = std(d,1) ./ sqrt(length(d));
    end

    y = y * 100;
    e = e * 100;

    % the plot
    xline(0.3, 'r', 'linewidth', 4, 'linestyle', ':')
    hold on
    plot(x,y,'k-', 'linewidth', 4);
    hold on
    plot(x, y+e, 'r-', 'linewidth', 2)
    plot(x, y-e, 'r-', 'linewidth', 2)
    xlabel('Context Similarity Threshold')
    ylabel('Prediction Accuracy (%)')
    grid on
    set(gca, 'fontsize', 18)
end % function
