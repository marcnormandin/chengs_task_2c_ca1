% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This code makes Figures S2B and S2C, which show histograms of occupancy.
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

tmp = load(fullfile(INPUT_FOLDER, 'ratemaps-2021_12_1_12_34_49_K1_d4.mat'));

USE_SMOOTHED_MAPS = false;
TTetrodes = tetrodes_compute_fraction_occupied(tmp.ratemaps, USE_SMOOTHED_MAPS);
TCalcium = calcium_compute_fraction_occupied(INPUT_FOLDER, USE_SMOOTHED_MAPS);
clear tmp

eliminateUnmatched = true;
TTetrodes = helper_add_group_labels_to_table(TTetrodes, load_sessions_to_groups_table(false), eliminateUnmatched);
TCalcium = helper_add_group_labels_to_table(TCalcium, load_sessions_to_groups_table(true), eliminateUnmatched);

%save(fullfile(OUTPUT_FOLDER, 'figure_S12A_trial_map_fraction_occupied.mat'), 'TTetrodes', 'TCalcium');

%% ARENA BINS VISITED
close all
hFig = figure('position', get(0, 'screensize'));

FRACTION = 0.6; % This means that it must be this fraction or greater for the sum

tetrodesTotal = length(TTetrodes.fractionOccupied);
tetrodesSumAbove = sum(TTetrodes.fractionOccupied > FRACTION);
tetrodesPercentAbove = tetrodesSumAbove / tetrodesTotal * 100;

calciumTotal = length(TCalcium.fractionOccupied);
calciumSumAbove = sum(TCalcium.fractionOccupied > FRACTION);
calciumPercentAbove = calciumSumAbove / calciumTotal * 100;

binEdges = 0:0.1:1; % bin width is 0.1

FONT_SIZE = 12;
ax = [];
ax(1) = subplot(1,2,1);
histogram(TTetrodes.fractionOccupied, binEdges);
xlabel('Fraction of Arena Bins Visited')
ylabel('Counts (1 per trial)')
grid on
title(sprintf('Tetrodes: Total = %d, Above %0.1f =  %d, Percent Above = %0.1f', tetrodesTotal, FRACTION, tetrodesSumAbove, tetrodesPercentAbove))
set(gca, 'fontsize', FONT_SIZE, 'fontweight', 'bold');
axis tight

ax(2) = subplot(1,2,2);
histogram(TCalcium.fractionOccupied, binEdges);
xlabel('Fraction of Arena Bins Visited')
ylabel('Counts (1 per trial)')
grid on
title(sprintf('Calcium: Total = %d, Above %0.1f =  %d, Percent Above = %0.1f', calciumTotal, FRACTION, calciumSumAbove, calciumPercentAbove))
set(gca, 'fontsize', FONT_SIZE, 'fontweight', 'bold');
sgtitle('Fraction of Arena Bins Visited', 'fontsize', FONT_SIZE, 'fontweight', 'bold');
axis tight

linkaxes(ax, 'xy')

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S2B_trial_map_fraction_occupied_histogram', {'png', 'svg', 'fig'});

%% Cups

tetrodesTotal = length(TTetrodes.fractionCupsOccupied);
tetrodesSumAbove = sum(TTetrodes.fractionCupsOccupied > FRACTION);
tetrodesPercentAbove = tetrodesSumAbove / tetrodesTotal * 100;

calciumTotal = length(TCalcium.fractionCupsOccupied);
calciumSumAbove = sum(TCalcium.fractionCupsOccupied > FRACTION);
calciumPercentAbove = calciumSumAbove / calciumTotal * 100;


hFig = figure('position', get(0, 'screensize'));

ax = [];
ax(1) = subplot(1,2,1);
histogram(TTetrodes.fractionCupsOccupied, binEdges);
xlabel('Fraction of Cup Bins Visited')
ylabel('Counts (1 per trial)')
grid on
title(sprintf('Tetrodes: Total = %d, Above %0.1f =  %d, Percent Above = %0.1f', tetrodesTotal, FRACTION, tetrodesSumAbove, tetrodesPercentAbove))
set(gca, 'fontsize', FONT_SIZE, 'fontweight', 'bold');
axis tight

ax(2) = subplot(1,2,2);
histogram(TCalcium.fractionCupsOccupied, binEdges);
xlabel('Fraction of Cup Bins Visited')
ylabel('Counts (1 per trial)')
grid on
title(sprintf('Calcium: Total = %d, Above %0.1f =  %d, Percent Above = %0.1f', calciumTotal, FRACTION, calciumSumAbove, calciumPercentAbove))
set(gca, 'fontsize', FONT_SIZE, 'fontweight', 'bold');
sgtitle('Fraction of Cup Bins Visited', 'fontsize', FONT_SIZE, 'fontweight', 'bold');
axis tight

linkaxes(ax, 'xy')

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S2C_trial_map_cup_fraction_occupied_histogram', {'png', 'svg', 'fig'});


%%
function [SessionsToGroups] = load_sessions_to_groups_table(isCalcium)
    if isCalcium
        SessionsToGroups = load_calcium_sessions_to_groups_table();
    else
        SessionsToGroups = load_tetrode_sessions_to_groups_table();        
    end
    
end % function
