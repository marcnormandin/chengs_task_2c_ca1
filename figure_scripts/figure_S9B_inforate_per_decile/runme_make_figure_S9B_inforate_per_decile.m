% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This code makes Figure S9B (inforate per decile). It needs for Figure S9A
% to be run first to make required Excels.
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = '../figure_S9A/local_output';

OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

% Read in the information rates previously calculated for Figure S9A
calData = readtable(fullfile(INPUT_FOLDER, 'information_content_calcium_data.xlsx'));
tetData = readtable(fullfile(INPUT_FOLDER, 'information_content_tetrodes_data.xlsx'));

%%
T = [calData; tetData];

% Recalc the decile edges
nBins = 10; % deciles
edges = linspace(0,1,nBins+1);
p = prctile(T.bestCorrelation, edges(2:end-1)*100);
pe = [0, p, 1];
T.binId = discretize(T.bestCorrelation, [0, p, 1]);
T.binId(T.bestCorrelation < 0) = 1;
% Remove the nans
T(isnan(T.binId),:) = [];


%% Save the results for statistics to Excel
% Sort the table for Isabel
T = sortrows(T, {'binId', 'groupId', 'animalName'});
ofn = fullfile(LOCAL_OUTPUT_FOLDER, 'caltet_spatial_information_rate_per_decile.xlsx');
writetable(T, ofn);

%% Organize the data for plotting
TT = [];
for binId = 1:10
    ifr = T.ics_avg(T.binId == binId,:);
    % Remove outliers for plotting
    ifr = rmoutliers(ifr);
    X = nan(length(ifr),2);
    X(:,1) = ifr;
    X(:,2) = binId;
    TT = [TT; X];
end

%% Make the figure
hFig = figure();
boxplot(TT(:,1), TT(:,2))
xlabel('Decile')
ylabel('Avg Spatial Information Content (bits/s)')
title(sprintf('Average Spatial Information Content (bits/s)\n(Calcium + Tetrode Data, Days 1,2,3)\nExtreme Outliers Removed'))
grid on
set(gca, 'fontweight', 'bold', 'fontsize', 18);

mulana_savefig(hFig, OUTPUT_FOLDER, "figure_S9B_inforate_per_decile", {'png', 'svg'});
