% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This code computes the population vector cumulative distributions per
% decile (Figure S12A).
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

% To save the excel for further processing
LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

% These are the folders created by the main analysis script to store the
% results
calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calData = load_data(calciumAnalysisResultsFolder);

%%
% We need to override the output folder so we don't possibly contaminate anything.
calData.analysisSettings.OUTPUT_FOLDER = LOCAL_OUTPUT_FOLDER;
if ~exist(calData.analysisSettings.OUTPUT_FOLDER,'dir')
    mkdir(calData.analysisSettings.OUTPUT_FOLDER);
end

calData.analysisSettings.OUTPUT_RESULTS_FOLDER = LOCAL_OUTPUT_FOLDER;
if ~exist(calData.analysisSettings.OUTPUT_RESULTS_FOLDER,'dir')
    mkdir(calData.analysisSettings.OUTPUT_RESULTS_FOLDER);
end

analysisSettings = calData.analysisSettings;
analysisInput = calData.analysisInput;
analysisResults = calData.analysisResults;

%% This is the function that was used to generate the pop Figure 4B
BestAligned = analysisResults.BestAligned;
nBins = 10; % deciles
tmpEdges = linspace(0,1,nBins+1);
p = prctile(BestAligned.bestCorrelation, tmpEdges(2:end-1)*100);
pe = [0, p, 1];
binIds = discretize(BestAligned.bestCorrelation, [0, p, 1]);
binIds(BestAligned.bestCorrelation < 0) = 1;
pe(1) = -inf;
% pe contains the edges we use for the best aligned correlation deciles
bounds = nan(nBins,2);
bounds(1:nBins, 1) = pe(1:end-1);
bounds(1:nBins, 2) = pe(2:end);
nBounds = size(bounds,1); % same as nBins


%%
F = [];
for iBound = 1:nBounds
    bestAlignedCorrelationMin = bounds(iBound,1);
    bestAlignedCorrelationMax = bounds(iBound,2);
    FB = compute_popvectors_across_contexts_unregistered_v3(analysisSettings, analysisResults.BestAligned, bestAlignedCorrelationMin, bestAlignedCorrelationMax);
    FB.iBound = ones(size(FB,1),1) * iBound;
    FB.bestAlignedCorrelationMin = ones(size(FB,1),1) * bestAlignedCorrelationMin;
    FB.bestAlignedCorrelationMax = ones(size(FB,1),1) * bestAlignedCorrelationMax;
    F = [F; FB];
end

%% Make the figure
hFig = figure('position', get(0, 'screensize'));
PP = 1;
QQ = 3; % 3 subplots
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroups = length(groupLabels);
for iGroup = 1:nGroups
    ll = {};
    for iBound = 1:nBounds
        ax(iGroup) = subplot(PP,QQ,iGroup);
        hold on
        groupLabel = groupLabels{iGroup};
        indRow = ismember(F.sessionName, groupLabel) & F.iBound == iBound;
        xData = F.cumdist_x{indRow};
        yData = F.cumdist_y{indRow};
       
        plot(xData, yData, 'linewidth', 2)
        xlabel('Normalized Dot Product')
        ylabel('Cumulative')
        legend(groupLabels, 'location', 'southoutside', 'orientation', 'horizontal');
        title(sprintf('%s', groupLabel));
        ll{end+1} = sprintf('Decile %d', iBound);
        grid on
        set(gca, 'fontweight', 'bold', 'fontsize', 18)
    end % iBound
    legend(ll, 'location', 'eastoutside', 'orientation', 'vertical');
end % iGroup

for i = 1:length(ax)
    set(0, 'currentfigure', hFig);
    set(hFig, 'currentaxes', ax(i));
    axis square
end

linkaxes(ax, 'xy')
sgtitle(sprintf('Population Similarity Across Context'))
mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S12A', {'png', 'svg'});


%% Added to export data for natcomms
NC = nan(651, size(F,1)); % each array is size 651
k = 1;
variableNames = cell(1,size(NC,2));
for iRow = 1:size(F,1)
    NC(:,k) = F.dp_across{iRow};
    variableNames{k} = sprintf('dot_product_across_%s_decile_%d', strrep(F.sessionName{iRow}, ' ', '_'), F.iBound(iRow));
    k = k + 1;
end % iGroup
NC = array2table(NC, 'VariableNames', variableNames)
writetable(NC, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_S12A.xlsx'), 'Sheet', 'figure_S12A');


%% Loading functions
function [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned)
    n_errors = 0;
    bad_ba_rows = [];
    for bestAligned_row_index = 1:size(BestAligned,1)
        animalName = BestAligned.animalName{bestAligned_row_index};
        sessionName = BestAligned.sessionName{bestAligned_row_index};
        cellName = BestAligned.cellName{bestAligned_row_index};
        bestAligned_cellDataInds = BestAligned.cellDataInds{bestAligned_row_index};
        [output] = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        mapsData_cellDataInds = output.cellDataInds;
        n_inds = length(bestAligned_cellDataInds);
        if n_inds ~= length(intersect(bestAligned_cellDataInds, mapsData_cellDataInds))
            %fprintf('Error with best aligned row index: %d\n', bestAligned_row_index);
            n_errors = n_errors + 1;
            bad_ba_rows(n_errors) = bestAligned_row_index;
        end
    end % bestAligned_row_index

end


function [s] = load_data(folder)
    [analysisResults, BestAligned, analysisSettings] = load_best_aligned(folder);
    %[analysisInput] = 
    [analysisInput] = load_analysis_input(folder, analysisSettings);

    % Verify the integrity of the data because we have to index between the tables.
%     [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned);
%     if n_errors > 0
%         error('Data compromised.\n');
%     end

    % Store it in a struct to make it easier to manage
    s = [];
    s.folder = folder;
    s.analysisSettings = analysisSettings;
    s.analysisResults = analysisResults;
    s.analysisInput = analysisInput;
    s.BestAligned = BestAligned;
end

function [analysisInput] = load_analysis_input(folder, analysisSettings)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    analysisInput = tmp.analysisInput;
end % function

function [analysisResults, BestAligned, analysisSettings] = load_best_aligned(folder)
    tmp = load(fullfile(folder, 'analysis_results.mat'));
    analysisSettings = tmp.analysisSettings;
    analysisResults = tmp.analysisResults;
    BestAligned = tmp.analysisResults.BestAligned;
    eliminateUnmatched = false;
    BestAligned = helper_add_group_labels_to_table(BestAligned, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function
