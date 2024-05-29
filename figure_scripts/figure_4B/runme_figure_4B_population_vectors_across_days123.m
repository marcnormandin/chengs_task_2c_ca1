% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes the paper's Figure 4B, which are the cumulative popvector
% curves for across and stability for days 1,2 and 3.
% The figure was post-processed to match the style of the paper.
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

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calData = load_data(calciumAnalysisResultsFolder);

%%
% We need to override the output folder so we don't contaminate anything.
calData.analysisSettings.OUTPUT_FOLDER = OUTPUT_FOLDER;
if ~exist(calData.analysisSettings.OUTPUT_FOLDER,'dir')
    mkdir(calData.analysisSettings.OUTPUT_FOLDER);
end

calData.analysisSettings.OUTPUT_RESULTS_FOLDER = OUTPUT_FOLDER;
if ~exist(calData.analysisSettings.OUTPUT_RESULTS_FOLDER,'dir')
    mkdir(calData.analysisSettings.OUTPUT_RESULTS_FOLDER);
end

analysisSettings = calData.analysisSettings;
analysisInput = calData.analysisInput;
analysisResults = calData.analysisResults;

%% This is the function that was used to generate the pop Figure 4B
F = compute_popvectors_across_contexts_unregistered_figure_4B(analysisSettings, analysisResults.BestAligned, analysisSettings.POPVECTORS_BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, OUTPUT_FOLDER);

%% Code added to save data in Excel as required by the journal
nRows = size(F,1);
A = zeros(651, nRows);
columnNames = cell(nRows,1);
for iRow = 1:nRows
    A(:,iRow) = F.dp_across{iRow};

    dayString = strrep(F.sessionName{iRow}, ' ', '_');
    if F.isStable(iRow) == true
        stabilityString = 'FI';
    else
        stabilityString = 'FS';
    end
    columnNames{iRow} = sprintf('dot_product_across_%s_%s', dayString, stabilityString);
end % iRow
T = array2table(A, 'VariableNames', columnNames);

writetable(T, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_4B.xlsx'));

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
    [analysisInput] = load_analysis_input(folder, analysisSettings);

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
