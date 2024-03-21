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
    mkdir(LOCAL_OUTPUT_FOLDER);
end

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');

fprintf('Loading the calcium dataset\n');
calDataset = load_dataset(calciumAnalysisResultsFolder);
fprintf('Loading the tetrode dataset\n');
tetDataset = load_dataset(tetrodesAnalysisResultsFolder);
datasets = {calDataset, tetDataset};

analysisSettings = cell(length(datasets),1);
for iDataset = 1:length(datasets)
    analysisSettings{iDataset} = datasets{iDataset}.analysisSettings;
end

%% Get the data to use for predictions
fprintf('Creating the prediction datasets\n');
predictionDatasets = cell(length(datasets),1);
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};
    predictionDatasets{iDataset} = organize_data_for_predictions(dataset);

    [SessionsToGroups] = load_sessions_to_groups_table(dataset.analysisSettings);
    [predictionDatasets{iDataset}] = helper_add_group_labels_to_table(predictionDatasets{iDataset}, SessionsToGroups, true);
end

%%
close all
clc

algoSettings = [];
algoSettings.kernelFunction = 'polynomial';
algoSettings.MIN_TRIALS = -inf;
algoSettings.MIN_CELLS = 2;
algoSettings.digTypes = {};

% Create a table for different predictions to be done
PS = [make_prediction_settings('Geo', 'Corr', true); ...
     make_prediction_settings('Wrong', 'Feat', false)];

for iRow = 1:size(PS,1)
    params = table2struct(PS(iRow,:));
    algoSettings.digTypes = {params.DigTypeA, params.DigTypeB};
    algoSettings.makeFigure = params.makeFigure;

    fprintf('Running prediction analysis\n');
    analysisResultsNoStability = process_allcells(predictionDatasets, analysisSettings, algoSettings);

    % Filter
    analysisResultsNoStability.R = analysisResultsNoStability.R0;
    analysisResultsNoStability.R(analysisResultsNoStability.R.nCellsUsed < algoSettings.MIN_CELLS,:) = [];
    analysisResultsNoStability.S = sortrows(analysisResultsNoStability.R, {'groupLabel', 'isCalcium'});

    % Make a string to embedd the dig types in the filenames.
    digTypesStr = join(algoSettings.digTypes, '_');
    digTypesStr = digTypesStr{1};

    % Save the accuracies (no stability)
    S = analysisResultsNoStability.S;
    S(~ismember(S.groupLabel, {'Day 1', 'Day 2', 'Day 3'}), :) = [];
    S.nDigs = cellfun(@(x)(length(x)), S.digIds);
    columnsToSave = {'animalName', 'sessionName', 'groupLabel', 'groupId', 'accuracy', 'nDigs', 'nCellsUsed', 'isCalcium'};
    S(:, ~ismember(S.Properties.VariableNames, columnsToSave)) = [];
    resultsFn = fullfile(LOCAL_OUTPUT_FOLDER, sprintf('train_predict_%s_same_day_using_centeroutangledeg_nostability_accuracies.xlsx', digTypesStr));
    writetable(S, resultsFn);
    fprintf('Data saved to: %s\n', resultsFn);

    % Make a figure (optional)
    if algoSettings.makeFigure == true
        [hFig] = make_figure(resultsFn);
        mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_2E_predictions', {'svg', 'png'});
    end
end

% Helper function
function T = make_prediction_settings(digTypeA, digTypeB, doMakeFigure)
    T = table({digTypeA}, {digTypeB}, doMakeFigure, 'VariableNames',{'DigTypeA', 'DigTypeB', 'makeFigure'});
end


%% Any Stability
function [analysisResults] = process_allcells(predictionDatasets, analysisSettings, algoSettings)
    R0 = compute_R_anycontext_anystability(predictionDatasets, analysisSettings, algoSettings);
    
    R = R0;
    R(R.nCellsUsed < algoSettings.MIN_CELLS,:) = [];

    analysisResults = [];
    analysisResults.R0 = R0;
    analysisResults.R = R;
    analysisResults.MIN_TRIALS = algoSettings.MIN_TRIALS;
    analysisResults.digTypes = algoSettings.digTypes;
    analysisResults.MIN_CELLS = algoSettings.MIN_CELLS;
end % function




%% FUNCTIONS
function [R] = compute_R_anycontext_anystability(predictionDatasets, analysisSettings, algoSettings)
    R = [];
    rk = 1;
    for iDataset = 1:length(predictionDatasets)
        predictionDataset = predictionDatasets{iDataset};
        
        [sessions] = get_sessions_from_table(predictionDataset);
        [SessionsToGroups] = load_sessions_to_groups_table(analysisSettings{iDataset});
        [sessions] = helper_add_group_labels_to_table(sessions, SessionsToGroups, true);
        badSessionInds = [];
        for iSession = 1:size(sessions,1)
            animalName = sessions.animalName{iSession};
            sessionName = sessions.sessionName{iSession};
            fprintf('Processing session data %d of %d: %s %s\n', iSession, size(sessions,1), animalName, sessionName);
        
            try
                [digIds, predictedDigIds, accuracy, nCellsUsed] = predict_session(predictionDataset, animalName, sessionName, algoSettings);
            
                R(rk).animalName = animalName;
                R(rk).sessionName = sessionName;
                R(rk).groupLabel = sessions.groupLabel{iSession};
                R(rk).groupId = sessions.groupId(iSession);
                R(rk).digIds = digIds;
                R(rk).predictedDigIds = predictedDigIds;
                R(rk).accuracy = accuracy;
                R(rk).nCellsUsed = nCellsUsed;
                R(rk).isCalcium = analysisSettings{iDataset}.IS_CALCIUM_DATA;
                rk = rk + 1;
            catch e
                badSessionInds(end+1) = iSession;
            end
        end
        disp(badSessionInds)
    end
    R = struct2table(R);
end % function

function [digIds, predictedDigIds, accuracy, nCellsUsed] = predict_session(predictionDataset, animalName, sessionName, algoSettings)
    MIN_TRIALS = algoSettings.MIN_TRIALS;
    digTypes = algoSettings.digTypes;

    predictionDataset(~ismember(predictionDataset.dig, digTypes),:) = [];

    sessionData = predictionDataset(ismember(predictionDataset.animalName, animalName) & ismember(predictionDataset.sessionName, sessionName), :);
    
    % We need a set of cells that are present for all of the trials
    cellNamesAll = sessionData.cellNames{1};
    nTrials = size(sessionData,1);

    if nTrials < MIN_TRIALS
        error('Less than %d trials', MIN_TRIALS)
    end

    nCells = length(cellNamesAll);
    centerOutAnglesDeg = nan(nTrials, nCells);

    for iTrial = 1:nTrials
        centerOutAnglesDeg(iTrial, :) = sessionData.centerOutAnglesDeg{iTrial};
    end
    badCellInds = find(any(isnan(centerOutAnglesDeg),1));

    centerOutAnglesDeg(:, badCellInds) = [];

    cellNames = cellNamesAll;
    cellNames(badCellInds) = [];

    nCellsUsed = length(cellNames);
    
    digMap = containers.Map();
    for iDig = 1:length(digTypes)
        digMap(digTypes{iDig}) = iDig;
    end
    digIds = nan(nTrials,1);
    for iTrial = 1:nTrials
        digIds(iTrial) = digMap(sessionData.dig{iTrial});
    end

    X = [sind(centerOutAnglesDeg), cosd(centerOutAnglesDeg)];
    Y = digIds;
    predictedDigIds = zeros(size(digIds));
    for iTrial = 1:nTrials
        fprintf('Predicting trial %d of %d\n', iTrial, nTrials);

        XTrain = X(setdiff(1:nTrials, iTrial),:);
        YTrain = Y(setdiff(1:nTrials, iTrial));

        SVMModel = fitcsvm(XTrain,YTrain,'Standardize',false,'KernelFunction',algoSettings.kernelFunction, 'KernelScale','auto');
    
        [yPredict,~] = predict(SVMModel,X(iTrial,:));
    
        predictedDigIds(iTrial) = yPredict;

    end % iTrial

    accuracy = sum(digIds == predictedDigIds) / nTrials;
end % function
