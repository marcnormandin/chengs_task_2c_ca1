% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure 6C.
%
% Note that the figures were post-processed in Adobe Illustator.
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

LOCAL_OUTPUT_FOLDER = './local_output'; % for Figure S5B
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calDataset = load_dataset(calciumAnalysisResultsFolder);
datasets = {calDataset};

%%
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
clc

params = [];
params.digTypes = {'Geo', 'Corr'};
params.stabilityType = 'stable'; % 'any', 'unstable', 'stable'
params.stabilityThreshhold = 0.3;
params.minCells = 2;
params.kernelFunction = 'polynomial';

dayPairs = nchoosek({'Day 1', 'Day 2', 'Day 3'}, 2);
nPairs = size(dayPairs,1);
RAll = [];
TAll = [];
for iPair = 1:nPairs
    params.trainDay = dayPairs{iPair,1};
    params.predictDay = dayPairs{iPair,2};
    
    R = compute_R_anycontext(predictionDatasets, analysisSettings, params);
    
    % Day Pair Totals
    x = R.accuracy * 100;
    nX = sum(isfinite(x)); % since we omit nan
    meanX = mean(x, 'all', 'omitnan');
    stdX = std(x, 1, 'all', 'omitnan');
    errX = stdX ./ sqrt(nX);

    tot = [];
    tot.trainGroupLabel = params.trainDay;
    tot.predictGroupLabel = params.predictDay;
    tot.nX = nX;
    tot.meanX = meanX;
    tot.stdX = stdX;
    tot.errX = errX;
    TAll = [TAll; struct2table(tot)];

    RAll = [RAll; R];

end % iPair

% Save the data
digsStr = join(params.digTypes, '_');
digsStr = digsStr{1};
RSave = RAll(:,ismember(RAll.Properties.VariableNames, {'animalName', 'trainGroupLabel', 'predictGroupLabel', 'accuracy', 'nPredictTrials', 'nPredictedTrials', 'isCalcium'}));
writetable(RSave, fullfile(LOCAL_OUTPUT_FOLDER, sprintf('train_predict_different_days_cells_%s_digs_%s_kernel_%s_minCells_%d.xlsx', params.stabilityType, digsStr, params.kernelFunction, params.minCells)));

%% Make the figure
clc
close all

hFig = figure();
bar(1:size(TAll,1), TAll.meanX)
hold on
errorbar(1:size(TAll,1), TAll.meanX, TAll.errX, TAll.errX, 'color', 'k', 'linestyle', 'none')
a = axis;
axis([a(1), a(2), 0, 100])
set(gca, 'fontweight', 'bold', 'fontsize', 18)
grid on
ylabel('Prediction Accuracy');
labels = cell(size(TAll,1),1);
for iRow = 1:size(TAll,1)
    l = sprintf('Train %s Predict %s', TAll.trainGroupLabel(iRow,:), TAll.predictGroupLabel(iRow,:));
    labels{iRow} = l;
end
xticks(1:size(TAll,1));
xticklabels(labels)
title('FI cell heading prediction over time')

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_6C', {'png', 'svg'});


%%
function [R] = compute_R_anycontext(predictionDatasets, analysisSettings, params)
    R = [];
    rk = 1;
    for iDataset = 1:length(predictionDatasets)
        predictionDataset = predictionDatasets{iDataset};
        
        [sessions] = get_sessions_from_table(predictionDataset);
        [SessionsToGroups] = load_sessions_to_groups_table(analysisSettings{iDataset});
        [sessions] = helper_add_group_labels_to_table(sessions, SessionsToGroups, true);
        badAnimalInds = [];
        animalNames = unique(sessions.animalName);
        nAnimals = length(animalNames);

        for iAnimal = 1:nAnimals
            animalName = animalNames{iAnimal};

            trainGroupLabel = params.trainDay;
            predictGroupLabel = params.predictDay;
            
            fprintf('\nProcessing %s %s %s\n', animalName, trainGroupLabel, predictGroupLabel);

            indTrain = find(ismember(sessions.animalName, animalName) & ismember(sessions.groupLabel, trainGroupLabel));
            indPredict = find(ismember(sessions.animalName, animalName) & ismember(sessions.groupLabel, predictGroupLabel));
            
            try
                trainSessionName = sessions.sessionName{indTrain};
                predictSessionName = sessions.sessionName{indPredict};

                [digIds, predictedDigIds, accuracy, scores, nPredictTrials, nPredictedTrials] = predict_different_day(predictionDataset, animalName, trainSessionName, predictSessionName, params);
            
                R(rk).animalName = animalName;
                R(rk).trainGroupLabel = trainGroupLabel;
                R(rk).predictGroupLabel = predictGroupLabel;
                R(rk).digIds = digIds;
                R(rk).predictedDigIds = predictedDigIds;
                R(rk).accuracy = accuracy;
                R(rk).nPredictTrials = nPredictTrials;
                R(rk).nPredictedTrials = nPredictedTrials;
                R(rk).scores = scores;
                R(rk).isCalcium = analysisSettings{iDataset}.IS_CALCIUM_DATA;
                rk = rk + 1;
            catch e
                fprintf('Error processing %s: %s\n%s', animalName, e.message, e.identifier);
                badAnimalInds(end+1) = iAnimal;
            end
        end
        %disp(badAnimalInds)
    end
    R = struct2table(R);
end % function

function [sessionData] = filter_sessiondata(sessionData, commonRegCellNames)
    nCommonRegCellNames = length(commonRegCellNames);
    
    regCellNames = sessionData.regCellNames{1}; % all trials have the same set of names
    cellNames = sessionData.cellNames{1}; % all trials have the same set of names
    commonCellNames = cell(nCommonRegCellNames,1);
    commonCellInds = nan(nCommonRegCellNames,1);
    for iCommon = 1:nCommonRegCellNames
        commonRegCellName = commonRegCellNames{iCommon};
        ind = find(ismember(regCellNames, commonRegCellName));
        if isempty(ind) || length(ind) ~= 1
            error('This should be impossible.');
        end
        % same ind is for the session cell name
        commonCellInds(iCommon) = ind;
        commonCellNames{iCommon} = cellNames{ind};
    end
    % Now that we have the inds relative to the session, we can filter the data
    
    nTrials = size(sessionData,1);
    for iTrial = 1:nTrials
        % get the originals
        %cellNames = sessionData.cellNames{iTrial};
        %regCellNames = sessionData.regCellNames{iTrial};
        maps = sessionData.maps{iTrial};
        centerOutAnglesDeg = sessionData.centerOutAnglesDeg{iTrial};
        bestCorrelations = sessionData.bestCorrelations{iTrial};
    
        % filter
        commonMaps = maps(:,:,commonCellInds);
        commonCenterOutAnglesDeg = centerOutAnglesDeg(commonCellInds);
        commonBestCorrelations = bestCorrelations(commonCellInds);
    
        % store
        sessionData.maps{iTrial} = commonMaps;
        sessionData.centerOutAnglesDeg{iTrial} = commonCenterOutAnglesDeg;
        sessionData.bestCorrelations{iTrial} = commonBestCorrelations;
    
        sessionData.cellNames{iTrial} = commonCellNames;
        sessionData.regCellNames{iTrial} = commonRegCellNames;
    end
end % function


function [X, Y, badCellInds] = prepare_svm_data(sessionData, params)
    nTrials = size(sessionData,1);

    nCells = length(sessionData.regCellNames{1});
    popData = nan(nTrials, nCells);
    for iTrial = 1:nTrials
        popData(iTrial, :) = sessionData.centerOutAnglesDeg{iTrial};
    end
    
    % find bad cell data (we need to eliminate the same cells from both
    % training and prediction so we need to filter AFTER, not now)
    badCellInds = find(any(isnan(popData),1));
    
    % map the strings to numbers
    digMap = containers.Map();
    for iDig = 1:length(params.digTypes)
        digMap(params.digTypes{iDig}) = iDig;
    end % iDig
    digIds = nan(nTrials,1);
    for iTrial = 1:nTrials
        digIds(iTrial) = digMap(sessionData.dig{iTrial});
    end

    X = popData;
    Y = digIds;
end

function [digIds, predictedDigIds, accuracy, scores, nPredictTrials, nPredictedTrials] = predict_different_day(predictionDataset, animalName, sessionNameTrain, sessionNamePredict, params)
    % We only want trials that had the desired dig types (each row is a
    % session trial)
    predictionDataset(~ismember(predictionDataset.dig, params.digTypes),:) = [];

    % Get the table data, which we will filter and sort after
    sessionDataTrain = predictionDataset(ismember(predictionDataset.animalName, animalName) & ismember(predictionDataset.sessionName, sessionNameTrain), :);
    sessionDataPredict = predictionDataset(ismember(predictionDataset.animalName, animalName) & ismember(predictionDataset.sessionName, sessionNamePredict), :);
    
    sessionDataTrain = filter_prediction_dataset(sessionDataTrain, params);


    % We need a common set of ordered names from train and predict. Each
    % trial from the same session has the same sequence of names, so we
    % just need one of them (from any rows).
    commonRegCellNames = intersect(sessionDataTrain.regCellNames{1}, sessionDataPredict.regCellNames{1});
    
    % Now keep only data from the common reg cell names (cell names must be
    % common to the training day and the prediction day). Overwrite.
    sessionDataTrain = filter_sessiondata(sessionDataTrain, commonRegCellNames);
    sessionDataPredict = filter_sessiondata(sessionDataPredict, commonRegCellNames);

    % Training data
    [XTrainFull, YTrainFull, badCellIndsTrainFull] = prepare_svm_data(sessionDataTrain, params);

    % A hack to get the dig ids
    [~, digIds, ~] = prepare_svm_data(sessionDataPredict, params);
    
    nPredictTrials = size(sessionDataPredict,1);
    predictedDigIds = nan(nPredictTrials,1);
    scores = nan(nPredictTrials,2);
    for iTrialPredict = 1:nPredictTrials % predict every trial on the prediction day
        % Predict data
        sessionDataPredictTrial = sessionDataPredict(iTrialPredict,:);
        [XPredict, ~, badCellIndsPredict] = prepare_svm_data(sessionDataPredictTrial, params);
        
        % find cells that are invalid in either training or predicting, so we
        % dont use them
        badCellInds = union(badCellIndsTrainFull, badCellIndsPredict);
    
        % Make a copy, and then we eliminate the invalid cells.
        XTrain = XTrainFull;
        YTrain = YTrainFull; % added

        XTrain(:,badCellInds) = [];
        XPredict(:,badCellInds) = [];

        fprintf('Predicting trial %d of %d\n', iTrialPredict, nPredictTrials);

        % At this point there may be no cells to be used
        if isempty(XTrain) || isempty(XPredict)
            continue;
        end
    
        % NOW since data has been filtered, apply the sin and cos
        XTrainSinCos = [cosd(XTrain), sind(XTrain)];
        XPredictSinCos = [cosd(XPredict), sind(XPredict)];

        SVMModel = fitcsvm(XTrainSinCos, YTrain,'Standardize',true,'KernelFunction',params.kernelFunction,...
        'KernelScale','auto', 'Prior', 'uniform');
    
        [yPredict,score] = predict(SVMModel,XPredictSinCos);
    
        predictedDigIds(iTrialPredict) = yPredict;
        scores(iTrialPredict,:) = score;
    
    end % iTrialPredict

    % Compute accuracy based on the trials that we could use
    validPredictionInds = find(isfinite(predictedDigIds));
    nValidPredictions = length(validPredictionInds);
    
    trueDigIds = digIds(validPredictionInds);
    validPredictedDigIds = predictedDigIds(validPredictionInds);

    accuracy = sum(trueDigIds == validPredictedDigIds) / nValidPredictions;

    nPredictedTrials = nValidPredictions;
end % function
