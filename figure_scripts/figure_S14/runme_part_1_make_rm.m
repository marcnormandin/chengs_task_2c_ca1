% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This part makes the "*_rm.mat" files needed by the other scripts which takes FOREVER.
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;


calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');

SPEED_THRESHOLD_CM_PER_S = 2;
SPIKE_THRESHOLD = 0;
NUM_TRIALS = 12;
RM_ROWS = [1:2:NUM_TRIALS, 2:2:NUM_TRIALS];

OUTPUT_FOLDER = "./local_output";
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

tmp = load(fullfile(calciumAnalysisResultsFolder, 'analysis_input.mat'));
MapsData = tmp.analysisInput.MapsData;
ANIMAL_NAMES ={'CMG129_CA1', 'CMG154_CA1', 'CMG161_CA1', 'CMG162_CA1', 'CMG169_CA1'};

for iDataset = 1:length(ANIMAL_NAMES)
    ANIMAL_NAME = ANIMAL_NAMES{iDataset};
    fn = fullfile(INPUT_FOLDER, 'calcium_storage', sprintf('%s.mat', ANIMAL_NAME));
    data = load(fn);

    movement = data.acfg.MovementData;
    ND = data.acfg.NeuralData;
    SD = data.acfg.ScopeData;

    T = [];
    k = 1;
    for iRow = 1:size(ND,1)
        fprintf('Processing %d/%d\n', iRow, size(ND,1));
        nd = table2struct(ND(iRow,:));
        % Get the timestamps associated with the row
        matchScope = table2struct(SD(ismember(SD.animalName, nd.animalName) & ismember(SD.sessionName, nd.sessionName) & SD.trialId == nd.trialId,:));
        matchMovement = table2struct(movement(ismember(movement.animalName, nd.animalName) & ismember(movement.sessionName, nd.sessionName) & movement.trialId == nd.trialId,:));
        
        tScope_ms = matchScope.timestamp_ms;
        tScope_ms = tScope_ms - tScope_ms(1);

        movement_t_ms = matchMovement.pos_t_ms;
        movement_t_ms = movement_t_ms - movement_t_ms(1);

        speed = interp1(movement_t_ms, matchMovement.smoothed_speed_cm_s, tScope_ms);
        speedInds = find(speed > SPEED_THRESHOLD_CM_PER_S);
        
        dt_s = median(diff(tScope_ms))/1000;
        
        % Use normalized inferred spikes
        y = nd.spikes;
        spikeMax = max(y, [], 'all');
        y = y ./ spikeMax;

        yInds = find(y>SPIKE_THRESHOLD);
        
        yValidInds = intersect(yInds, speedInds);
        
        m1 = sum(y(yValidInds)) / (length(yValidInds) * dt_s);
    
        % binary
        m2 = length(yValidInds) / (length(speedInds) * dt_s);
        
        % magnitude 
        m3 = sum(y(yValidInds)) / (length(speedInds) * dt_s);

        % maximum rate
        m4 = max(nd.spikes(yValidInds), [], 'all', 'omitnan') / dt_s;
        if isempty(m4)
            m4 = nan;
        end
    
        % store
        T(k).groupName = nd.groupName;
        T(k).animalName = nd.animalName;
        T(k).sessionName = nd.sessionName;
        T(k).trialId = nd.trialId;
        T(k).cellId = nd.cellId;
        T(k).spikeMax = spikeMax;
        T(k).m1 = m1;
        T(k).m2 = m2;
        T(k).m3 = m3;
        T(k).m4 = m4;
        k = k + 1;
    end % iRow
    T = struct2table(T);

    [sessions] = get_sessions_and_cellnames_from_table(T);
    
    RM = [];
    kRM = 1;
    
    for iRow = 1:size(sessions, 1)
        fprintf('Processing RM %d/%d\n', iRow, size(sessions,1));
        s = table2struct(sessions(iRow,:));
    
        matches = T(ismember(T.animalName, s.animalName) & ismember(T.sessionName, s.sessionName) & T.cellId == s.cellId,:);
        isExcluded = true(size(matches,1),1);
        for iMatch = 1:size(matches,1)
            findInd = find(ismember(MapsData.animalName, matches.animalName{iMatch}) & ismember(MapsData.sessionName, matches.sessionName{iMatch}) ...
                & MapsData.trialId == matches.trialId(iMatch) & MapsData.cellId == matches.cellId(iMatch));
            if length(findInd) == 1
                if MapsData.passedCelllikeSpatialFootprintFilter(findInd)
                    isExcluded(iMatch) = false;
                end
            end
        end
    
        matches(isExcluded,:) = [];
    
        if size(matches,1) <= 2
            continue;
        end
        
        RM(kRM).groupName = matches.groupName{1};
        RM(kRM).animalName = matches.animalName{1};
        RM(kRM).sessionName = matches.sessionName{1};
        RM(kRM).cellId = matches.cellId(1);
    
        fieldNames = {'m1', 'm2', 'm3', 'm4'};
        for iType = 1:length(fieldNames)
            fieldName = fieldNames{iType};
            rm = nan(NUM_TRIALS, NUM_TRIALS);
            for iMatchA = 1:size(matches,1)
                trialIdA = matches.trialId(iMatchA);
                for iMatchB = 1:size(matches,1)
                    trialIdB = matches.trialId(iMatchB);
                    if trialIdA == trialIdB
                        continue;
                    end
    
                    % 
                    d = abs(matches.(fieldName)(iMatchA) - matches.(fieldName)(iMatchB));
                    tA = find(RM_ROWS == trialIdA);
                    tB = find(RM_ROWS == trialIdB);
                    rm(tA, tB) = d;
                    rm(tB, tA) = d; % overlaped...
                end
            end
    
            RM(kRM).(fieldName) = rm;
        end % iType
        kRM = kRM + 1;
    end % iRow
    RM = struct2table(RM);

    RMA = [];
    kRMA = 1;
    [sessions] = get_sessions_from_table(RM);
    for iSession = 1:size(sessions,1)
        s = table2struct(sessions(iSession,:));
    
        matches = RM( ismember(RM.animalName, s.animalName) & ismember(RM.sessionName, s.sessionName),:);
    
        fieldNames = {'m1', 'm2', 'm3', 'm4'};
        RMA(kRMA).animalName = s.animalName;
        RMA(kRMA).sessionName = s.sessionName;
    
        for iField = 1:length(fieldNames)
            fieldName = fieldNames{iField};
            rma = nan(NUM_TRIALS, NUM_TRIALS, size(matches,1));
            for iMatch = 1:size(matches,1)
                rma(:,:,iMatch) = matches.(fieldName){iMatch};
            end
            rma = mean(rma, 3, 'omitnan');
            
            RMA(kRMA).(fieldName) = rma;
        end
        kRMA = kRMA + 1;
    end
    RMA = struct2table(RMA);


    close all
    
    for iFieldName = 1:length(fieldNames)
        fieldName = fieldNames{iFieldName};
    
        hFig = figure('position', get(0, 'screensize'));
        for iSession = 1:size(RMA,1)
            subplot(1,size(RMA,1),iSession)
            imagesc(RMA.(fieldName){iSession});
            title(sprintf('%s %s', fieldName, RMA.sessionName{iSession}))
            axis equal tight
            xticks(1:NUM_TRIALS)
            xticklabels(arrayfun(@(x)(num2str(x)), RM_ROWS, 'UniformOutput', false))
            yticks(1:NUM_TRIALS)
            yticklabels(arrayfun(@(x)(num2str(x)), RM_ROWS, 'UniformOutput', false))
        end
        mulana_savefig(hFig, OUTPUT_FOLDER, sprintf('%s_rm_%s', ANIMAL_NAME, fieldName), {'png', 'svg'})
    end

    save(fullfile(OUTPUT_FOLDER, sprintf('%s_rm.mat', ANIMAL_NAME)), 'RM', 'RMA', 'NUM_TRIALS', 'RM_ROWS', 'SPEED_THRESHOLD_CM_PER_S', 'SPIKE_THRESHOLD');
end % iDataset

function [sessions] = get_sessions_and_cellnames_from_table(T)
    % Returns a unique table of animal name / session name pairs
    animalIndex = find(ismember(T.Properties.VariableNames, 'animalName'));
    sessionIndex = find(ismember(T.Properties.VariableNames, 'sessionName'));
    cellIdIndex = find(ismember(T.Properties.VariableNames, 'cellId'));
    sessions = unique(T(:,[animalIndex, sessionIndex, cellIdIndex]), 'rows');
end
