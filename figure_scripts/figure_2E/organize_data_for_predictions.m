% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024

function [R] = organize_data_for_predictions(dataset)
    AVERAGE_OVER_MULITPLE_FIELDS = false;
    PLACEFIELD_PERCENTILE_THRESHOLD = 95;
    R = [];

    MapsData = dataset.MapsData;
    BA = dataset.BestAligned;
    
    % get an example map. All maps of the same dataset need to have the
    % same dimensions.
    mExample = MapsData.map{1};
    mapWidth = size(mExample,2);
    mapHeight = size(mExample,1);
    clear mExample
    
    [sessions] = get_sessions_from_table(MapsData);
    for iSession = 1:size(sessions,1)
        animalName = sessions.animalName{iSession};
        sessionName = sessions.sessionName{iSession};
    
        fprintf('Processing %d of %d: %s %s\n', iSession, size(sessions,1), animalName, sessionName);
    
        SessionMapsData = MapsData(ismember(MapsData.animalName, animalName) & ismember(MapsData.sessionName, sessionName), :);
        BASession = BA(ismember(BA.animalName, animalName) & ismember(BA.sessionName, sessionName), :);
    
        % Meta data for the session and trials
        TrialsData = unique(SessionMapsData(:, ismember(SessionMapsData.Properties.VariableNames, {'animalName', 'sessionName', 'trialId', 'contextId', 'dig'})), 'rows');
        nTrials = size(TrialsData,1);
    
        % Get the cell names
        cellNames = unique(MapsData.cellName(ismember(MapsData.animalName, animalName) & ismember(MapsData.sessionName, sessionName)));
        nCells = length(cellNames);

        % Get the registered cell names
        IS_CALCIUM_DATA = any(ismember(MapsData.Properties.VariableNames, 'registeredCellName'));
        if IS_CALCIUM_DATA
            regCellNames = cell(nCells,1);
            for iCell = 1:nCells
                cellName = cellNames{iCell};
                imatch = find(ismember(MapsData.cellName, cellName), 1, 'first');
                if isempty(imatch)
                    error('No registered cell name found for %s\n', cellName);
                end
                regCellNames{iCell} = MapsData.registeredCellName{imatch};
            end
        end


        % Get the correlation values in the same order as the cell names
        bestCorrelations = nan(1,nCells);
        for iCell = 1:nCells
            cellName = cellNames{iCell};
            ind = find(ismember(BASession.cellName, cellName));
            if ~isempty(ind)
                bestCorrelations(iCell) = BASession.bestCorrelation(ind);
            end
        end
    
        TrialsData.maps = cell(nTrials,1);
        TrialsData.cellNames = cell(nTrials,1); % this will be the same for each row, but better
        TrialsData.centerOutAnglesDeg = cell(nTrials,1);
        TrialsData.bestCorrelations = cell(nTrials,1);

        if IS_CALCIUM_DATA
            TrialsData.regCellNames = cell(nTrials,1); % this will be the same for each row, but better
        end
    
        % For each trial, put the maps in order of the cell names. If a map
        % isn't present, then fill it with nan's so we can store in an array.
        for iTrial = 1:nTrials
            trialId = TrialsData.trialId(iTrial);
            %contextId = TrialsData.contextId(iTrial); not needed
            MapsTrialData = SessionMapsData(SessionMapsData.trialId == trialId,:);
    
            trialMaps = nan(mapHeight, mapWidth, nCells);
            for iCell = 1:nCells
                cellName = cellNames{iCell};
                iMatch = find(ismember(MapsTrialData.cellName, cellName));
                if ~isempty(iMatch)
                    trialMaps(:,:,iCell) = MapsTrialData.map{iMatch};
                end
            end

            
    
            % Now compute the center out angle
            trialCenterOutAngles = nan(nCells,1);
            trialRadialDistance = nan(nCells,1);
            for iCell = 1:size(trialMaps,3)
                m = squeeze(trialMaps(:,:,iCell));
                if all(isnan(m),'all')
                    continue; % skip if no map data
                end
    
                centerX = size(m,2)/2; % all maps are the same, but just in case
                centerY = size(m,1)/2;
                if AVERAGE_OVER_MULITPLE_FIELDS
                    [aAngle, ax, ay] = ml_alg_placemap_center_out_angle_mulitple_fields(m, PLACEFIELD_PERCENTILE_THRESHOLD);
                else
                    % This will just get the angle to the field with the largest
                    % area.
                    [aAngle, ax, ay] = ml_alg_placemap_center_out_angle(m, PLACEFIELD_PERCENTILE_THRESHOLD);            
                end
                trialRadialDistance(iCell) = sqrt((ax - centerX)^2 + (ay - centerY)^2);
                trialCenterOutAngles(iCell) = aAngle;
            end
    
    
            TrialsData.maps{iTrial} = trialMaps;
            TrialsData.centerOutAnglesDeg{iTrial} = trialCenterOutAngles;
            TrialsData.radialDistance{iTrial} = trialRadialDistance;
            TrialsData.cellNames{iTrial} = cellNames;
            TrialsData.bestCorrelations{iTrial} = bestCorrelations;
            if IS_CALCIUM_DATA
                TrialsData.regCellNames{iTrial} = regCellNames;
            end
        end
        
    
        if isempty(R)
            R = TrialsData;
        else
            R = [R; TrialsData];
        end
    end
end % function
