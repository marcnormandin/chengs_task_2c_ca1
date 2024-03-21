% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [TMissing] = helper_mapsdata_get_missing_trials_as_zeromaps(T)
    % All maps should have this size
    mapSize = size(T.map{1});

    % Get all of the sessions present in the main map table
    [Sessions] = get_sessions_from_table(T);
    Sessions.numTrials = nan(size(Sessions,1),1);
    for iSession = 1:size(Sessions,1)
        animalName = Sessions.animalName{iSession};
        sessionName = Sessions.sessionName{iSession};
        TSession = T(ismember(T.animalName, animalName) & ismember(T.sessionName, sessionName),:);
        numTrials = max(TSession.trialId);
        Sessions.numTrials(iSession) = numTrials;
    end

    TMissing = [];
    % Select one of them
    for iSession = 1:size(Sessions,1)
        animalName = Sessions.animalName{iSession};
        sessionName = Sessions.sessionName{iSession};
        % Get the trial info
        SessionTrialInfo = get_session_trial_info(T, animalName, sessionName);
        % Get all the map data associated with the session
        TSession = T(ismember(T.animalName, animalName) & ismember(T.sessionName, sessionName),:);
        % Get all of the cell names present
        cellNames = unique(TSession.cellName);
        nCells = length(cellNames);
        for iCell = 1:nCells
            cellName = cellNames{iCell};
            TSessionCell = TSession(ismember(TSession.cellName, cellName),:);
            % Check if the cell has missing maps
            trialIdsPresent = TSessionCell.trialId(:);
            trialIdsWanted = SessionTrialInfo.trialId;
            trialIdsMissing = setdiff(trialIdsWanted, trialIdsPresent);
            if ~isempty(trialIdsMissing)
                for iMissing = 1:length(trialIdsMissing)
                    k = length(TMissing) + 1;
                    TMissing(k).animalName = animalName;
                    TMissing(k).sessionName = sessionName;
                    TMissing(k).dayNum = TSessionCell.dayNum(1);
                    TMissing(k).cellName = cellName;
                    TMissing(k).cellId = TSessionCell.cellId(1);
                    TMissing(k).animalSessionCellName = TSessionCell.animalSessionCellName{1};

                    trialIdMissing = trialIdsMissing(iMissing);
                    ind = find(SessionTrialInfo.trialId == trialIdMissing);

                    TMissing(k).trialId = trialIdMissing;
                    TMissing(k).contextId = SessionTrialInfo.contextId(trialIdMissing);
                    TMissing(k).map = zeros(mapSize);
                    TMissing(k).dig = SessionTrialInfo.dig{ind};
                    TMissing(k).mfr = nan;
                    TMissing(k).pfr = nan;
                    TMissing(k).infoRate = nan;
                    TMissing(k).infoSpike = nan;
                end
            end
        end
    end % iSession
    TMissing = struct2table(TMissing);
end % function



function [B] = get_session_trial_info(T, animalName, sessionName)
TSession = T(ismember(T.animalName, animalName) & ismember(T.sessionName, sessionName),:);
udigs = unique(TSession.dig);
digToIdMap = containers.Map(udigs, 1:length(udigs));
idToDigMap = containers.Map(1:length(udigs), udigs);
A = unique([TSession.trialId, TSession.contextId, TSession.dayNum, cellfun(@(x)digToIdMap(x), TSession.dig)], 'rows');
B = array2table(A(:,[1,2,3]), 'VariableNames', {'trialId', 'contextId', 'dayNum'});
digIds = A(:,4);

digs = arrayfun(@(x)sprintf('%s',idToDigMap(x)), digIds, 'UniformOutput', false);
B.dig = digs;
end % function

