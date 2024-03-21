% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R, success] = compute_normal_bfo90_table(PerCell, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION)

    BFO90_ROTATIONS_TO_CHECK = [0, 90, 180, 270];

    fieldSets = helper_load_bfo90_fieldsets_to_record();

    % Get the unique sessions for which we will compute one BF90 each
    SessionTable = get_sessions_from_table(PerCell);
    numSessions = size(SessionTable,1);

    R = [];
    success = true;
    wb = [];
    k = 1;

    wb = waitbar(0, 'Starting', 'name', 'Normal BFO90');

    fprintf('Processing (%d) Normal BFO90: ', numSessions);

    for iSession = 1:numSessions
        animalName = SessionTable.animalName{iSession};
        sessionName = SessionTable.sessionName{iSession};

        waitbar(iSession/numSessions, wb, sprintf('Progress: %d %% (%s %s)', floor(iSession/numSessions*100), strrep(animalName, '_', ' '), sessionName))

        R(k).animalName = animalName;
        R(k).sessionName = sessionName;

        % Now get all of the PerCell data associated with the current session
        PerCellSession = PerCell(ismember(PerCell.animalName, animalName) & ismember(PerCell.sessionName, sessionName), :);

        for iField = 1:length(fieldSets)
            fieldNamesWanted = fieldSets(iField).fieldNames;

            [hc_group, prob_group, totalSamples, numCells, numSamplesExcluded] = compute_bfo_using_percell_table(PerCellSession, BFO90_ROTATIONS_TO_CHECK, fieldNamesWanted, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION);

            R(k).(sprintf('%s_hist', fieldSets(iField).name)) = hc_group;
            R(k).(sprintf('%s_prob', fieldSets(iField).name)) = prob_group;
            R(k).(sprintf('%s_num_samples', fieldSets(iField).name)) = totalSamples;
            R(k).(sprintf('%s_num_cells', fieldSets(iField).name)) = numCells;
            
            R(k).(sprintf('%s_num_samples_excluded', fieldSets(iField).name)) = numSamplesExcluded;
            
            R(k).(sprintf('%s_excluded_prob', fieldSets(iField).name)) = numSamplesExcluded ./ (numSamplesExcluded+totalSamples) * 100;
        end
        
        R(k).filterMinMaps = BFO90_FILTER_MIN_MAPS;
        R(k).filterMinCorrelation = BFO90_FILTER_MIN_TRIAL_CORRELATION;
        R(k).rotations = BFO90_ROTATIONS_TO_CHECK;
        
        k = k + 1;
    end
    
    if ~isempty(wb)
        close(wb);
    end
    fprintf(' done!\n');

    R = struct2table(R);

end % function






