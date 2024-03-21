% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R, success] = compute_normal_bfo180_similarity_table(PerCell180, BFO180_FILTER_MIN_MAPS, BFO180_FILTER_MIN_TRIAL_CORRELATION)

    PerCell180.animalSessionCellId = nan(size(PerCell180,1),1);
    %PerCell180.animalAnimalSessionCellName = cell(size(PerCell180,1),1);

    % Add a unique id so that we can do paired stats tests
    uniqueAnimalSessionCellNames = unique(PerCell180.animalSessionCellName);
    for iRow = 1:size(PerCell180,1)
       an = PerCell180.animalSessionCellName{iRow};
       anid = find(ismember(uniqueAnimalSessionCellNames, an));
       PerCell180.animalSessionCellId(iRow) = anid;
    end

    fieldSets = helper_load_bfo180_fieldsets_to_record();

    % Get the unique sessions for which we will compute one BF90 each
    SessionTable = get_sessions_from_table(PerCell180);
    numSessions = size(SessionTable,1);

    R = [];
    success = true;
    wb = [];
    k = 1;

    wb = waitbar(0, 'Starting', 'name', 'Normal BF180 Similarity');

    fprintf('Processing (%d) Normal BF180 Similarity: ', numSessions);

    for iSession = 1:numSessions
        animalName = SessionTable.animalName{iSession};
        sessionName = SessionTable.sessionName{iSession};

        waitbar(iSession/numSessions, wb, sprintf('Progress: %d %% (%s %s)', floor(iSession/numSessions*100), strrep(animalName, '_', ' '), sessionName))

        R(k).animalName = animalName;
        R(k).sessionName = sessionName;

        % Now get all of the PerCell180 data associated with the current session
        PerCell180Session = PerCell180(ismember(PerCell180.animalName, animalName) & ismember(PerCell180.sessionName, sessionName), :);

        for iField = 1:length(fieldSets)
            fieldNamesWanted = fieldSets(iField).fieldNames;

            [averageCorrelations, totalSamples, numCells, numSamplesExcluded, animalSessionCellIds, animalSessionCellNames] = compute_bfo180_similarity_using_percell_table(PerCell180Session, fieldNamesWanted, BFO180_FILTER_MIN_MAPS);

            R(k).(sprintf('%s_averageCorrelations', fieldSets(iField).name)) = averageCorrelations;

            R(k).(sprintf('%s_num_samples', fieldSets(iField).name)) = totalSamples;
            R(k).(sprintf('%s_num_cells', fieldSets(iField).name)) = numCells;
            R(k).(sprintf('%s_num_samples_excluded', fieldSets(iField).name)) = numSamplesExcluded;
            
            R(k).(sprintf('%s_animal_session_cell_ids', fieldSets(iField).name)) = animalSessionCellIds;
            
            R(k).(sprintf('%s_animal_session_cell_names', fieldSets(iField).name)) = animalSessionCellNames;
        end
        
        R(k).filterMinMaps = BFO180_FILTER_MIN_MAPS;
        R(k).filterMinCorrelation = BFO180_FILTER_MIN_TRIAL_CORRELATION;
        R(k).rotations = [0, 180];
        
        k = k + 1;
    end
    
    if ~isempty(wb)
        close(wb);
    end
    fprintf(' done!\n');

    R = struct2table(R);

end % function






