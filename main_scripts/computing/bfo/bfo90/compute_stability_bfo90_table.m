% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R, success] = compute_stability_bfo90_table(PerCell, BestAligned, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION)

    BFO90_ROTATIONS_TO_CHECK = [0, 90, 180, 270];


    fieldSets = helper_load_bfo90_fieldsets_to_record();

    % We want matches of animalName, sessionName, and cellName between the
    % tables because they might not be the same.
    MatchTable = get_match_table_for_cells(PerCell, BestAligned);
    

    % Get the matching tables
    %cellTable = get_sessions_and_cellnames_from_table(PerCell)
    PerCellMatched = PerCell(MatchTable.TableARowIndex,:);
    BestAlignedMatched = BestAligned(MatchTable.TableBRowIndex,:);

    % Categorize stability. Add a column to PerCell.
    PerCellMatched.stable = BestAlignedMatched.bestCorrelation(:) >= BESTALIGNED_STABILITY_THRESHOLD_CRITERIA;
    
%     PerCellMatched = PerCell;
%     PerCellMatched.stable = PerCellMatched.isStable;

    % Get the unique sessions for which we will compute one BF90 each
    SessionTable = get_sessions_from_table(PerCellMatched);
    numSessions = size(SessionTable,1);

    R = [];
    success = true;
    wb = [];
    k = 1;

    wb = waitbar(0, 'Starting', 'name', 'Stable/Unstable BFO90');

    fprintf('Processing (%d) Stable/Unstable BFO90: ', numSessions);

    for iSession = 1:numSessions
        animalName = SessionTable.animalName{iSession};
        sessionName = SessionTable.sessionName{iSession};



        waitbar(iSession/numSessions, wb, sprintf('Progress: %d %% (%s %s)', floor(iSession/numSessions*100), strrep(animalName, '_', ' '), sessionName))

        R(k).animalName = animalName;
        R(k).sessionName = sessionName;

        % Now get all of the PerCell data associated with the current session
        PerCellMatchedSession = PerCellMatched(ismember(PerCellMatched.animalName, animalName) & ismember(PerCellMatched.sessionName, sessionName), :);

        % Unstable
        for iField = 1:length(fieldSets)
            fieldNamesWanted = fieldSets(iField).fieldNames;

            [hc_group, prob_group, totalSamples, numCells, numSamplesExcluded] = compute_bfo_using_percell_table(PerCellMatchedSession(PerCellMatchedSession.stable==0,:), BFO90_ROTATIONS_TO_CHECK, fieldNamesWanted, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION);

            R(k).(sprintf('%s_%s_hist', 'unstable', fieldSets(iField).name)) = hc_group;
            R(k).(sprintf('%s_%s_prob', 'unstable', fieldSets(iField).name)) = prob_group;
            R(k).(sprintf('%s_%s_num_samples', 'unstable', fieldSets(iField).name)) = totalSamples;
            R(k).(sprintf('%s_%s_num_cells', 'unstable', fieldSets(iField).name)) = numCells;
            
            R(k).(sprintf('%s_%s_num_samples_excluded', 'unstable', fieldSets(iField).name)) = numSamplesExcluded;
            R(k).(sprintf('%s_%s_excluded_prob', 'unstable', fieldSets(iField).name)) = numSamplesExcluded ./ (numSamplesExcluded+totalSamples) * 100;
        end

        % Stable
        for iField = 1:length(fieldSets)
            fieldNamesWanted = fieldSets(iField).fieldNames;

            [hc_group, prob_group, totalSamples, numCells, numSamplesExcluded] = compute_bfo_using_percell_table(PerCellMatchedSession(PerCellMatchedSession.stable==1,:), BFO90_ROTATIONS_TO_CHECK, fieldNamesWanted, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION);

            R(k).(sprintf('%s_%s_hist', 'stable', fieldSets(iField).name)) = hc_group;
            R(k).(sprintf('%s_%s_prob', 'stable', fieldSets(iField).name)) = prob_group;
            R(k).(sprintf('%s_%s_num_samples', 'stable', fieldSets(iField).name)) = totalSamples;
            R(k).(sprintf('%s_%s_num_cells', 'stable', fieldSets(iField).name)) = numCells;
            
            R(k).(sprintf('%s_%s_num_samples_excluded', 'stable', fieldSets(iField).name)) = numSamplesExcluded;
            R(k).(sprintf('%s_%s_excluded_prob', 'stable', fieldSets(iField).name)) = numSamplesExcluded ./ (numSamplesExcluded+totalSamples) * 100;
        end

        % Stable or Unstable
        for iField = 1:length(fieldSets)
            fieldNamesWanted = fieldSets(iField).fieldNames;

            [hc_group, prob_group, totalSamples, numCells, numSamplesExcluded] = compute_bfo_using_percell_table(PerCellMatchedSession, BFO90_ROTATIONS_TO_CHECK, fieldNamesWanted, BFO90_FILTER_MIN_MAPS, BFO90_FILTER_MIN_TRIAL_CORRELATION);

            R(k).(sprintf('%s_%s_hist', 'either', fieldSets(iField).name)) = hc_group;
            R(k).(sprintf('%s_%s_prob', 'either', fieldSets(iField).name)) = prob_group;
            R(k).(sprintf('%s_%s_num_samples', 'either', fieldSets(iField).name)) = totalSamples;
            R(k).(sprintf('%s_%s_num_cells', 'either', fieldSets(iField).name)) = numCells;
            
            R(k).(sprintf('%s_%s_num_samples_excluded', 'either', fieldSets(iField).name)) = numSamplesExcluded;
            R(k).(sprintf('%s_%s_excluded_prob', 'either', fieldSets(iField).name)) = numSamplesExcluded ./ (numSamplesExcluded+totalSamples) * 100;
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






