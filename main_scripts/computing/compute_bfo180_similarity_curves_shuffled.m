% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_bfo180_similarity_curves_shuffled(analysisSettings, analysisInput)

    % If no shuffles are wanted, return an empty table
    if isempty(analysisSettings.NORMAL_BFO180_SIMILARITY_NSHUFFLES) || analysisSettings.NORMAL_BFO180_SIMILARITY_NSHUFFLES <= 0
        R = [];
        return
    end

    DaysToSessions = load_sessions_to_groups_table(analysisSettings);

    nShuffles = analysisSettings.NORMAL_BFO180_SIMILARITY_NSHUFFLES;

    MapsData = analysisInput.MapsData;
    % Add the day labels
    MapsData = helper_add_group_labels_to_table(MapsData, DaysToSessions, true);
    
    dayLabels = unique(MapsData.groupLabel);
    nDayLabels = length(dayLabels);

    
    nMaps = size(MapsData,1);

    shuffledRows = nan(nShuffles, nMaps);
    for iShuffle = 1:nShuffles
        shuffledRows(iShuffle,:) = randperm(nMaps, nMaps);
    end
    shuffledRows = reshape(shuffledRows, nMaps, nShuffles);

    TShuffles = [];

    for iShuffle = 1:nShuffles
        fprintf('Computing shuffle distribution: %d / %d\n', iShuffle, nShuffles);

        ShuffledMapsData = MapsData;
        
        % This will shuffles maps across days, but then the curves are the
        % same per day.
        %ShuffledMapsData.map = ShuffledMapsData.map(shuffledRows(:,iShuffle));
        
        % Shuffle maps only within the same day.
        for iDay = 1:nDayLabels
            dayLabel = dayLabels{iDay};
            dayInds = find(ismember(MapsData.groupLabel, dayLabel));
            randInds = randperm(length(dayInds), length(dayInds));
            
            ShuffledMapsData.map(reshape(dayInds, length(dayInds), 1)) = ShuffledMapsData.map(reshape(dayInds(randInds), length(dayInds),1));
        end
        

        % Compute the PerCell table used for cumulative similarity using 0, 180
        % There will be one row per unique cell
        [NormalPerCell180, success] = compute_percell_for_table(ShuffledMapsData, analysisSettings.NORMAL_BFO180_ROTATIONS_DEG, analysisSettings.NORMAL_BFO180_CONTEXTS_TO_REFLECT_MAPS, analysisSettings.NORMAL_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
        if ~success
            fprintf('Error while processing PerCell correlations\n');
        end

        % Compute the BFO180 similarity (average best correlations)
        [NormalBFO180Similarity, success] = compute_normal_bfo180_similarity_table(NormalPerCell180, analysisSettings.NORMAL_BFO180_FILTER_MIN_MAPS, analysisSettings.NORMAL_BFO180_FILTER_MIN_TRIAL_CORRELATION);
        if ~success
            fprintf('Error while processing Normal BFO180 Similarity\n');
        end

        [T] = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, NormalBFO180Similarity);
        T.shuffleId = iShuffle * ones(size(T,1),1);

        if isempty(TShuffles)
            TShuffles = T;
        else
            TShuffles = cat(1, TShuffles, T);
        end

    end % iShuffle


    dayLabels = unique(TShuffles.dayLabel);
    nDayLabels = length(dayLabels);
    fieldNames = {'all_averageCorrelations', 'all_averageCorrelations', 'context2_averageCorrelations', 'within_averageCorrelations', 'across_averageCorrelations'};
    nFieldNames = length(fieldNames);
    R = [];
    k = 1;
    for iDay = 1:nDayLabels
        dayLabel = dayLabels{iDay};
        for iField = 1:nFieldNames
            fieldName = fieldNames{iField};

            % All shuffles asssociated with the given day
            S = TShuffles(ismember(TShuffles.dayLabel, dayLabel),:);

            nRows = size(S,1);
            x = [];
            for iRow = 1:nRows
                x = cat(1, x, S.(fieldName){iRow});
            end

            R(k).dayLabel = dayLabel;
            R(k).(sprintf('%s_samples', fieldName)) = x;


            [uzt,czt] = ml_alg_cumdist(x);
            czt = reshape(czt, length(czt), 1);
            R(k).(sprintf('%s_x', fieldName)) = uzt;
            R(k).(sprintf('%s_y', fieldName)) = czt;
        end % iField
        k = k + 1;
    end % iDay
    R = struct2table(R);

end % function
