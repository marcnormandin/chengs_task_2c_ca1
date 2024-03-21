% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [analysisResults] = main_analysis_compute_tables(analysisSettings, analysisInput)

    isCalciumData = analysisSettings.IS_CALCIUM_DATA;
    
    MapsData = analysisInput.MapsData;
    MapsSquareData = analysisInput.MapsSquareData;
        
    

    %% Compute the PerCell table used for the NORMAL BFO 90
    % There will be one row per unique cell
    [NormalPerCell90, success] = compute_percell_for_table(MapsSquareData, analysisSettings.NORMAL_BFO90_ROTATIONS_DEG, analysisSettings.NORMAL_BFO90_CONTEXTS_TO_REFLECT_MAPS, analysisSettings.NORMAL_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    if ~success
        fprintf('Error while processing PerCell correlations\n');
    end
    
    
    
 

    %% Compute the PerCell table used for cumulative similarity using 0, 180
    % There will be one row per unique cell
    [NormalPerCell180, success] = compute_percell_for_table(MapsData, analysisSettings.NORMAL_BFO180_ROTATIONS_DEG, analysisSettings.NORMAL_BFO180_CONTEXTS_TO_REFLECT_MAPS, analysisSettings.NORMAL_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    if ~success
        fprintf('Error while processing PerCell correlations\n');
    end



    %% Compute the BFO180 similarity (average best correlations)
    [NormalBFO180Similarity, success] = compute_normal_bfo180_similarity_table(NormalPerCell180, analysisSettings.NORMAL_BFO180_FILTER_MIN_MAPS, analysisSettings.NORMAL_BFO180_FILTER_MIN_TRIAL_CORRELATION);
    if ~success
        fprintf('Error while processing Normal BFO180 Similarity\n');
    end



    %% Compute the best aligned maps and their associated correlation
    % There will be one row per unique cell
    if strcmpi(analysisSettings.BESTALIGNED_VERSION, 'v2')
        [BestAligned, success] = ml_alg_compute_best_aligned_map_table_v2(MapsData, analysisSettings.BESTALIGNED_COMPARISON_METHOD, false);
    elseif strcmpi(analysisSettings.BESTALIGNED_VERSION, 'v3')
        [BestAligned, success] = ml_alg_compute_best_aligned_map_table_v3(MapsData, analysisSettings.BESTALIGNED_COMPARISON_METHOD, false);
    else
        error('Unsupported best aligned algorithm version. Must be v2 or v3.');
    end
    if ~success
        fprintf('Error while processing Best Aligned Maps\n');
    end
    
    
    %% Compute the best aligned counts for stable versus unstable cells
    BestAligned.isStable = BestAligned.bestCorrelation >= analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA;

    % Counts
    BestAlignedCounts = get_sessions_from_table(BestAligned);
    numStable = zeros(size(BestAlignedCounts,1),1);
    numUnstable = zeros(size(BestAlignedCounts,1),1);
    total = zeros(size(BestAlignedCounts,1),1);
    percentStable = zeros(size(BestAlignedCounts,1),1);

    for i = 1:size(BestAlignedCounts,1)
        sessionData = BestAligned(ismember(BestAligned.animalName, BestAlignedCounts.animalName{i}) & ismember(BestAligned.sessionName, BestAlignedCounts.sessionName{i}),:);
        numStable(i) = sum(sessionData.isStable == true);
        numUnstable(i) = sum(sessionData.isStable == false);
        total(i) = numStable(i) + numUnstable(i);
        percentStable(i) = numStable(i) / total(i) * 100;
    end

    BestAlignedCounts.numStable = numStable;
    BestAlignedCounts.numUnstable = numUnstable;
    BestAlignedCounts.total = total;
    BestAlignedCounts.percentStable = percentStable;



    %% Compute all of the rate matrices
    % There will be one row per unique cell
    if ~isCalciumData
        [RateMatrices, success] = compute_rate_matrix_table(MapsData, analysisSettings.RATE_MATRIX_RATETYPE, analysisSettings.RATE_MATRIX_NUM_TRIALS);
        if ~success
            fprintf('Error while processing Rate Matrices (%s)\n', analysisSettings.RATE_MATRIX_RATETYPE);
        end
    end



    %% Compute the BFO90 using BestAligned being stable or unstable. This allows for one context to be reflected
    %  if the option is set in the settings. This can be used with tetrode
    %  or calcium data.
    animalSessionCellNamesStable = BestAligned.animalSessionCellName(BestAligned.isStable == true);
    animalSessionCellNamesUnstable = BestAligned.animalSessionCellName(BestAligned.isStable == false);
    [PerCell90Stable, success] = compute_percell_for_table(MapsSquareData(ismember(MapsSquareData.animalSessionCellName, animalSessionCellNamesStable),:), analysisSettings.STABILITY_BFO90_ROTATIONS_DEG, [], analysisSettings.STABILITY_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    PerCell90Stable.isStable = true(size(PerCell90Stable,1),1);
    if ~success
        fprintf('Error while processing PerCell correlations for stable\n');
    end
    [PerCell90Unstable, success] = compute_percell_for_table(MapsSquareData(ismember(MapsSquareData.animalSessionCellName, animalSessionCellNamesUnstable),:), analysisSettings.STABILITY_BFO90_ROTATIONS_DEG, analysisSettings.STABILITY_BFO90_REFLECT_UNSTABLE_CONTEXT, analysisSettings.STABILITY_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    if ~success
        fprintf('Error while processing PerCell correlations for unstable\n');
    end
    PerCell90Unstable.isStable = false(size(PerCell90Unstable,1),1);
    StabilityPerCell90 = cat(1, PerCell90Stable, PerCell90Unstable);
    
    % Compute the PER animal stability data. Note that the returned table
    % will have one BFO90 result PER SESSION of each animal. The group
    % computation code, later on, will combine sessions of animals to compute each
    % day's per animal average BFO90.
    [StabilityBFO90, success] = compute_stability_bfo90_table_v2(StabilityPerCell90, analysisSettings.STABILITY_BFO90_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO90_FILTER_MIN_TRIAL_CORRELATION);
    if ~success
        fprintf('Error while processing Stability BF90\n');
    else
        % Create the xls file for Isabel
        X = helper_add_group_labels_to_table(StabilityBFO90, load_sessions_to_groups_table(analysisSettings), true);
        X(~ismember(X.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
        if analysisSettings.IS_CALCIUM_DATA
            fn = fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'calcium_stability_bfo90_animal_withinacross_days123.xlsx');
        else
            fn = fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'tetrodes_stability_bfo90_animal_withinacross_days123.xlsx');
        end 
        writetable(X, fn)
    end
    
    %% STABILITY BF180 to be used for the cumulative distributions
    [StabilityPerCell180, StabilityBFO180Similarity] =  helper_compute_stability_bfo180(analysisSettings, analysisInput.MapsData, BestAligned);
    
    % Registered cells that are classified on a single day.
    if isCalciumData
        groupLabels = {'Day 1', 'Day 2', 'Day 3'};
        for iGroup = 1:length(groupLabels)
            groupLabel = groupLabels{iGroup};
            [StabilityPerCell180_Reg, StabilityBFO180Similarity_Reg] = helper_compute_stability_reg_1day_bfo180(analysisSettings, analysisInput.CellRegTable, analysisInput.MapsData, BestAligned, groupLabel);
            
            analysisResults.Stability180_Reg(iGroup).StabilityPerCell180_Reg = StabilityPerCell180_Reg;
            analysisResults.Stability180_Reg(iGroup).StabilityBFO180Similarity_Reg = StabilityBFO180Similarity_Reg;
            analysisResults.Stability180_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL = groupLabel;
            
            % No need to save since I give Celia the mat files.
        end % iGroup
    end

   
    %% Compute stability per cell and bfo90 for calcium data, and registered cells, where the classification is based on a given day.
    if isCalciumData
        groupLabels = {'Day 1', 'Day 2', 'Day 3'};
        for iGroup = 1:length(groupLabels)
            groupLabel = groupLabels{iGroup};
            [StabilityPerCell90_Reg, StabilityBFO90_Reg] = compute_stability_reg_1day_BFO90_per_animal(analysisSettings, analysisInput, BestAligned, groupLabel);

            % Create the xls file for Isabel
            X = helper_add_group_labels_to_table(StabilityBFO90_Reg, load_sessions_to_groups_table(analysisSettings), true);
            X(~ismember(X.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
            writetable(X, fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, sprintf('calcium_stability_bfo90_reg_%s_animal_withinacross_days123.xlsx', strrep(groupLabel, ' ', '_'))));
            
            % Now store it
            analysisResults.Stability90_Reg(iGroup).StabilityPerCell90_Reg = StabilityPerCell90_Reg;
            analysisResults.Stability90_Reg(iGroup).StabilityBFO90_Reg = StabilityBFO90_Reg;
            analysisResults.Stability90_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL = groupLabel; 
        end
    end


    %% Compute the Normal BFO90
    [NormalBFO90, success] = compute_normal_bfo90_table(NormalPerCell90, analysisSettings.NORMAL_BFO90_FILTER_MIN_MAPS, analysisSettings.NORMAL_BFO90_FILTER_MIN_TRIAL_CORRELATION);
    if ~success
        fprintf('Error while processing Stability BF90\n');
    end

    %% Apply the filtering to tables if applicable
    if ~isCalciumData
        RateMatricesFiltered = RateMatrices(RateMatrices.numRates >= analysisSettings.RATE_MATRIX_FILTER_MIN_TRIALS,:);
    end
    
    
    %% Compute the shuffles
    NormalBFO180SimilarityShuffled = compute_bfo180_similarity_curves_shuffled(analysisSettings, analysisInput);

    
   
    
    %% STORE
    analysisResults.NormalPerCell90 = NormalPerCell90;
    analysisResults.NormalPerCell180 = NormalPerCell180;
    
    analysisResults.NormalBFO90 = NormalBFO90;
    analysisResults.NormalBFO180Similarity = NormalBFO180Similarity;
    analysisResults.NormalBFO180SimilarityShuffled = NormalBFO180SimilarityShuffled;
    
    analysisResults.BestAligned = BestAligned;
    analysisResults.BestAlignedCounts = BestAlignedCounts;
    
    % One of the many versions of stability because it keeps being asked to
    % be modified.
    analysisResults.StabilityPerCell90 = StabilityPerCell90;
    analysisResults.StabilityBFO90 = StabilityBFO90;
    analysisResults.StabilityPerCell180 = StabilityPerCell180;
    analysisResults.StabilityBFO180Similarity = StabilityBFO180Similarity;
    
    if isCalciumData
        %analysisResults.StabilityPerCell90_Reg = StabilityPerCell90_Reg;
        %analysisResults.StabilityBFO90_Reg = StabilityBFO90_Reg;
        
        %analysisResults.StabilityPerCell180_Reg = StabilityPerCell180_Reg;
        %analysisResults.StabilityBFO180Similarity_Reg = StabilityBFO180Similarity_Reg;
    end
        
    if ~isCalciumData
        analysisResults.RateMatrices = RateMatrices;
        analysisResults.RateMatricesFiltered = RateMatricesFiltered;
    end

end % function




    
function [StabilityPerCell180, StabilityBFO180Similarity] =  helper_compute_stability_bfo180(analysisSettings, MapsData, BestAligned)
% Compute the BFO180 using BestAligned being stable or unstable. This allows for one context to be reflected
%  if the option is set in the settings. This can be used with tetrode
%  or calcium data.
animalSessionCellNamesStable = BestAligned.animalSessionCellName(BestAligned.isStable == true);
animalSessionCellNamesUnstable = BestAligned.animalSessionCellName(BestAligned.isStable == false);

[PerCell180Stable, success] = compute_percell_for_table(MapsData(ismember(MapsData.animalSessionCellName, animalSessionCellNamesStable),:), analysisSettings.STABILITY_BFO180_ROTATIONS_DEG, [], analysisSettings.STABILITY_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
PerCell180Stable.isStable = true(size(PerCell180Stable,1),1);
if ~success
    fprintf('Error while processing PerCell correlations for stable\n');
end
[PerCell180Unstable, success] = compute_percell_for_table(MapsData(ismember(MapsData.animalSessionCellName, animalSessionCellNamesUnstable),:), analysisSettings.STABILITY_BFO180_ROTATIONS_DEG, analysisSettings.STABILITY_BFO180_REFLECT_UNSTABLE_CONTEXT, analysisSettings.STABILITY_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
if ~success
    fprintf('Error while processing PerCell correlations for unstable\n');
end
PerCell180Unstable.isStable = false(size(PerCell180Unstable,1),1);
StabilityPerCell180 = cat(1, PerCell180Stable, PerCell180Unstable);
    

% Compute the BFO180 similarity (average best correlations)
[StabilityBFO180Similarity_stable, success] = compute_normal_bfo180_similarity_table(StabilityPerCell180(StabilityPerCell180.isStable == true,:), analysisSettings.STABILITY_BFO180_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO180_FILTER_MIN_TRIAL_CORRELATION);
if ~success
    fprintf('Error while processing Stability BFO180 (Stable) Similarity\n');
end
[StabilityBFO180Similarity_unstable, success] = compute_normal_bfo180_similarity_table(StabilityPerCell180(StabilityPerCell180.isStable == false,:), analysisSettings.STABILITY_BFO180_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO180_FILTER_MIN_TRIAL_CORRELATION);
if ~success
    fprintf('Error while processing Stability BFO180 (Stable) Similarity\n');
end
StabilityBFO180Similarity.stable = StabilityBFO180Similarity_stable;
StabilityBFO180Similarity.unstable = StabilityBFO180Similarity_unstable;
end % function


%% This only applies to Calcium data since we need to have registered cells (across days)
function [StabilityPerCell180, StabilityBFO180Similarity] =  helper_compute_stability_reg_1day_bfo180(analysisSettings, CellRegTable, MapsData, BestAligned, STABILITY_CLASSIFICATION_GROUP_LABEL)
if ~analysisSettings.IS_CALCIUM_DATA
    error('This function is only applicable for calcium data since it requires registered cells.\n');
end

% % Compute the BFO180 using BestAligned being stable or unstable. This allows for one context to be reflected
% %  if the option is set in the settings. This can be used with tetrode
% %  or calcium data.
% animalSessionCellNamesStable = BestAligned.animalSessionCellName(BestAligned.isStable == true);
% animalSessionCellNamesUnstable = BestAligned.animalSessionCellName(BestAligned.isStable == false);

SessionsToGroups = load_sessions_to_groups_table(analysisSettings);

% Classify stablity the cells based on a given day (from the settings)
StabilityClassifications = classify_stability_registered_cells_on_specified_day(BestAligned, CellRegTable, SessionsToGroups, STABILITY_CLASSIFICATION_GROUP_LABEL, analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);

% Extract the stable and unstable registered names
registeredCellNamesStable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == true);
registeredCellNamesUnstable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == false);

    

[PerCell180Stable, success] = compute_percell_for_table(MapsData(ismember(MapsData.registeredCellName, registeredCellNamesStable),:), analysisSettings.STABILITY_BFO180_ROTATIONS_DEG, [], analysisSettings.STABILITY_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
PerCell180Stable.isStable = true(size(PerCell180Stable,1),1);
if ~success
    fprintf('Error while processing PerCell correlations for stable\n');
end
[PerCell180Unstable, success] = compute_percell_for_table(MapsData(ismember(MapsData.registeredCellName, registeredCellNamesUnstable),:), analysisSettings.STABILITY_BFO180_ROTATIONS_DEG, analysisSettings.STABILITY_BFO180_REFLECT_UNSTABLE_CONTEXT, analysisSettings.STABILITY_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
if ~success
    fprintf('Error while processing PerCell correlations for unstable\n');
end
PerCell180Unstable.isStable = false(size(PerCell180Unstable,1),1);
StabilityPerCell180 = cat(1, PerCell180Stable, PerCell180Unstable);
    

% Compute the BFO180 similarity (average best correlations)
[StabilityBFO180Similarity_stable, success] = compute_normal_bfo180_similarity_table(StabilityPerCell180(StabilityPerCell180.isStable == true,:), analysisSettings.STABILITY_BFO180_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO180_FILTER_MIN_TRIAL_CORRELATION);
if ~success
    fprintf('Error while processing Stability BFO180 (Stable) Similarity\n');
end
[StabilityBFO180Similarity_unstable, success] = compute_normal_bfo180_similarity_table(StabilityPerCell180(StabilityPerCell180.isStable == false,:), analysisSettings.STABILITY_BFO180_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO180_FILTER_MIN_TRIAL_CORRELATION);
if ~success
    fprintf('Error while processing Stability BFO180 (Stable) Similarity\n');
end
StabilityBFO180Similarity.stable = StabilityBFO180Similarity_stable;
StabilityBFO180Similarity.unstable = StabilityBFO180Similarity_unstable;
end % function

