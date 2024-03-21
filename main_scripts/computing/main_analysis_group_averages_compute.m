% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [analysisResultsGroupAverages] = main_analysis_group_averages_compute(analysisSettings, analysisResults)
    %%
    % To be made populated
    analysisResultsGroupAverages = [];
    
    % Just pass it on since it is already averaged across the days.
    analysisResultsGroupAverages.NormalBFO180SimilarityShuffled = analysisResults.NormalBFO180SimilarityShuffled;



    %% Plot BFO90 ALL Days 1,2,3
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);

    % COMPUTE
    [MeanTable, ErrorTable] = bfotable_compute_group_averages(analysisResults.NormalBFO90, SessionsToGroups);

    % STORE
    analysisResultsGroupAverages.NormalBFO90.MeanTable = MeanTable;
    analysisResultsGroupAverages.NormalBFO90.ErrorTable = ErrorTable;
    analysisResultsGroupAverages.NormalBFO90.SessionsToGroups = SessionsToGroups;


    %% Stability BFO90 (per Animal)

    % COMPUTE
    [MeanTable, ErrorTable] = bfotable_compute_group_averages(analysisResults.StabilityBFO90, SessionsToGroups);

    % STORE
    analysisResultsGroupAverages.StabilityBFO90.MeanTable = MeanTable;
    analysisResultsGroupAverages.StabilityBFO90.ErrorTable = ErrorTable;

    
    %% Stability BFO90 (per Animal) registered and classified based on a single day.
    % Can only be used with calcium since that is the only data that is
    % registered.
    % COMPUTE
    if analysisSettings.IS_CALCIUM_DATA
        numGroups = length(analysisResults.Stability90_Reg);
        for iGroup = 1:numGroups
            groupLabel = analysisResults.Stability90_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL;
            
            % Get the results to process
            StabilityBFO90_Reg = analysisResults.Stability90_Reg(iGroup).StabilityBFO90_Reg;
            
            % Compute
            [MeanTable, ErrorTable] = bfotable_compute_group_averages(StabilityBFO90_Reg, SessionsToGroups);

            % STORE
            analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).MeanTable = MeanTable;
            analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).ErrorTable = ErrorTable;
            analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL = groupLabel;
        end % iGroup
    end
    
    
    %% Stability BFO90 PER CELL
    % Do it per cell, but for individual mice first (as a diagnostic).
    animalNames = unique(analysisResults.StabilityPerCell90.animalName);
    numAnimals = length(animalNames);
    for iAnimal = 1:numAnimals
        % Use only data for a single animal
        animalName = animalNames{iAnimal};
        X = analysisResults.StabilityPerCell90;
        X(~ismember(X.animalName, animalName), :) = [];
        
        [StabilityBFO90PerCell, success] = compute_bfo_stability_table_per_cell( analysisSettings, X, analysisResults.BestAligned, SessionsToGroups);
        [StableMeanTable, StableErrorTable] = bfotable_compute_group_averages(StabilityBFO90PerCell(StabilityBFO90PerCell.isStable==true,:), SessionsToGroups);
        [UnstableMeanTable, UnstableErrorTable] = bfotable_compute_group_averages(StabilityBFO90PerCell(StabilityBFO90PerCell.isStable==false,:), SessionsToGroups);
    end

    % This is the plot using all of the animals.
    [StabilityBFO90PerCell, success] = compute_bfo_stability_table_per_cell(analysisSettings, analysisResults.StabilityPerCell90, analysisResults.BestAligned, SessionsToGroups);
    [StableMeanTable, StableErrorTable] = bfotable_compute_group_averages(StabilityBFO90PerCell(StabilityBFO90PerCell.isStable==true,:), SessionsToGroups);
    [UnstableMeanTable, UnstableErrorTable] = bfotable_compute_group_averages(StabilityBFO90PerCell(StabilityBFO90PerCell.isStable==false,:), SessionsToGroups);

    analysisResultsGroupAverages.StabilityBFO90PerCell.StableMeanTable = StableMeanTable;
    analysisResultsGroupAverages.StabilityBFO90PerCell.StableErrorTable = StableErrorTable;
    analysisResultsGroupAverages.StabilityBFO90PerCell.UnstableMeanTable = UnstableMeanTable;
    analysisResultsGroupAverages.StabilityBFO90PerCell.UnstableErrorTable = UnstableErrorTable;
    


    %% Rate Matrices
    if ~analysisSettings.IS_CALCIUM_DATA
        % COMPUTE AVERAGES
        SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
        RateMatricesAveraged = compute_per_group_animal_averaged_rate_difference_matrix(SessionsToGroups, analysisResults.RateMatricesFiltered, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL);

        % STORE
        analysisResultsGroupAverages.RateMatrices = RateMatricesAveraged;
    end

    %% Rate bars
    if ~analysisSettings.IS_CALCIUM_DATA
        % COMPUTE AND STORE
        analysisResultsGroupAverages.RateMatricesHistogram = compute_group_average_rate_differences_within_across_for_table(RateMatricesAveraged);
    end


    %% Number of stable vs unstable cells

    % COMPUTE
    BestAlignedCounts = analysisResults.BestAlignedCounts;

    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = []; % only use days 1,2,3

    uniquePairs = unique_groups_sessions_to_groups(SessionsToGroups);
    % For each group, get the animal and sessions that are members of it
    GroupCounts = [];
    k = 1;
    for iGroup = 1:size(uniquePairs,1)
       g = SessionsToGroups(SessionsToGroups.groupId == uniquePairs{iGroup,1},:);
       % For the given group, get indices into the BestAlignedCounts
       matchTable = get_match_table_by_sessions(g, BestAlignedCounts);
       
       if ~isempty(matchTable)

           gCounts = BestAlignedCounts(matchTable.TableBRowIndex,:);

           % Compute the stats
           GroupCounts(k).groupId = g.groupId(1); % all will be the same
           GroupCounts(k).groupLabel = g.groupLabel{1}; % all will be the same
           GroupCounts(k).meanPercentStable = mean(gCounts.percentStable);
           GroupCounts(k).errorPercentStable = std(gCounts.percentStable) ./ sqrt(size(gCounts,1));
           k = k + 1;
       end
    end
    GroupCounts = struct2table(GroupCounts);

    % STORE
    analysisResultsGroupAverages.StableCells = GroupCounts;



    %% STABILITY rate matrices and histograms for days 1,2,3
    if ~analysisSettings.IS_CALCIUM_DATA
        % Groups
        SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
        SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];

        % Compute
        [RMAverageStable, RMAverageUnstable, RMAverageStableWithinAcross, RMAverageUnstableWithinAcross] = compute_stability_animal_averaged_rate_matrices_within_across(SessionsToGroups, analysisResults.RateMatricesFiltered, analysisResults.BestAligned, ...
            analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL);

        % STORE
        analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStable = RMAverageStable;
        analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstable = RMAverageUnstable;
        analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStableWithinAcross = RMAverageStableWithinAcross;
        analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstableWithinAcross = RMAverageUnstableWithinAcross;
    end

    
    %% SAVE THE DATA
    % All the bfo180 similarity since we need all of it for plotting
    analysisResultsGroupAverages.NormalBFO180Similarity = analysisResults.NormalBFO180Similarity;
    analysisResultsGroupAverages.StabilityBFO180Similarity = analysisResults.StabilityBFO180Similarity;
    
    if analysisSettings.IS_CALCIUM_DATA
        analysisResultsGroupAverages.Stability180_Reg = analysisResults.Stability180_Reg;
    end
    
    %% Compute the kstests
    [R] = compute_stats_stability_bfo180_sim_curves_days123(analysisSettings, analysisResultsGroupAverages);
    save(fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'kstests_stability_bfo180_sim_curves_days123.mat'), 'analysisSettings', 'R');
end % function

