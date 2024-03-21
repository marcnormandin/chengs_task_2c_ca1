% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [RMAverageStable, RMAverageUnstable, RMAverageStableWithinAcross, RMAverageUnstableWithinAcross] = ...
    compute_stability_animal_averaged_rate_matrices_within_across(SessionsToGroups, RateMatrices, BestAligned, ...
    BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD, RATE_MATRICES_PER_ANIMAL_NORMALIZATION_METHOD)


    % Get cells that are in both the RateMatrices and the BestAligned tables
    MatchTable = get_match_table_for_cells(RateMatrices, BestAligned);

    % Make a new table for cells
    RateMatricesMatch = RateMatrices(MatchTable.TableARowIndex,:);
    RateMatricesMatch.stability = BestAligned.bestCorrelation(MatchTable.TableBRowIndex) >= BESTALIGNED_STABILITY_THRESHOLD_CRITERIA;
    RateMatricesMatch.bestCorrelation = BestAligned.bestCorrelation(MatchTable.TableBRowIndex);

    % Separate into stable and unstable 
    RateMatricesMatchStable = RateMatricesMatch(RateMatricesMatch.stability == true,:);
    RateMatricesMatchUnstable = RateMatricesMatch(RateMatricesMatch.stability == false,:);

    % Compute the averages
    % previous was: compute_per_day_animal_averaged_rate_difference_matrix
    RMAverageStable = compute_per_group_animal_averaged_rate_difference_matrix(SessionsToGroups, RateMatricesMatchStable, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD, RATE_MATRICES_PER_ANIMAL_NORMALIZATION_METHOD);
    RMAverageUnstable = compute_per_group_animal_averaged_rate_difference_matrix(SessionsToGroups, RateMatricesMatchUnstable, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD, RATE_MATRICES_PER_ANIMAL_NORMALIZATION_METHOD);


    % Within and Across
    RMAverageStableWithinAcross = [];
    if ~isempty(RMAverageStable)
        RMAverageStableWithinAcross = compute_group_average_rate_differences_within_across_for_table(RMAverageStable);
    end
    
    RMAverageUnstableWithinAcross = [];
    if ~isempty(RMAverageUnstable)
        RMAverageUnstableWithinAcross = compute_group_average_rate_differences_within_across_for_table(RMAverageUnstable);
    end
