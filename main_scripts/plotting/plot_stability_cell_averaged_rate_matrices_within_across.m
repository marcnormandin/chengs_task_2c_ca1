% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig1, hFig2, hFig3, hFig4] = plot_stability_cell_averaged_rate_matrices_within_across(DaysToSessions, RateMatrices, BestAligned, ....
    BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD)
    % Get cells that are in both the RateMatrices and the BestAligned tables
    MatchTable = get_match_table_for_cells(RateMatrices, BestAligned);

    % Make a new table for cells
    RateMatricesMatch = RateMatrices(MatchTable.TableARowIndex,:);
    RateMatricesMatch.stability = BestAligned.bestCorrelation(MatchTable.TableBRowIndex) > BESTALIGNED_STABILITY_THRESHOLD_CRITERIA;
    RateMatricesMatch.bestCorrelation = BestAligned.bestCorrelation(MatchTable.TableBRowIndex);

    % Separate into stable and unstable 
    RateMatricesMatchStable = RateMatricesMatch(RateMatricesMatch.stability == true,:);
    RateMatricesMatchUnstable = RateMatricesMatch(RateMatricesMatch.stability == false,:);

    % Compute the averages
    
    RMAverageStable = compute_per_day_cell_averaged_rate_difference_matrix(DaysToSessions, RateMatricesMatchStable, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD);
    RMAverageUnstable = compute_per_day_cell_averaged_rate_difference_matrix(DaysToSessions, RateMatricesMatchUnstable, RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD);

    % Plot
    hFig1 = plot_per_day_cell_averaged_rate_difference_matrix_table(RMAverageStable);
    sgtitle(sprintf('Stable (Per cell (%s))', RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD));

    hFig2 = plot_per_day_cell_averaged_rate_difference_matrix_table(RMAverageUnstable);
    sgtitle(sprintf('Unstable (Per cell (%s))', RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD));


    % Within and Across
    RMAverageStableWithinAcross = compute_average_rate_differences_within_across_for_table(RMAverageStable);
    RMAverageUnstableWithinAcross = compute_average_rate_differences_within_across_for_table(RMAverageUnstable);

    % Plot
    hFig3 = plot_average_rate_differences_within_across_bars_table(RMAverageStableWithinAcross);
    sgtitle(sprintf('Stable (Per cell (%s))', RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD));

    hFig4 = plot_average_rate_differences_within_across_bars_table(RMAverageUnstableWithinAcross);
    sgtitle(sprintf('Unstable (Per cell (%s))', RATE_MATRICES_PER_CELL_NORMALIZATION_METHOD));
end % function

