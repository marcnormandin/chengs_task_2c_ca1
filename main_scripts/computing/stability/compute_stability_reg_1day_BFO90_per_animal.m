% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Compute the BFO90 per (animal, session) pairs later to be used to average
% over animal per day.
function [StabilityPerCell90, StabilityBFO90] = compute_stability_reg_1day_BFO90_per_animal(analysisSettings, analysisInput, BestAligned, STABILITY_CLASSIFICATION_GROUP_LABEL)
% Only to be used with calcium data since we need to have cells registered
% across days/sessions.

    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);

    % Classify stablity the cells based on a given day (from the settings)
    StabilityClassifications = classify_stability_registered_cells_on_specified_day(BestAligned, analysisInput.CellRegTable, SessionsToGroups, STABILITY_CLASSIFICATION_GROUP_LABEL, analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);

    % Extract the stable and unstable registered names
    registeredCellNamesStable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == true);
    registeredCellNamesUnstable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == false);

    % Simply notation
    MapsSquareData = analysisInput.MapsSquareData;

    % Compute for only stable per (animal, session)
    [PerCell90Stable, success] = compute_percell_for_table(MapsSquareData(ismember(MapsSquareData.registeredCellName, registeredCellNamesStable),:), analysisSettings.STABILITY_BFO90_ROTATIONS_DEG, [], analysisSettings.STABILITY_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    PerCell90Stable.isStable = true(size(PerCell90Stable,1),1);
    if ~success
        fprintf('Error while processing PerCell correlations for stable\n');
    end
    % Compute for only unstable per (animal, session)
    [PerCell90Unstable, success] = compute_percell_for_table(MapsSquareData(ismember(MapsSquareData.registeredCellName, registeredCellNamesUnstable),:), analysisSettings.STABILITY_BFO90_ROTATIONS_DEG, analysisSettings.STABILITY_BFO90_REFLECT_UNSTABLE_CONTEXT, analysisSettings.STABILITY_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    if ~success
        fprintf('Error while processing PerCell correlations for unstable\n');
    end
    PerCell90Unstable.isStable = false(size(PerCell90Unstable,1),1);

    % Combine the results
    StabilityPerCell90 = cat(1, PerCell90Stable, PerCell90Unstable);

    % Compute the PER animal stability data. Note that the returned table
    % will have one BFO90 result PER SESSION of each animal. The group
    % computation code, later on, will combine sessions of animals to compute each
    % day's per animal average BFO90.
    [StabilityBFO90, success] = compute_stability_bfo90_table_v2(StabilityPerCell90, analysisSettings.STABILITY_BFO90_FILTER_MIN_MAPS, analysisSettings.STABILITY_BFO90_FILTER_MIN_TRIAL_CORRELATION);
    if ~success
        fprintf('Error while processing Stability BF90\n');
    end

end % function
    