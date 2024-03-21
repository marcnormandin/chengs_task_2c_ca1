% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [registeredCellNamesStable, registeredCellNamesUnstable] =  get_stable_unstable_registered_cell_names_classified_1day(analysisSettings, CellRegTable, BestAligned, STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA)
    % This function returns a set of stable registered cell names and a set
    % of unstable registered cell names as classified on one day (given).
    
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);

    % Classify stablity the cells based on a given day (from the settings)
    StabilityClassifications = classify_stability_registered_cells_on_specified_day(BestAligned, CellRegTable, SessionsToGroups, STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);

    % Extract the stable and unstable registered names
    registeredCellNamesStable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == true);
    registeredCellNamesUnstable = StabilityClassifications.registeredCellName(StabilityClassifications.isStable == false);
end % function
