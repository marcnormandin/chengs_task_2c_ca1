% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [analysisSettings] = analysis_settings_load()
    analysisSettings = [];
    
    % This should never be set manually here, as the shuffled pipeline will
    % set them.
    analysisSettings.USE_SHUFFLED_MAPS_AT_INITIALIZATION = false;
    analysisSettings.NUM_SHUFFLES = 0;
    
    
    % TYPE OF DATA
    analysisSettings.IS_CALCIUM_DATA = false; % default. Calcium run script overrides this.
    
    % Global value
    BFO90_FILTER_MIN_MAPS = 2;

    % NORMAL BFO 90 settings
    analysisSettings.NORMAL_BFO90_ROTATIONS_DEG = [0, 90, 180, 270];
    analysisSettings.NORMAL_BFO90_FILTER_MIN_MAPS = BFO90_FILTER_MIN_MAPS;
    analysisSettings.NORMAL_BFO90_FILTER_MIN_TRIAL_CORRELATION = 0.3;
    analysisSettings.NORMAL_BFO90_CONTEXTS_TO_REFLECT_MAPS = [];
    analysisSettings.NORMAL_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION = -inf;

    % NORMAL BFO 180 settings
    analysisSettings.NORMAL_BFO180_ROTATIONS_DEG = [0, 180];
    analysisSettings.NORMAL_BFO180_FILTER_MIN_MAPS = BFO90_FILTER_MIN_MAPS;
    analysisSettings.NORMAL_BFO180_FILTER_MIN_TRIAL_CORRELATION = -inf; % No threshold
    analysisSettings.NORMAL_BFO180_CONTEXTS_TO_REFLECT_MAPS = [];
    analysisSettings.NORMAL_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION = -inf;

    % USED FOR THE SHUFFLED DISTRIBUTION.
    analysisSettings.NORMAL_BFO180_SIMILARITY_NSHUFFLES = 1; % or 0 to not perform any
    
    % Best Aligned settings
    analysisSettings.BESTALIGNED_COMPARISON_METHOD = 'correlation';
    analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA = 0.3;
    analysisSettings.BESTALIGNED_VERSION = 'v2'; % 'v2' or 'v3'
    
    % BFO 90 for Stability
    analysisSettings.STABILITY_BFO90_ROTATIONS_DEG = [0, 90, 180, 270];
    analysisSettings.STABILITY_BFO90_FILTER_MIN_MAPS = BFO90_FILTER_MIN_MAPS;
    analysisSettings.STABILITY_BFO90_FILTER_MIN_TRIAL_CORRELATION = -inf; % No threshold
    analysisSettings.STABILITY_BFO90_CONTEXTS_TO_REFLECT_MAPS = [];
    analysisSettings.STABILITY_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION = -inf;
    
    % STABILITY BFO 180 settings. Used for cumulative distributions.
    analysisSettings.STABILITY_BFO180_ROTATIONS_DEG = [0, 180];
    analysisSettings.STABILITY_BFO180_FILTER_MIN_MAPS = BFO90_FILTER_MIN_MAPS;
    analysisSettings.STABILITY_BFO180_FILTER_MIN_TRIAL_CORRELATION = -inf; % No threshold
    analysisSettings.STABILITY_BFO180_CONTEXTS_TO_REFLECT_MAPS = [];
    analysisSettings.STABILITY_BFO180_MIN_DIFFERENCE_FOR_MAX_CORRELATION = -inf;
    analysisSettings.STABILITY_BFO180_REFLECT_UNSTABLE_CONTEXT = []; % which context to reflect for unstable cells

    
    % This is the group label that the classification of stable vs unstable
    % will be used for. So we can classifiy only on a single day.
    analysisSettings.STABILITY_CLASSIFICATION_GROUP_LABEL = 'Day 3';
    analysisSettings.POPVECTORS_BESTALIGNED_STABILITY_THRESHOLD_CRITERIA = analysisSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA; % We dont have to use what is set for best aligned, we can use a separate one for the popvectors (if desired)
    analysisSettings.POPVECTOR_WITHIN_STABILITY_CLASSIFICATION_GROUP_LABEL = 'Day 1';
    analysisSettings.POPVECTOR_ACROSS_STABILITY_CLASSIFICATION_GROUP_LABEL = 'Day 3';

        
    analysisSettings.STABILITY_BFO90_REFLECT_UNSTABLE_CONTEXT = []; % which context to reflect for unstable cells

    
    % STABILITY BFO 90 PER CELL
    analysisSettings.BFO90_PER_CELL_MINIMUM_DATA = 4; % number of trial comparisons.
    
    % RATE MATRIX SETTINGS (only uses for tetrodes).
    analysisSettings.RATE_MATRIX_RATETYPE = 'mfrs';
    analysisSettings.RATE_MATRIX_NUM_TRIALS = 12;
    analysisSettings.RATE_MATRIX_FILTER_MIN_TRIALS = 4;
    % For averaging
    analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL = 'none';
    analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL = 'none';

    % BEHAVIOUR SETTINGS
    analysisSettings.ELIMINATE_NO_DIGS = true;
    

    
    % CALCIUM SPECIFIC
    analysisSettings.CALCIUM_APPLY_CELLLIKE_SPATIAL_FOOTPRINT_FILTER = true;
    analysisSettings.CALCIUM_ICS_FILTER_PERCENTILE = [];
    analysisSettings.CALCIUM_USE_ONLY_MANUAL_SFP_SCORES = false;
    
    % TETRODE SPECIFIC
    analysisSettings.TETRODES_INFORATE_FILTER_PERCENTILE = [];
    

end % function
