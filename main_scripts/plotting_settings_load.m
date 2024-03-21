% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [plottingSettings] = plotting_settings_load(analysisSettings)
    plottingSettings = [];
    plottingSettings.NORMAL_BFO90_SHOW_EXCLUDED_COLUMN = false;
    plottingSettings.STABILITY_BFO90_SHOW_EXCLUDED_COLUMN = false;
    plottingSettings.RATE_MATRIX_CLIM = [];  % empty to use automatic ones

    plottingSettings.NORMAL_BFO180_SIMILARITY_SHUFFLED_SHOW = true;
    
    plottingSettings.OUTPUT_FOLDER = []; %workspace.OUTPUT_RESULTS_FOLDER;
    plottingSettings.FIGURE_FORMATS_TO_SAVE = {'png', 'svg', 'fig'};
    if analysisSettings.IS_CALCIUM_DATA
        plottingSettings.FIGURE_PREFIX = 'calcium_';
    else
        plottingSettings.FIGURE_PREFIX = 'tetrodes_';
    end
end % function
