% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function main_analysis_group_averages_plot(analysisSettings, plottingSettings, analysisResultsGroupAverages)
    plot_and_save_normal_bfo90_all(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    
    plot_and_save_stability_bfo90_animal_withinacross(analysisSettings, analysisResultsGroupAverages, plottingSettings);

    % These don't require cells to be registered so will work with calcium
    % and tetrode data.
    plot_and_save_stability_bfo180_sim_curves_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    
    if analysisSettings.IS_CALCIUM_DATA
        plot_and_save_stability_bfo180_reg_sim_curves_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    end
    
    if analysisSettings.IS_CALCIUM_DATA
        plot_and_save_stability_bfo90_reg_animal_withinacross(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    end
    

    
    %plot_and_save_stability_bfo90_reflected_withinacross(analysisSettings, analysisResultsGroupAverages, plottingSettings);

    plot_and_save_rate_matrices_and_hists_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    plot_and_save_bfo180_similarity_curves_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    plot_and_save_stable_cell_percentages_hist_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    plot_and_save_stable_rate_matrices_and_hists_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings);
    
end % function

%%
% For the plot we only want days 1,2,3
function plot_and_save_normal_bfo90_all(analysisSettings, analysisResultsGroupAverages, plottingSettings) 
    MeanTable = analysisResultsGroupAverages.NormalBFO90.MeanTable;
    ErrorTable = analysisResultsGroupAverages.NormalBFO90.ErrorTable;

    MeanTable(~ismember(MeanTable.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
    ErrorTable(~ismember(ErrorTable.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];

    hFig = bfo90table_plot_all_prob(MeanTable, ErrorTable, plottingSettings.NORMAL_BFO90_SHOW_EXCLUDED_COLUMN);
    
    fnPrefix = sprintf("%snormal_bfo90_all_days123", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);
end

function plot_and_save_stability_bfo90_animal_withinacross(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    % This plots per animal
    
    % PLOT 1
    MeanTable = analysisResultsGroupAverages.StabilityBFO90.MeanTable;
    ErrorTable = analysisResultsGroupAverages.StabilityBFO90.ErrorTable;
    MeanTable(~ismember(MeanTable.groupId, [1,2,3]),:) = [];
    ErrorTable(~ismember(ErrorTable.groupId, [1,2,3]),:) = [];

    
    [hFig] = stabilitybfo90table_plot_withinacross_prob(MeanTable, ErrorTable, plottingSettings.STABILITY_BFO90_SHOW_EXCLUDED_COLUMN);
    fnPrefix = sprintf("%sstability_bfo90_animal_withinacross_days123", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);

    % PLOT 2
    % For paper we want day 3
    MeanTable = analysisResultsGroupAverages.StabilityBFO90.MeanTable;
    ErrorTable = analysisResultsGroupAverages.StabilityBFO90.ErrorTable;
    MeanTable(MeanTable.groupId ~= 3, :) = [];
    ErrorTable(ErrorTable.groupId ~= 3, :) = [];
    [hFig] = stabilitybfo90table_plot_withinacross_prob(MeanTable, ErrorTable, plottingSettings.STABILITY_BFO90_SHOW_EXCLUDED_COLUMN);
    fnPrefix = sprintf("%sstability_bfo90_animal_withinacross_day3", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);
end % function

function plot_and_save_stability_bfo90_reg_animal_withinacross(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    % This plots per animal REGISTERED
    numGroups = length(analysisResultsGroupAverages.StabilityBFO90_Reg);
    for iGroup = 1:numGroups
        
        STABILITY_CLASSIFICATION_GROUP_LABEL = analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL;
        
        % PLOT 1
        MeanTable = analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).MeanTable;
        ErrorTable = analysisResultsGroupAverages.StabilityBFO90_Reg(iGroup).ErrorTable;
        MeanTable(~ismember(MeanTable.groupId, [1,2,3]),:) = [];
        ErrorTable(~ismember(ErrorTable.groupId, [1,2,3]),:) = [];


        [hFig] = stabilitybfo90table_plot_withinacross_prob(MeanTable, ErrorTable, plottingSettings.STABILITY_BFO90_SHOW_EXCLUDED_COLUMN);
        fnPrefix = sprintf("%sstability_bfo90_reg_%s_animal_withinacross_days123", plottingSettings.FIGURE_PREFIX, strrep(STABILITY_CLASSIFICATION_GROUP_LABEL, ' ', '_'));
        sgtitle(sprintf('REGISTERED (Classified: %s)', STABILITY_CLASSIFICATION_GROUP_LABEL));
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);

        % PLOT 2
        % For paper we want day 3
        MeanTable = analysisResultsGroupAverages.StabilityBFO90_Reg.MeanTable;
        ErrorTable = analysisResultsGroupAverages.StabilityBFO90_Reg.ErrorTable;
        MeanTable(MeanTable.groupId ~= 3, :) = [];
        ErrorTable(ErrorTable.groupId ~= 3, :) = [];
        [hFig] = stabilitybfo90table_plot_withinacross_prob(MeanTable, ErrorTable, plottingSettings.STABILITY_BFO90_SHOW_EXCLUDED_COLUMN);
        sgtitle(sprintf('REGISTERED (Classified: %s)', STABILITY_CLASSIFICATION_GROUP_LABEL));
        fnPrefix = sprintf("%sstability_bfo90_reg_%s_animal_withinacross_days3", plottingSettings.FIGURE_PREFIX, strrep(STABILITY_CLASSIFICATION_GROUP_LABEL, ' ', '_'));
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end % iGroup
    
end % function


function plot_and_save_rate_matrices_and_hists_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    if ~analysisSettings.IS_CALCIUM_DATA
        % Plotting
        RateMatricesAveraged = analysisResultsGroupAverages.RateMatrices;
        % Plot only days 1,2,3
        RateMatricesAveraged(~ismember(RateMatricesAveraged.groupId, [1,2,3]),:) = [];
        
        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RateMatricesAveraged, plottingSettings.RATE_MATRIX_CLIM);
        
        fnPrefix = sprintf("%srate_matrices_averaged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end

    % Rate bars
    if ~analysisSettings.IS_CALCIUM_DATA
        RateMatricesAveragedHistogram = analysisResultsGroupAverages.RateMatricesHistogram;
        hFig = plot_average_rate_differences_within_across_bars_table(RateMatricesAveragedHistogram);        
        fnPrefix = sprintf("%srate_hist_averaged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end

end % function

function plot_and_save_bfo180_similarity_curves_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];
    
    shuffledTable = analysisResultsGroupAverages.NormalBFO180SimilarityShuffled;
    if plottingSettings.NORMAL_BFO180_SIMILARITY_SHUFFLED_SHOW ~= true
        shuffledTable = [];
    end

    [hFig, ax] = plot_bfo180_average_similarity_curves_per_period(SessionsToGroups, analysisResultsGroupAverages.NormalBFO180Similarity, shuffledTable);

    fnPrefix = sprintf("%sbfo180_similarity_curves_averaged_days123", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);
end % function

function plot_and_save_stable_cell_percentages_hist_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    GroupCounts = analysisResultsGroupAverages.StableCells;
    colours = colour_matrix_blue_red_yellow_purple();
    
    hFig = figure();
    [b] = barplot_with_errors(GroupCounts.meanPercentStable'./100,GroupCounts.errorPercentStable'./100, colours);
    set(gca, 'xtick', 1:size(GroupCounts,1));
    set(gca, 'xticklabels', GroupCounts.groupLabel)
    ylabel('Stable Cell / Total Ratio')
    grid on
    
    fnPrefix = sprintf("%sbestaligned_stable_counts_hist_averaged_days123", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);
end % function

function plot_and_save_stable_rate_matrices_and_hists_averaged_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    if ~analysisSettings.IS_CALCIUM_DATA
        RMAverageStable = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStable;
        RMAverageUnstable = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstable;
        RMAverageStableWithinAcross = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStableWithinAcross;
        RMAverageUnstableWithinAcross = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstableWithinAcross;

        % Plot
        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RMAverageStable, plottingSettings.RATE_MATRIX_CLIM);
        sgtitle(sprintf('Stable (Per animal (%s, %s))', analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));    
        fnPrefix = sprintf("%sstable_rate_matrices_averaged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);


        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RMAverageUnstable, plottingSettings.RATE_MATRIX_CLIM);
        sgtitle(sprintf('Unstable (Per animal (%s, %s))', analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sunstable_rate_matrices_averaged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);

        % Hist 1
        hFig = plot_average_rate_differences_within_across_bars_table(RMAverageStableWithinAcross);
        sgtitle(sprintf('Stable (Per animal (%s, %s))', analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sstable_rate_matrices_averaged_hist_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);


        % Hist 2
        hFig = plot_average_rate_differences_within_across_bars_table(RMAverageUnstableWithinAcross);
        sgtitle(sprintf('Unstable (Per animal (%s, %s))', analysisSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, analysisSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sunstable_rate_matrices_averaged_hist_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end % not calcium
end % function




%% STABILITY CUMULATIVE CURVE FUNCTIONS
function plot_and_save_stability_bfo180_reg_sim_curves_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];
    
    numGroups = length(analysisResultsGroupAverages.Stability180_Reg);
    for iGroup = 1:numGroups
        STABILITY_CLASSIFICATION_GROUP_LABEL = analysisResultsGroupAverages.Stability180_Reg(iGroup).STABILITY_CLASSIFICATION_GROUP_LABEL;
        
        StabilityBFO180Similarity_Reg = analysisResultsGroupAverages.Stability180_Reg(iGroup).StabilityBFO180Similarity_Reg;
        
        hFig = plot_stability_bfo180_average_similarity_curves_per_period(SessionsToGroups, StabilityBFO180Similarity_Reg);
        sgtitle(sprintf('Cumulative Similarity\n(REGISTERED) (Classified on %s)', STABILITY_CLASSIFICATION_GROUP_LABEL))

        fnPrefix = sprintf("%sstability_reg_%s_bfo180_similarity_curves_averaged_days123", plottingSettings.FIGURE_PREFIX, strrep(STABILITY_CLASSIFICATION_GROUP_LABEL, ' ', '_'));
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end % iGroup
end % function

function plot_and_save_stability_bfo180_sim_curves_days123(analysisSettings, analysisResultsGroupAverages, plottingSettings)
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];

    hFig = plot_stability_bfo180_average_similarity_curves_per_period(SessionsToGroups, analysisResultsGroupAverages.StabilityBFO180Similarity);
    
    fnPrefix = sprintf("%sstability_bfo180_similarity_curves_averaged_days123", plottingSettings.FIGURE_PREFIX);
    ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
    close(hFig);
end % function

function [hFig] = plot_stability_bfo180_average_similarity_curves_per_period(DaysToSessions, StabilityBFO180Similarity)
    TStable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, StabilityBFO180Similarity.stable);
    TUnstable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, StabilityBFO180Similarity.unstable);
    
    uniqueDays = unique(TStable.dayNum);
    numDays = length(uniqueDays);

    hFig = figure('position', get(0, 'screensize'));

    typeNamesToPlot = {'within', 'across'};

        coloursYellowPurple = [
            0.929, 0.694, 0.125; ...
            0.494, 0.184, 0.556];

    %     coloursBlueRed = [...
    %             0, 0.4470, 0.7410; ...
    %             0.8500, 0.3250, 0.0980];

    colours = coloursYellowPurple;

    ax = [];
    for iDay = 1:numDays
        ax(iDay) = subplot(1,numDays,iDay);
        
        dayNum = uniqueDays(iDay);
        
        help_the_helper(typeNamesToPlot, colours, true, TStable, dayNum);
        help_the_helper(typeNamesToPlot, colours, false, TUnstable, dayNum);
        
        % Make the legend
        l = cell(2*length(typeNamesToPlot),1);
        for i = 1:length(typeNamesToPlot)
            x = typeNamesToPlot{i};
            l{i} = sprintf('stable %s', x);
        end
        for i = 1:length(typeNamesToPlot)
            x = typeNamesToPlot{i};
            l{i+length(typeNamesToPlot)} = sprintf('unstable %s', x);
        end

        legend(l, 'location', 'southoutside', 'orientation', 'vertical')
        grid on
        grid minor
    end
    linkaxes(ax, 'xy')
    sgtitle(sprintf('Cumulative Similarity\n(NO REGISTRATION USED)'))
end % function

function help_the_helper(typeNamesToPlot, colours, isStable, T, dayNum)
    DayData = T(T.dayNum == dayNum,:);

    
    for iType = 1:length(typeNamesToPlot)
       x = DayData.(sprintf('%s_x', typeNamesToPlot{iType})){:};
       y = DayData.(sprintf('%s_y', typeNamesToPlot{iType})){:};

       if isStable
        plot(x,y,'-', 'color', colours(iType,:), 'linewidth', 2)
       else
        plot(x,y,':', 'color', colours(iType,:), 'linewidth', 2)
       end
       
       hold on
       xlabel('Correlation')
       ylabel('Cumulative')
       axis equal tight
    end
    title(sprintf('%s', DayData.dayLabel{1}))
end % function
