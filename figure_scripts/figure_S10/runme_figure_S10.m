% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure S10 (BF90 per context and per stability type)
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

% Calcium Results
calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calData = load(fullfile(calciumAnalysisResultsFolder, 'analysis_results.mat'));

% Process
run_script(calData.analysisSettings, calData.analysisResults.StabilityPerCell90, OUTPUT_FOLDER)

function run_script(analysisSettings, StabilityPerCell90, OUTPUT_FOLDER)

    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    eliminateUnmatched = true;
    
    
    StabilityPerCell90 = helper_add_group_labels_to_table(StabilityPerCell90, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
    StabilityPerCell90 = sortrows(StabilityPerCell90, 'groupId');
    StabilityPerCell90(~ismember(StabilityPerCell90.groupLabel, groupLabels),:) = [];
    
    % Add the histograms per cell since the percell doesn't have it.
    fieldNames2 = {'context1', 'context2'};
    
    % Init
    T = StabilityPerCell90;
    for iFn = 1:length(fieldNames2)
        T.(sprintf('%s_hist', fieldNames2{iFn})) = zeros(size(T,1),4);
        T.(sprintf('%s_prob', fieldNames2{iFn})) = zeros(size(T,1),4);
    end
        
    
    for iRow = 1:size(T,1)
        for iFn = 1:length(fieldNames2)
            fn = fieldNames2{iFn};
            inds = T.(sprintf('vind_%s', fn)){iRow};
            hc = histcounts(inds, 1:5);
    
            n = length(inds);
    
            T.(sprintf('%s_hist', fn))(iRow,:) = hc;
    
            T.(sprintf('%s_prob', fn))(iRow,:) = hc ./ n;
    
        end
    end
    
    TAll = T;
    TAll = sortrows(TAll, 'groupId');
    
    RUnstable = helper_reformat_table(TAll(TAll.isStable==false,:));
    RStable = helper_reformat_table(TAll(TAll.isStable==true,:));
    
    % Now save each table a sheet of an excel file
%     outputFilename = fullfile(OUTPUT_FOLDER, 'figure_S8_percelldata.xlsx');
%     writetable(RUnstable, outputFilename,'Sheet','UNSTABLE');
%     writetable(RStable, outputFilename,'Sheet','STABLE');
    
    %%
    R = [];
    R.Unstable = RUnstable;
    R.Stable = RStable;
    
    days = {'Day 1', 'Day 2', 'Day 3'};
    stabilities = {'Stable', 'Unstable'};
    cellTypes = {'FI', 'FS'};
    contexts = {'A', 'B'};
    angles = [0, 90, 180, 270];
    
    MAll = [];
    SEAll = [];
    for iDay = 1:length(days)
        day = days{iDay};
        for iStability = 1:length(stabilities)
            stability = stabilities{iStability};
            
            M = zeros(4,2);
            SE = zeros(4,2);
    
            for iContext = 1:2
                context = contexts{iContext};
                
                for iAngle = 1:4
                    angle = angles(iAngle);
    
                    r = R.(stability);
                    r = r(ismember(r.day, day) & ismember(r.context, context) & r.rotationDeg == angle, :);
    
                    x = mean(r.prob, 'all', 'omitnan');
                    s = std(r.prob,1,'all','omitnan');
                    n = size(r,1);
                    se = s ./ sqrt(n);
    
                    M(iAngle, iContext) = x;
                    SE(iAngle, iContext) = se;
                end
            end
    
            MAll{iDay,iStability} = M;
            SEAll{iDay, iStability} = SE;
        end
    end
    
    
    hFig = figure('position', get(0, 'screensize'));
    k = 1;
    for iDay = 1:3
        for iStability = 1:2
            ax(k) = subplot(3,2,k);
            k = k + 1;
    
            % Convert from fraction to percent
            Y = MAll{iDay,iStability} * 100;
            E = SEAll{iDay,iStability} * 100;
    
            colours4 = colour_matrix_blue_red_yellow_purple();
            b = barplot_with_errors(Y', E', colours4);
        
    
            xticks(1:4);
            xticklabels({['0' char(176)], ['90' char(176)], ['180' char(176)], ['270' char(176)]});
            
            grid on
            ax = gca;
            ax.FontSize = 15;
                        
            title(sprintf('%s %s', cellTypes{iStability}, days{iDay}))
            legend({'Context A', 'Context B'})
        end
        
        
    end
    linkaxes(ax, 'xy');
    
    fnPrefix = 'figure_S10_BFO90_per_context_stability';
    ml_savefig(hFig, OUTPUT_FOLDER, fnPrefix, {'png', 'svg', 'fig'});
end % function

function R = helper_reformat_table(TAll)
    R = [];
    nMS = size(TAll,1);
    for iMS = 1:nMS
        % Process a row
        TRow = TAll(iMS,:);
        R1 = helper_reformat_row(TRow, 1);
        R2 = helper_reformat_row(TRow, 2);
        R = [R; R1; R2];
    end
end % function
 
function R = helper_reformat_row(TRow, contextId)
    R = [];

    fn1 = sprintf('context%d_prob', contextId);
    fn2 = sprintf('context%d_hist', contextId);

    nAngles = length(TRow.rotations);

    for k = 1:nAngles
        R(k).subject = TRow.animalName{1};
        R(k).day = TRow.groupLabel{1};
        R(k).animalSessionCellName = TRow.animalSessionCellName{1};
        if contextId == 1
            R(k).context = 'A';
        else
            R(k).context = 'B';
        end
        R(k).rotationDeg = TRow.rotations(k);
        R(k).prob = TRow.(fn1)(k);
        R(k).hist = TRow.(fn2)(k);
    end
    R = struct2table(R);
end
