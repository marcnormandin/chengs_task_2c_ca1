% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes figure S1, which is the dig matrices.
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

T = readtable(fullfile(INPUT_FOLDER, 'digs.xlsx'));

%% The data has multiple notations for digs. e.g. Feat and F and F(MAYBE),etc
T(cellfun(@(x)isempty(x), T.secondDig),:) = [];
T(cellfun(@(x)isempty(x), T.firstDig),:) = [];

delInds = [];
for iRow = 1:size(T,1)
    if strcmpi(T.secondDig{iRow}, 'None')
        delInds = [delInds, iRow];
    end
end
T(delInds,:) = [];

unique(T.firstDig)
unique(T.secondDig)

M = zeros(4,4);
for iRow = 1:size(T,1)
    sArray = split(T.secondDig{iRow},',');
    priorIndex = get_index(T.firstDig{iRow});

    for iSecond = 1:length(sArray)
        s = sArray(iSecond);
        sIndex = get_index(s);

        M(priorIndex,sIndex) = M(priorIndex,sIndex) + 1;
        priorIndex = sIndex;
    end
end


%%
MS = M(2:4,1:4);
hFig = figure;
heatmap(MS, 'XLabel', 'Next Dig', 'YLabel', 'Current Dig')
colormap jet

title('Subsequent Dig Analysis')
ax = gca;
ax.XData = ["C", "F", "G", "W"];
ax.YData = ["F", "G", "W"];

mulana_savefig(hFig, OUTPUT_FOLDER,'figure_S1_digmatrices', {'png', 'svg'});

%%
function index = get_index(s)
    if strcmpi(s, 'C') || strcmpi(s, 'Corr')
        index = 1;
    elseif strcmpi(s, 'F') || strcmpi(s, 'Feat')
        index = 2;
    elseif strcmpi(s, 'G') || strcmpi(s, 'Geo')
        index = 3;
    elseif strcmpi(s, 'W') || strcmpi(s, 'Wrong')
        index = 4;
    else
        error('%s is not valid', s);
    end
end % function

