% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [X] = helper_percell90_add_stats(X)

% Add distributions per cell
numCells = size(X,1);
numRotations = 4;
X.hc_within = zeros(numCells, numRotations);
X.prob_within = zeros(numCells, numRotations);
X.samples_within = zeros(numCells,1);

X.hc_across = zeros(numCells, numRotations);
X.prob_across = zeros(numCells, numRotations);
X.samples_across = zeros(numCells,1);


X.hc_any = zeros(numCells, numRotations);
X.prob_any = zeros(numCells, numRotations);
X.samples_any = zeros(numCells,1);

for iCell = 1:numCells
    % within
    vind = [X.vind_context1{iCell}, X.vind_context2{iCell}];
    hc = histcounts(vind, 1:(numRotations+1));
    p = hc ./ sum(hc);
    X.hc_within(iCell,:) = hc;
    X.prob_within(iCell,:) = p;
    X.samples_within(iCell) = sum(hc);
    
    % across
    vind = X.vind_different{iCell};
    hc = histcounts(vind, 1:(numRotations+1));
    p = hc ./ sum(hc);
    X.hc_across(iCell,:) = hc;
    X.prob_across(iCell,:) = p;
    X.samples_across(iCell) = sum(hc);
    
    % any
    vind = X.vind_any{iCell};
    hc = histcounts(vind, 1:(numRotations+1));
    p = hc ./ sum(hc);
    X.hc_any(iCell,:) = hc;
    X.prob_any(iCell,:) = p;
    X.samples_any(iCell) = sum(hc);
end
end % function
