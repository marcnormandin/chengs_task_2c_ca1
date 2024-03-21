% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_whatever_we_need_per_map(occupancyMap, perimeterPadding_bins)
R = [];

% Compute the fraction of bins visited (entire arena)
BW1 = occupancyMap > 0; % whether or not the bin was visited at least once
BW1 = BW1(1+perimeterPadding_bins:(end-perimeterPadding_bins),:);
BW1 = BW1(1:end, 1+perimeterPadding_bins:(end-perimeterPadding_bins));
occupied = sum(BW1, "all");
total = numel(BW1);
R.fractionOccupied = occupied / total;

% Now compute how much of the cup bins were occupied
BW2 = occupancyMap  > 0; % don't apply padding because it will screw up the cup locations
ARENA_WIDTH_CM = 20;
ARENA_HEIGHT_CM = 30;
ARENA_WIDTH_BIN = size(occupancyMap,2);
ARENA_HEIGHT_BIN = size(occupancyMap,1);

BWCUPS = zeros(size(occupancyMap));
CUP_RADIUS_CM = 2.25; % Marc/I measured these.
CUP_CENTER_CM = 2.9 + CUP_RADIUS_CM;
% The cups all have the same size and offset, but we can have a different
% bin scaling in x and y (but they are actually the same)

CUP_CENTERS_CM = [ ...
    CUP_CENTER_CM, CUP_CENTER_CM;
    CUP_CENTER_CM + ARENA_WIDTH_CM/2, CUP_CENTER_CM;
    CUP_CENTER_CM, CUP_CENTER_CM + ARENA_HEIGHT_CM/2;
    CUP_CENTER_CM + ARENA_WIDTH_CM/2, CUP_CENTER_CM + ARENA_HEIGHT_CM/2];

%2023-09-28. Added a binary array for if a cup was visited at all or not.
CUPS_VISITED = false(1, size(CUP_CENTERS_CM,1));

% for every bin of the map
for i = 1:size(BWCUPS,1)
    for j = 1:size(BWCUPS,2)
        % create a small grid for the bin
        cx = (j-1) * ARENA_WIDTH_CM / ARENA_WIDTH_BIN + ARENA_WIDTH_CM / ARENA_WIDTH_BIN/2;
        cy = (i-1) * ARENA_WIDTH_CM / ARENA_WIDTH_BIN + ARENA_WIDTH_CM / ARENA_WIDTH_BIN/2;
        [XX,YY] = meshgrid(...
            linspace(cx - ARENA_WIDTH_CM / ARENA_WIDTH_BIN/2, cx + ARENA_WIDTH_CM / ARENA_WIDTH_BIN/2, 11), ...
            linspace(cy - ARENA_HEIGHT_CM / ARENA_HEIGHT_BIN/2, cy + ARENA_HEIGHT_CM / ARENA_HEIGHT_BIN/2, 11));
        
        % for every point of the small grid
        for k = 1:numel(XX)
            xx = XX(k);
            yy = YY(k);
            % for every cup
            for iCup = 1:size(CUP_CENTERS_CM,1)
                % check the distance, if grid point is within the radius of
                % the cup the mark it
                d = sqrt( (xx - CUP_CENTERS_CM(iCup,1))^2 + (yy - CUP_CENTERS_CM(iCup,2))^2);
                if d <= CUP_RADIUS_CM
                    BWCUPS(i,j) = 1;

                    % overwriting is fine
                    CUPS_VISITED(iCup) = true;
                end
            end
        end
    end
end

% Now compute how many bins the cups represent
numCupBins = sum(BWCUPS,'all');
numCupBinsOccupied = sum(BWCUPS .* BW2, 'all');
R.fractionCupsOccupied = numCupBinsOccupied / numCupBins;

% Now just put the fraction of cups that were visted (0,1,2,3,4) as a
% fraction
R.numCupsVisited = sum(CUPS_VISITED);

end % function
