% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [randMaps, randIndices] = ml_util_shuffle_table_maps(maps)
    % This will shuffle a table of maps with no repetitions.

    numMaps = size(maps,1);

    randIndices = randperm(numMaps);
    randMaps = maps(randIndices);


    % Testing.
    % iMap = 1;
    % 
    % figure
    % subplot(1,2,1)
    % imagesc(maps{iMap})
    % subplot(1,2,2)
    % imagesc(randMaps{iMap})
end % function
