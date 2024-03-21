% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [J] = ml_alg_entropy_prep_image(I, NUM_BINS)
    % J will be discrete valued
    minI = min(I, [], 'all', 'omitnan');
    maxI = max(I, [], 'all', 'omitnan');
    
    J = (I - minI) ./ (maxI - minI) * NUM_BINS;
    
    J = floor(J);
end
