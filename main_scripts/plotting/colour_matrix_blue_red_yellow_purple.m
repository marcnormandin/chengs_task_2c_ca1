% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [colours] = colour_matrix_blue_red_yellow_purple()
    colours = cat(1, colour_matrix_blue_red(), colour_matrix_yellow_purple());
end
