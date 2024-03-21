% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_average_rate_differences_within_across_bars(T)
    R = compute_average_rate_differences_within_across_for_table(T);

    hFig = plot_average_rate_differences_within_across_bars_table(R);

end % function
