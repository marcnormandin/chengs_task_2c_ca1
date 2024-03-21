% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [b] = barplot_with_errors(Y,E, colours)
    Y = Y';
    E = E';
    
    x = 1:size(Y,1);
    
    numGroups = size(Y,1);
    numBars = size(Y,2);
    
    b = bar(x, Y, 'grouped', 'FaceColor','flat');
    for k = 1:size(Y,2)
        b(k).CData = colours(k,:);
    end
    
    hold on
    groupWidth = min(0.8, numBars / (numBars + 1.5));
    for i = 1:numBars
        xx = (1:numGroups) - groupWidth/2 + (2*i-1)*groupWidth / (2*numBars);
        errorbar(xx, Y(:,i), E(:,i), 'k.', 'LineWidth', 2);
    end
end % function
