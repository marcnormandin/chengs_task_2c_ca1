% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [h] = plot_cell_rate_difference_matrix_by_context(rateMatrix, rateMatrixTrialIds, rateMatrixContextIds)

    if size(rateMatrix,1) ~= size(rateMatrix,2)
        error('Rate matrix must be square.');
    end
    
    rateMatrixNumTrials = size(rateMatrix,1);
    
    % Make labels for the rows and columns
    rateMatrixLabels = cell(rateMatrixNumTrials,1);
    for i = 1:rateMatrixNumTrials
        rateMatrixLabels{i} = sprintf('C%d-T%0.2d', rateMatrixContextIds(i), rateMatrixTrialIds(i));
    end
    
    numContext1 = sum(rateMatrixContextIds==1);
    numContext2 = sum(rateMatrixContextIds==2);
    
    % Matrix plot
    Z = nan(size(rateMatrix,1)+1, size(rateMatrix,2)+1);
    Z(1:size(rateMatrix,1), 1:size(rateMatrix,2)) = rateMatrix;
    rectangleColour = [0, 0, 0];
    
    h = pcolor(Z);
    shading flat
    axis equal tight
    colormap jet
    colorbar
    xticks(1.5:size(rateMatrix,1)+0.5) % center the labels
    yticks(1.5:size(rateMatrix,2)+0.5)
    xticklabels(rateMatrixLabels);
    yticklabels(rateMatrixLabels);
    xtickangle(90)
    set(gca, 'ydir', 'reverse')
    % Draw a rectangle aroung only context 1 with context 1
    rectangle('position', [1, 1, numContext1, numContext1], 'Curvature', 0, 'LineWidth', 4, 'EdgeColor', rectangleColour)
    % Draw a rectangle aroung only context 2 with context 2
    rectangle('position', [1+numContext1, 1+numContext1, numContext2, numContext2], 'Curvature', 0, 'LineWidth', 4, 'EdgeColor', rectangleColour)
end % function
