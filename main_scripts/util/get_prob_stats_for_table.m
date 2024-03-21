% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T,E] = get_prob_stats_for_table(R, probNames)

    T = [];
    E = [];
    
    for iField = 1:length(probNames)
        fn = probNames{iField};
        
        % We have to treat the excluded data differently so check
        if contains(fn, 'excluded')
            X = R.(fn);
            % Isabel said that if we have all zeros for the angles, then do
            % not include it.
            indEmpty = find(all(X == 0, 2));
            if ~isempty(indEmpty)
                X(indEmpty,:) = [];
            end
            y = mean(X, 1);
            e = std(X, 0, 1) ./ sqrt(size(X,1));
        else
            X = R.(fn);
            % Isabel said that if we have all zeros for the angles, then do
            % not include it.
            indEmpty = find(all(X == 0, 2));
            if ~isempty(indEmpty)
                X(indEmpty,:) = [];
            end

            y = mean(X, 1);
            % Normalize just in case because if a sessions data has no data
            % (all zeros, then it screws up), or eliminate.
            y = y ./ sum(y, 'all') * 100;

            e = std(X, 0, 1) ./ sqrt(size(X,1));
        end
        
        E.(fn) = e;
        T.(fn) = y;
    end
    T = struct2table(T);
    E = struct2table(E);
end % function
