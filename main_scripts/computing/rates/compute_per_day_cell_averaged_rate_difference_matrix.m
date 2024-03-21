% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = compute_per_day_cell_averaged_rate_difference_matrix(DaysToSessions, RateMatrices, perCellNormalizationMethod)
    uniqueDays = unique(DaysToSessions.dayNum);
    numDays = length(uniqueDays);
    
    T = [];
    k = 1;
    for iDay = 1:numDays
        dayNum = uniqueDays(iDay);
        Day = DaysToSessions(DaysToSessions.dayNum == dayNum,:);

        MatchTable = get_match_table_by_sessions(Day, RateMatrices);

        RM = RateMatrices(MatchTable.TableBRowIndex,:);
        rateMatrices = [];
        for iM = 1:size(RM,1)
            rm = RM.rateMatrix{iM};
            
            if strcmpi(perCellNormalizationMethod, 'minmax')
                s = (rm - min(rm,[], 'all')) ./ (max(rm, [], 'all') - min(rm, [], 'all'));
            elseif strcmpi(perCellNormalizationMethod, 'none')
                s = rm;
            elseif strcmpi(perCellNormalizationMethod, 'totalsum')
                s = rm ./ nansum(rm, 'all');
            else
                error('Invalid normalizationMethod');
            end
            rateMatrices = cat(3, rateMatrices, s);
        end

        rateMatrixTrialIds = RM.rateMatrixTrialIds(1,:);
        rateMatrixContextIds = RM.rateMatrixContextIds(1,:);

        dayRateMatrix = nanmean(rateMatrices, 3);
        % Can also include the standard error
        
        T(k).dayNum = dayNum;
        T(k).dayLabel = Day.dayLabel{1}; % should all be the same
        
        T(k).rateMatrixMean = dayRateMatrix;
        
        T(k).perCellNormalizationMethod = perCellNormalizationMethod;
        T(k).numCellsAveraged = size(rateMatrices,3);
        
        T(k).numAnimalsAveraged = length(Day.animalName);
        T(k).animalNames = Day.animalName;
        T(k).sessionNames = Day.sessionName;
        T(k).rateMatrixTrialIds = rateMatrixTrialIds;
        T(k).rateMatrixContextIds = rateMatrixContextIds;
        
        k = k + 1;
    end
    T = struct2table(T);
end % function
