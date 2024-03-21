% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_average_rate_differences_within_across_for_table(T)
    % Input should be a table that was created with one of the
    % compute average rate matrix tables. Each DAY should have ONE rate matrix.
    
    uniqueDays = unique(T.dayNum);
    numDays = length(uniqueDays);

    R = [];
    k = 1;
    for iDay = 1:numDays
        dayNum = uniqueDays(iDay);
        DayData = T(T.dayNum == dayNum, :);
        rateMatrixMean = DayData.rateMatrixMean{1}; 
        rateMatrixContextIds = DayData.rateMatrixContextIds(1,:);

        withinData = [];
        acrossData = [];
        for i = 1:length(rateMatrixContextIds)
            for j = 1+1:length(rateMatrixContextIds)
                x = rateMatrixMean(i,j);
                if rateMatrixContextIds(i) ~= rateMatrixContextIds(j)
                    acrossData = cat(1, acrossData, x);
                else
                    withinData = cat(1, withinData, x);
                end
            end
        end

        withinMean = nanmean(withinData, 'all');
        withinError = nanstd(withinData, 0, 'all') ./ sqrt(length(withinData));

        acrossMean = nanmean(acrossData, 'all');
        acrossError = nanstd(acrossData, 0, 'all') ./ sqrt(length(acrossData));

        R(k).dayNum = dayNum;
        R(k).dayLabel = DayData.dayLabel{1};
        R(k).withinMean = withinMean;
        R(k).withinError = withinError;
        R(k).acrossMean = acrossMean;
        R(k).acrossError = acrossError;

        k = k + 1;
    end
    R = struct2table(R);
end % function
