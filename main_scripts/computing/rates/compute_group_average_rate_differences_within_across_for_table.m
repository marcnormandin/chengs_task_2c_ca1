% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_group_average_rate_differences_within_across_for_table(RateMatricesAveraged)
    numGroups = size(RateMatricesAveraged,1);
    
    R = [];
    k = 1;
    for iGroup = 1:numGroups
        groupId = RateMatricesAveraged.groupId(iGroup);
        groupLabel = RateMatricesAveraged.groupLabel{iGroup};
        rateMatrixMean = RateMatricesAveraged.rateMatrixMean{iGroup}; 
        rateMatrixContextIds = RateMatricesAveraged.rateMatrixContextIds(iGroup,:);

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
        withinError = nanstd(withinData, 1, 'all') ./ sqrt(length(withinData));

        acrossMean = nanmean(acrossData, 'all');
        acrossError = nanstd(acrossData, 1, 'all') ./ sqrt(length(acrossData));

        R(k).groupId = groupId;
        R(k).groupLabel = groupLabel;
        R(k).withinMean = withinMean;
        R(k).withinError = withinError;
        R(k).acrossMean = acrossMean;
        R(k).acrossError = acrossError;

        k = k + 1;
    end
    R = struct2table(R);
end % function
