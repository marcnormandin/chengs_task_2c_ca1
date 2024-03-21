% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, NormalBFO180Similarity)

    uniqueDays = unique(DaysToSessions.dayNum);
    numDays = length(uniqueDays);
    
    T = [];
    k = 1;
    for iDay = 1:numDays
        dayNum = uniqueDays(iDay);
        Day = DaysToSessions(DaysToSessions.dayNum == dayNum,:);

        MatchTable = get_match_table_by_sessions(Day, NormalBFO180Similarity);
        
        if ~isempty(MatchTable)

            DayData = NormalBFO180Similarity(MatchTable.TableBRowIndex,:);

            T(k).dayNum = dayNum;
            T(k).dayLabel = Day.dayLabel{1}; % all should be the same

            % Get every field that ends in 'averageCorrelations'
            matchingInds = find(contains(DayData.Properties.VariableNames, 'averageCorrelations'));
            fieldNames = DayData.Properties.VariableNames(matchingInds);
            numNames = length(matchingInds);
            for iName = 1:numNames
                cellData = DayData.(fieldNames{iName});
                allCorrs = [];
                for j = 1:length(cellData)
                   allCorrs = cat(1, allCorrs, cellData{j});
                end
                T(k).(fieldNames{iName}) = allCorrs;

                [uz,cz] = ml_alg_cumdist(allCorrs);

                uz = reshape(uz, length(uz), 1);
                cz = reshape(cz, length(cz), 1);

                s = split(fieldNames{iName}, '_');
                prefix = s{1};
                T(k).(sprintf('%s_x', prefix)) = uz;
                T(k).(sprintf('%s_y', prefix)) = cz;
            end

            k = k + 1;
        end % ~isempty(MatchTable)
    end
    T = struct2table(T);
end % function

    