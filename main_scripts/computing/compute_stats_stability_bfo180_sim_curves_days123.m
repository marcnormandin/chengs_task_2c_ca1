% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_stats_stability_bfo180_sim_curves_days123(analysisSettings, analysisResultsGroupAverages)
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];
    DaysToSessions = SessionsToGroups;

    StabilityBFO180Similarity = analysisResultsGroupAverages.StabilityBFO180Similarity;
    
    TStable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, StabilityBFO180Similarity.stable);
    TUnstable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, StabilityBFO180Similarity.unstable);
    
    TStable.isStable = true(size(TStable,1),1);
    TUnstable.isStable = false(size(TUnstable,1),1);
    
    T = cat(1, TStable, TUnstable);
    T(:, ~ismember(T.Properties.VariableNames, {'dayNum', 'dayLabel', 'isStable', 'within_averageCorrelations', 'across_averageCorrelations'})) = [];
    
    
    dayLabels = unique(T.dayLabel);
    nDays = length(dayLabels);
    
    stability = [true, false];
    types = {'within_averageCorrelations', 'across_averageCorrelations'};
    
    results = [];
    k = 1;
        
    for iDay = 1:nDays
        dayLabel = dayLabels{iDay};
        for iStable1 = 1:2
            for iStable2 = 1:2
                for iType1 = 1:2
                    for iType2 = 1:2
                        type1 = types{iType1};
                        type2 = types{iType2};
                        stable1 = stability(iStable1);
                        stable2 = stability(iStable2);
                        
                        ind1 = find(ismember(T.dayLabel, dayLabel) & T.isStable == stable1);
                        ind2 = find(ismember(T.dayLabel, dayLabel) & T.isStable == stable2);
                        
                        x1 = T.(type1)(ind1);
                        x1 = x1{1};
                        
                        x2 = T.(type2)(ind2);
                        x2 = x2{1};
                        
                        results(k).dayLabel = dayLabel;
                        if stable1
                            s1 = 'stable';
                        else
                            s1 = 'unstable';
                        end
                        if stable2
                            s2 = 'stable';
                        else
                            s2 = 'unstable';
                        end
                        
                        [H,p] = kstest2(x1, x2);
                        
                        results(k).curve_1 = sprintf('%s_%s', type1(1:6), s1);
                        results(k).curve_2 = sprintf('%s_%s', type2(1:6), s2);
                        results(k).kstest2_H = H;
                        results(k).kstest2_p = p;
                        results(k).num_samples_1 = length(x1);
                        results(k).num_samples_2 = length(x2);
                        k = k + 1;
                    end
                end
            end
        end
    end
    R = struct2table(results);
end % function

    
    
    
    
    