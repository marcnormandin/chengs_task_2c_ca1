% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = helper_add_group_labels_to_table(T, SessionsToGroups, eliminateUnmatched)
    % This is a helper function to be used with tables that have an animal
    % name and session name. It will add the group labels and ids.
    
    MatchTable = get_match_table_by_sessions(T, SessionsToGroups);    

    % defaults
    T.groupId = nan(size(T,1),1);
    T.groupLabel = cell(size(T,1),1);
    
    T.groupId(MatchTable.TableARowIndex) = SessionsToGroups.groupId(MatchTable.TableBRowIndex);
    T.groupLabel(MatchTable.TableARowIndex) = SessionsToGroups.groupLabel(MatchTable.TableBRowIndex);

    if eliminateUnmatched
        T(isnan(T.groupId),:) = []; % eliminate data that is not part the sessions
    end
end % function
