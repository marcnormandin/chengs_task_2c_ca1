% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [output] = compute_cell_rate_difference_matrix_by_context(firingRates, trialIds, contextIds, rateMatrixNumTrials)
    if rateMatrixNumTrials <= 0 || length(rateMatrixNumTrials) ~= 1
        error('rateMatrixNumTrials must be greater than zero');
    end
    
    if length(unique(trialIds)) ~= length(trialIds)
        error('Trial ids are repeated!');
    end
    
    if mod(rateMatrixNumTrials, 2) == 1
        error('rateMatrixNumTrials must be even.');
    end
        
    rateMatrixContextIds = repmat([1,2], rateMatrixNumTrials/2,1); % will only work if context is 1 or 2
    rateMatrixContextIds = rateMatrixContextIds(:); % 1 1 1 1 1 1 2 2 2 2 2 2
    rateMatrixContextIds = reshape(rateMatrixContextIds, 1, length(rateMatrixContextIds));
    
    % Sort the rate matrix rows and cols by context, so that we can see
    % patterns based on context more easily.
    % Get the unique combinations of (trial id, context id)
    urows = unique([trialIds, contextIds], 'rows');
    [surows, ~] = sortrows(urows, 2);
    
    context1Pairs = surows(surows(:,2)==1,:);
    context2Pairs = surows(surows(:,2)==2,:);
    
    hasContext1 = ~isempty(context1Pairs);
    hasContext2 = ~isempty(context2Pairs);
    
    if hasContext1
        % Make sure that all trial ids are uniformly even or odd
        tids = context1Pairs(:,1);
        context1IsEven = all(mod(tids,2)==0);
    end
    
    if hasContext2
        % Make sure that all trial ids are uniformly even or odd
        tids = context2Pairs(:,1);
        context2IsEven = all(mod(tids,2)==0);
    end
    
    % possibililty A: Context 1 has odd-valued trials and Context 2 has
    % even-valued trials.
    rateMatrixTrialIdsA = [1:2:rateMatrixNumTrials, 2:2:rateMatrixNumTrials]; 
        
    % possibility B: Context 1 has even-valued trials and Context 2 has
    % odd-valued trials.
    rateMatrixTrialIdsB = [2:2:rateMatrixNumTrials, 1:2:rateMatrixNumTrials];
    
    rateMatrixTrialIds = nan;
    
    if hasContext1 && hasContext2
        if context1IsEven && ~context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        elseif ~context1IsEven && context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        else
            error('Invalid contexts do not contain uniformly even or odd trial ids')
        end
    elseif hasContext1
        if context1IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        else
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        end
    elseif hasContext2
        if context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        else
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        end
    else
        error('Something is messed up with the context ids')
    end
        
    % Create and compute elements for the matrix
    rateMatrix = nan(rateMatrixNumTrials, rateMatrixNumTrials);
    numRates = length(firingRates);
    for k1 = 1:numRates
      for k2 = 1:numRates
          trial1 = trialIds(k1); % Get the actual trial id
          trial2 = trialIds(k2); % Get the actual trial id
          
          if trial1 == trial2 % We dont care about trial with itself
              continue;
          end

          rateDiff = abs( firingRates(k1) - firingRates(k2) );

          rmrow = find(rateMatrixTrialIds == trial1);
          rmcol = find(rateMatrixTrialIds == trial2);
          
          if length(rmrow) ~= 1 || length(rmcol) ~= 1
              error('Data error. A trial maps to more than one row.');
          end
          
          rateMatrix(rmrow, rmcol) = rateDiff;
      end
    end
        
    output.numRates = numRates;
    output.rateMatrix = rateMatrix;
    output.rateMatrixTrialIds = rateMatrixTrialIds;
    output.rateMatrixContextIds = rateMatrixContextIds;
end % function
