function  final_val=evalEachRxn(model,rxn,ruleArray,geneRatio,f_and,f_or)
    % ruleArray is number of cells that are parsed with 'or' rules
    % and inside each cell there is another cell array than contains genes
    % which are seperated by and rule.
   
    if ((numel(ruleArray)==1) && (numel(ruleArray{1})==1))
        % this is for those genes which are catalized by one gene
        if (isKey(geneRatio,ruleArray{1}{1}))
            final_val = geneRatio(ruleArray{1}{1});
        else
            final_val = NaN;
        end
    elseif ((numel(ruleArray)==1) && (numel(ruleArray{1})>1)) %% this is the case when they are catlized by and only
        andVal = [];
        for j = 1:numel(ruleArray{1})
            % get value and apply operation
            if (isKey(geneRatio,ruleArray{1}{j}))
                andVal = [andVal; geneRatio(ruleArray{1}{j})];
            end
            
        end
        % this script is for ignoring nan values
        andVal = andVal(~isnan(andVal));
        if ~isempty(andVal)
            final_val = f_and(andVal);
            if numel(find(isnan(andVal)))>0
                disp('please check andVal contains nan (can not use geomean)')
            end
        else
            final_val = NaN;
        end
    else
        % this is for those catlized by or as well as and
        onVal = [];
        for i = 1:numel(ruleArray)
            andVal = [];
            for j = 1:numel(ruleArray{i})
                %%get value and apply operation
                if (isKey(geneRatio,ruleArray{i}{j}))
                    andVal = [andVal;geneRatio(ruleArray{i}{j})];
                end
            end
            % this script is for ignoring nan values
            andVal = andVal(~isnan(andVal));
            if (~isempty(andVal))
                % apply function of andRule
                
                onVal = [onVal;f_and(andVal)];
                if numel(find(isnan(andVal)))>0
                    disp('please check andVal contains nan (can not use geomean)')
                end
            end
        end
         % this script is for ignoring nan values
        onVal = onVal(~isnan(onVal));
        if (~isempty(onVal))
            final_val=f_or(onVal);
        else
            final_val=NaN;
        end
        
    end
end