function  newRuleArray=newRuleEachRxn(model,rxn,ruleArray,geneRatio)

    if ((numel(ruleArray)==1) && (numel(ruleArray{1})==1)) % rxns associated to one gene
        newRuleArray = ruleArray;
    elseif ((numel(ruleArray)==1) && (numel(ruleArray{1})>1)) %% this is the case when they are catlized by and only   
        % exclude low expressed genes
        newRuleArray = ruleArray;
    else
       % this is for those catlized by or as well as and 
       ggg = {};
       for i = 1:numel(ruleArray)
            % this for loop for seperation with or
            gg = {};
            for j = 1:numel(ruleArray{i})
                % this for loop is seperation with and
                if ((isKey(geneRatio,ruleArray{i}{j})) && (~isequal(geneRatio(ruleArray{i}{j}),0)))
                    gg{end+1} = ruleArray{i}{j};
                end
            end
            if ~isempty(gg)
             ggg{end+1} = gg;
            end
       end
       newRuleArray = ggg;
    end
end