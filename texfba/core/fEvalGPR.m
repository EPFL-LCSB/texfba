function [rxns, rxnExp, newgrRules, optionOut] = fEvalGPR(model, genes, ...
    geneRatio, f_and, f_or)
% Evaluate GRP based on gene expression data
%
% USAGE:
%
%       [rxns, rxnExp, newgrRules, optionOut] = fEvalGPR(model, genes, geneRatio, f_and, f_or)
%
% INPUTS:
%    model:           model with TFA structure
%    genes:           vector of genes in the model
%    geneRatio:       vector of corresponding gene ratios
%
%
% OPTIONAL INPUTS:
%    f_and:           (default = @min)
%    f_or:            (default = @max)
%
% OUTPUTS:
%    rxns:            mo
%    rxnExp:          mo
%    newgrRules:      new GPR rules modified from the original GPR (low
%                     expressed genes are excluded from complex GPR)
%    optionOut:       mo
%
%
% .. Author:
% Vikash Pandey 2016
% Anush Chiappino-Pepe 2017
%

if (nargin < 4)
    f_and = @min;
end
if (nargin < 5)
    f_or = @max;
end

% crate map from gene index to ratio
geneRatio_map = containers.Map(genes, geneRatio);

[GPRrules, ~] = parseGPRModify(model);
% extract info:
GPRrxns = keys(GPRrules)';
GPRdecomp = values(GPRrules)';
[y,r]=ismember(GPRrxns,model.rxns);
GPR4comparsion = model.grRules(r);

[rxns, rxnExp, newgrRules] = getRxnval(model, GPRrules, geneRatio_map, ...
    f_and, f_or);

vv=GPRrules.values;
ids=[];
genes={};

for i = 1:numel(newgrRules)
    if ~isequal(vv{i},newgrRules{i})
        ids = [ids;i];
    end
    for j = 1:numel(newgrRules{i})
        for k = 1:numel(newgrRules{i}{j})
            genes{end+1} = newgrRules{i}{j}{k};
        end
    end
end

optionOut.diffIds = ids;
optionOut.genes = unique(genes);
optionOut.GPRrules = GPRrules.values;

end


function [rxns, rxExp, newgrRules] = getRxnval(model, parseGPRs, ...
    geneRatio_map, f_and, f_or)
% model
% parseGPRs: reactions genes GPR
% geneRatio_map : geneRatio_map
rxns = parseGPRs.keys;
rxExp = nan(numel(rxns),1);
rxnWeight = nan(numel(rxns),1);
parseRules = parseGPRs.values;
newgrRules = {};
for r = 1:numel(rxns)
    newgrRule = newRuleEachRxn(model,rxns{r},parseRules{r},geneRatio_map);
    newgrRules{end+1} = newgrRule;
    rxExp(r) = evalEachRxn(model,rxns{r},parseRules{r},geneRatio_map,f_and,f_or);
end
end