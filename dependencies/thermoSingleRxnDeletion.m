function [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution, impactTasks] = thermoSingleRxnDeletion(model, method, rxnList, verbFlag, flagTasks, essThr, indNF)
% Performs single reaction deletion analysis using TFA, tMOMA or tlinearMOMA
%
% USAGE:
%
%    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution, impactTasks] = thermoSingleRxnDeletion(model, method, rxnList, verbFlag, flagTasks, essThr)
%
% INPUT:
%    model:           TFA model structure including reaction names
%
% OPTIONAL INPUTS:
%    method:          Either 'TFA', 'tMOMA', or 'tlMOMA' (Default = 'TFA')
%    rxnList:         List of reactions to be deleted (Default = all reactions)
%    verbFlag:        Verbose output (Default = false)
%    flagTasks:       Determine which BBBs cannot be produced upon knockout
%                     (default = false)
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for flagTasks (default = 0.1)
%
% OUTPUTS:
%    grRatio:         Computed growth rate ratio between deletion strain and wild type
%    grRateKO:        Deletion strain growth rates (1/h)
%    grRateWT:        Wild type growth rate (1/h)
%    hasEffect:       Does a reaction deletion affect anything
%    delRxn:          Deleted reaction
%    fluxSolution:    TFA/tMOMA/tlMOMA fluxes for `KO` strains
%    impactTasks:     list of BBBs and its production upon KO of rxnList
%
% .. Authors:
%       - Richard Que 12/04/2009 Based on singleGeneDeletion.m written by Markus Herrgard
%       - Karthik Raman 06/28/2017 Based on github.com/RamanLab/FastSL
%       - Anush Chiappino-Pepe 10/8/2017 thermo and checkBBBTasks
tic
if (nargin < 2)
    method = 'TFA';
end
if (nargin < 3)
    rxnList = model.rxns;
else
    if (isempty(rxnList))
        rxnList = model.rxns;
    end
end
if (nargin < 4)
    verbFlag = false;
end
if (nargin < 5)
    flagTasks = false;
end
if (nargin < 6)
    essThr = 0.1;
end
if (nargin < 7)
    indNF = getAllVar(model,{'NF'});
end

nDelRxns = length(rxnList);
impactTasks=cell(nDelRxns,2);
impactTasks(:,1)=rxnList;

if isempty(indNF)
    model = addNetFluxVariables(model);
    indNF = getAllVar(model,{'NF'});
end
% if strcmp(method,'tMOMA')
%     minmax = runMinMax(model,model.varNames(indNF));
%     minmax =roundn(minmax,-9);
%     model.var_lb(indNF) = minmax(:,1);
%     model.var_ub(indNF) = minmax(:,2);
% end

solWT = optimizeThermoModel(model); %adapt for minNorm
grRateWT = solWT.val;

if isempty(grRateWT) || grRateWT==0
    error('infeasible WT model: check bounds from minmax')
end

% Identify reactions that do not carry a flux in solWT; none of these can be lethal
% Identify J_z, the set of reactions that do not carry a flux in TFA
solWTtfa = optimizeThermoModel(model,true);
Jzrxns = model.rxns(solWTtfa.x(indNF)==0);

grRateKO = ones(nDelRxns, 1)*grRateWT;
hasEffect = true(nDelRxns, 1);
fluxSolution = repmat(solWT.x(indNF), 1, nDelRxns);
delRxn = columnVector(rxnList);
if (verbFlag)
    fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n', 'No', 'Perc', 'Name', 'Growth rate', 'Rel. GR');
end
showprogress(0,'Thermo single reaction deletion analysis in progress ...');
for i = 1:nDelRxns
    showprogress(i/nDelRxns);
    if ismember(rxnList{i}, Jzrxns)
	% If the reaction carries no flux in WT, deleting it cannot affect
	% the flux solution. Assign WT solution without solving LP.
        solKO = solWT;
        hasEffect(i) = false;
    else
        modelDel = model;
        modelDel.var_ub(ismember(model.varNames,strcat('F_',rxnList{i}))) = 0;
        modelDel.var_ub(ismember(model.varNames,strcat('R_',rxnList{i}))) = 0;
        modelDel.var_lb(ismember(model.varNames,strcat('F_',rxnList{i}))) = 0;
        modelDel.var_lb(ismember(model.varNames,strcat('R_',rxnList{i}))) = 0;
        modelDel.var_ub(ismember(model.varNames,strcat('NF_',rxnList{i}))) = 0;
        modelDel.var_lb(ismember(model.varNames,strcat('NF_',rxnList{i}))) = 0;
        switch method
            case 'tlMOMA'
                solKO = tlinearMOMA(model, modelDel, 'max'); %not yet adapted
            case 'tMOMA'
                solKO = tMOMA(model, modelDel, 'max', false, true, 'cplex'); %true minNorm
            otherwise
                solKO = optimizeThermoModel(modelDel);
        end
    end
    if (solKO.stat == 1)
        grRateKO(i) = solKO.val;
        fluxSolution(:, i) = solKO.x(indNF);
        if flagTasks
            if grRateKO(i) < essThr*solWT.val
                [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
            else
                impactTasks{i,2} = [];
            end
        end
    else
        grRateKO(i) = NaN;
        fluxSolution(:,i) = nan(length(indNF),1);
        if flagTasks
            [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
        end
    end
    if (verbFlag)
        fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n', i, 100*i/nDelRxns, rxnList{i}, grRateKO(i), grRateKO(i)/grRateWT*100);
    end
end

grRatio = grRateKO/grRateWT;
toc
end
