function [newModel, indUSE] = UpExp(model, reacsInd, UpMin, indFvar, indRvar)
% Integrate constraints for upregulated genes unidirectional
%
% USAGE:
%
%       [newModel, indUSE] = UpExp(model, reacsInd, UpMin, indFvar, indRvar)
%
% INPUTS:
%    model:           model with TFA structure
%    reacsInd:        reaction indexes in the .rxns field
%    UpMin:           lower bound required for upregulated rxns
%
%
% OPTIONAL INPUTS:
%    indFvar:         indexes of forward fluxes (default = find them)
%    indRvar:         indexes of reverse fluxes (default = find them)
%
% OUTPUTS:
%    newModel:        model with constraints
%    indUSE:          indexes of use variables associated to upregulated
%                     reactions (tagged with 'UP_')
%
% Daniel F. Hernandez 2016
% Anush Chiappino-Pepe 2017
%


if ~exist('indFvar', 'var') || isempty(indFvar)
    [indFvar,~] = getAllVar(model,{'F'});
end
if ~exist('indRvar', 'var') || isempty(indRvar)
    [indRvar,~] = getAllVar(model,{'R'});
end

% define use variables
intTag = {'UP'};
numVars = size(model.A,2);
indUSE = zeros(length(reacsInd),1);
for i=1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, '_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1)=0;
    indUSE(i,1)=numVars+i;
end

% define contraints for MILP
numCons = size(model.A,1);
% - F_rxn - R_rxn + UpMin*UP_rxn < 0
% if UP_rxn = 1 => F_rxn + R_rxn > UpMin 
% if UP_rxn = 0 => F_rxn + R_rxn > 0
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'UP1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indFvar(reacsInd(i))) = -1;
    model.A(numCons+i, indRvar(reacsInd(i))) = -1;
    model.A(numCons+i, indUSE(i)) = UpMin(i);
    model.constraintType(numCons+i) = {'<'};
end

% Adder Var
[numCons, numVars] = size(model.A);
model.varNames(numVars+1,1) = {'all_High'}; % name in UnifiedExpressionVar.m
model.var_ub(numVars+1,1) = Inf;
model.var_lb(numVars+1,1) = 0;
model.vartypes(numVars+1,1) = {'C'};
model.f(numVars+1,1) = 0;

% Adder Constraint
% sum(UP_rxn) = sum_ALLUSEUPGE
model.constraintNames(numCons+1) = {'sum_allHigh'};
model.rhs(numCons+1) = 0;
model.A(numCons+1, numVars+1) = -1;
model.A(numCons+1, indUSE) = 1;
model.constraintType(numCons+1) = {'='};

newModel = model;