function [newModel, indUSE] = DownExp(model, reacsInd, lowMax, indFvar, indRvar)
% Integrate constraints for downregulated genes unidirectional
%
% USAGE:
%
%       [newModel, indUSE] = DownExp(model, reacsInd, lowMax, indFvar, indRvar)
%
% INPUTS:
%    model:           model with TFA structure
%    reacsInd:        reaction indexes in the .rxns field
%    lowMax:          upper bound required for downregulated rxns
%
%
% OPTIONAL INPUTS:
%    indFvar:         indexes of forward fluxes (default = find them)
%    indRvar:         indexes of reverse fluxes (default = find them)
%
% OUTPUTS:
%    newModel:        model with constraints
%    indUSE:          indexes of use variables associated to downregulated
%                     reactions (tagged with 'DOWN_')
%
% Daniel F. Hernandez 2016
% Anush Chiappino-Pepe 2017 - redefinition of function
%

if ~exist('indFvar', 'var') || isempty(indFvar)
    [indFvar,~] = getAllVar(model,{'F'});
end

if ~exist('indRvar', 'var') || isempty(indRvar)
    [indRvar,~] = getAllVar(model,{'R'});
end

%reacsInd; lowMax;%THESE TWO MUST BE SAME SIZE AND THE F AND R VARSinds MUST
%BE MATCHED UP %NOTHING LESS THAN BASAL VALUES HERE!

% define use variables
intTag = {'DOWN'};
numVars = size(model.A,2);
indUSE = zeros(length(reacsInd),1);
for i=1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, '_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indUSE(i,1) = numVars+i;
end

% define contraints for MILP
numCons = size(model.A,1);
% F_rxn + R_rxn + 1000*DOWN_rxn < lowMax + 1000
% if DOWN_rxn = 1 => F_rxn + R_rxn < lowMax 
% if DOWN_rxn = 0 => F_rxn + R_rxn < lowMax + 1000
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'DOWN1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = lowMax(i) + 1000;
    model.A(numCons+i, indFvar(reacsInd(i))) = 1;
    model.A(numCons+i, indRvar(reacsInd(i))) = 1;
    model.A(numCons+i, indUSE(i)) = 1000;
    model.constraintType(numCons+i) = {'<'};
end

% [numCons,~] = size(model.A);
% % F_rxn + R_rxn + lowMax*DOWN_rxn < lowMax
% % if DOWN_rxn = 1 => F_rxn + R_rxn < 0
% % if DOWN_rxn = 0 => F_rxn + R_rxn < lowMax
% for i=1:length(reacsInd)
%     model.constraintNames(numCons+i) = strcat({'DOWN2_'}, model.rxns(reacsInd(i)));
%     model.rhs(numCons+i) = lowMax(i);
%     model.A(numCons+i, indFvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indRvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indUSE(i)) = lowMax(i);
%     model.constraintType(numCons+i) = {'>'};
% end

% Adder Var
[numCons, numVars] = size(model.A);
model.varNames(numVars+1,1) = {'all_Low'}; % name in UnifiedExpressionVar.m
model.var_ub(numVars+1,1) = Inf;
model.var_lb(numVars+1,1) = 0;
model.vartypes(numVars+1,1) = {'C'};
model.f(numVars+1,1) = 0;

% Adder Constraint
% sum(DOWN_rxn) = ALLUSEDOWNGE
model.constraintNames(numCons+1) = {'sum_allLow'};
model.rhs(numCons+1) = 0;
model.A(numCons+1, numVars+1) = -1;
model.A(numCons+1, indUSE) = 1;
model.constraintType(numCons+1) = {'='};

newModel = model;

