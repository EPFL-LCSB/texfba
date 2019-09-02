function newModel = UnifiedExpVar(model)
% Define consistency score as objective function
%
% USAGE:
%
%       newModel = UnifiedExpVar(model)
%
% INPUTS:
%    model:           model with TFA structure after integration of gene 
%                     expression constraints with UpExp.m / DownExp.m
%
%
% OUTPUTS:
%    newModel:        model with consistency score as objective function
%
% Daniel F. Hernandez 2016
% Anush Chiappino-Pepe 2017
%

[numCons, numVars] = size(model.A);
model.varNames{numVars+1,1} = 'Total_UpsAndDowns';
model.var_ub(numVars+1,1) = Inf;
model.var_lb(numVars+1,1) = 0;
model.vartypes{numVars+1,1} = 'C';
model.f(numVars+1,1) = 1;

[~,a] = ismember({'all_High'},model.varNames);
[~,b] = ismember({'all_Low'},model.varNames);

% sum(UP_rxn) + sum(DOWN_rxn) = TotalUpsAndDowns
model.constraintNames{numCons+1} = 'SumExpCons';
model.rhs(numCons+1) = 0;
model.A(numCons+1, numVars+1) = -1;
model.A(numCons+1, [a,b]) = 1;
model.constraintType{numCons+1} = '=';

newModel=model;