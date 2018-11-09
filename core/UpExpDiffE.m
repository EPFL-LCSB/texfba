function [newModel, indUSE] = UpExpDiffE(model, reacsInd, UpMin, indFvar, indRvar)
% Integrate constraints for upregulated genes bidirectional
%
% USAGE:
%
%       [newModel, indUSE] = UpExp(model, reacsInd, UpMin, indFvar, indRvar)
%
% INPUTS:
%    model:           model with TFA structure
%    reacsInd:        reaction indexes in the .rxns field
%    UpMin:           lower bound required for upregulated rxns (two
%                     columns for R_ and F_ vars)
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
% Anush Chiappino-Pepe 2017 - redefinition of function
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
    model.f(numVars+i,1) = 0;
    indUSE(i,1) = numVars+i;
end

% add auxiliary vars
numVars = numVars+length(reacsInd); % sometimes the function size does not report the updated value! this is a temporary solution
% numVars = size(model.A,2);
indUSEF = zeros(length(reacsInd),1);
for i=1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, 'f_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indUSEF(i,1) = numVars+i;
end

numVars = numVars+length(reacsInd);
% numVars = size(model.A,2);
indUSER = zeros(length(reacsInd),1);
for i=1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, 'r_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indUSER(i,1) = numVars+i;
end

% define contraints for MILP
numCons = size(model.A,1);
% - F_rxn + abs(UpMin(i,2))*UPf_rxn < 0
% if UPf_rxn = 1 => F_rxn > abs(LowMax(i,2))
% if UPf_rxn = 0 => F_rxn > 0
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'UPF1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indFvar(reacsInd(i))) = -1;
    model.A(numCons+i, indUSEF(i)) = UpMin(i,2);
    model.constraintType(numCons+i) = {'<'};
end

numCons = numCons+length(reacsInd);
% numCons = size(model.A,1);
% - R_rxn + abs(UpMin(i,1))*UPr_rxn < 0
% if UPr_rxn = 1 => R_rxn > abs(UpMin(i,1))
% if UPr_rxn = 0 => R_rxn > 0
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'UPR1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indRvar(reacsInd(i))) = -1;
    model.A(numCons+i, indUSER(i)) = abs(UpMin(i,1));
    model.constraintType(numCons+i) = {'<'};
end

numCons = numCons+length(reacsInd);
% numCons = size(model.A,1);
% UPr_rxn + UPf_rxn = UP_rxn
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'UPFR1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indUSER(i)) = -1;
    model.A(numCons+i, indUSEF(i)) = -1;
    model.A(numCons+i, indUSE(i)) = 1;
    model.constraintType(numCons+i) = {'='};
end

% % define contraints for MILP
% [numCons,~] = size(model.A);
% % F_rxn - 1000*UPf_rxn < UpMin(i,2)
% % if UPf_rxn = 1 => F_rxn < UpMin(i,2) + 1000
% % if UPf_rxn = 0 => F_rxn < UpMin(i,2)
% for i=1:length(reacsInd)
%     model.constraintNames(numCons+i) = strcat({'UPF2_'}, model.rxns(reacsInd(i)));
%     model.rhs(numCons+i) = UpMin(i,2);
%     model.A(numCons+i, indFvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indUSEF(i)) = -1000;
%     model.constraintType(numCons+i) = {'<'};
% end
% 
% [numCons,~] = size(model.A);
% % R_rxn - 1000*UPr_rxn < abs(UpMin(i,1))
% % if UPr_rxn = 1 => R_rxn < abs(LowMax(i,1)) + 1000
% % if UPr_rxn = 0 => R_rxn < abs(LowMax(i,1))
% for i=1:length(reacsInd)
%     model.constraintNames(numCons+i) = strcat({'UPR2_'}, model.rxns(reacsInd(i)));
%     model.rhs(numCons+i) = abs(UpMin(reacsInd(i),1));
%     model.A(numCons+i, indRvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indUSER(i)) = -1000;
%     model.constraintType(numCons+i) = {'<'};
% end

% redefine sum_ALLUSEUPGE
% sum(UP_rxn) = sum_ALLUSEUPGE
if ismember({'sum_allHigh'}, model.constraintNames)
    [~,a] = ismember({'sum_allHigh'},model.constraintNames);
    model.A(a,indUSE) = 1;
else
    error('make sure the constraint sum_allHigh is added in UpExp.m')
end

newModel = model;