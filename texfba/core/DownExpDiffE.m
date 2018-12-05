function [newModel, indUSE] = DownExpDiffE(model, reacsInd, LowMax, indFvar, indRvar)
% Integrate constraints for downregulated genes bidirectional
%
% USAGE:
%
%       [newModel, indUSE] = DownExp(model, reacsInd, lowMax, indFvar, indRvar)
%
% INPUTS:
%    model:           model with TFA structure
%    reacsInd:        reaction indexes in the .rxns field
%    lowMax:          lower bound required for upregulated rxns (two
%                     columns for R_ and F_ vars)
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

%reacsInd; LowMax;%THESE TWO MUST BE SAME SIZE AND THE F AND R VARSinds MUST BE MATCHED UP %NOTHING LESS THAN BASAL VALUES HERE!

% define use variables
intTag = {'DOWN'};
numVars = size(model.A,2);
indUSE = zeros(length(reacsInd),1);
for i = 1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, '_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indUSE(i,1) = numVars+i;
end

% add auxiliary vars
numVars = numVars+length(reacsInd);
% numVars = size(model.A,2);
indUSEF = zeros(length(reacsInd),1);
for i = 1:length(reacsInd)
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
for i = 1:length(reacsInd)
    model.varNames(numVars+i,1) = strcat(intTag, 'r_', model.rxns(reacsInd(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indUSER(i,1) = numVars+i;
end

% define contraints for MILP
numCons = size(model.A,1);
% F_rxn + 1000*DOWNf_rxn < 1000 + LowMax(i,2)
% if DOWNf_rxn = 1 => F_rxn < LowMax(i,2)
% if DOWNf_rxn = 0 => F_rxn < 1000 + LowMax(i,2)
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'DOWNF1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 1000 + LowMax(i,2);
    model.A(numCons+i, indFvar(reacsInd(i))) = 1;
    model.A(numCons+i, indUSEF(i)) = 1000;
    model.constraintType(numCons+i) = {'<'};
end

numCons = numCons+length(reacsInd);
% numCons = size(model.A,1);
% R_rxn + 1000*DOWNr_rxn < 1000 + abs(LowMax(i,1))
% if DOWNr_rxn = 1 => R_rxn < abs(LowMax(i,1))
% if DOWNr_rxn = 0 => R_rxn < 1000 + abs(LowMax(i,1))
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'DOWNR1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 1000 + abs(LowMax(i,1));
    model.A(numCons+i, indRvar(reacsInd(i))) = 1;
    model.A(numCons+i, indUSER(i)) = 1000;
    model.constraintType(numCons+i) = {'<'};
end

numCons = numCons+length(reacsInd);
% numCons = size(model.A,1);
% DOWNr_rxn + DOWNf_rxn = DOWN_rxn
for i=1:length(reacsInd)
    model.constraintNames(numCons+i) = strcat({'DOWNFR1_'}, model.rxns(reacsInd(i)));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indUSER(i)) = -1;
    model.A(numCons+i, indUSEF(i)) = -1;
    model.A(numCons+i, indUSE(i)) = 1;
    model.constraintType(numCons+i) = {'='};
end

% % define contraints for MILP
% [numCons,~] = size(model.A);
% % F_rxn + 1000*DOWNf_rxn < LowMax(i,2)
% % if DOWNf_rxn = 1 => F_rxn < LowMax(i,2)
% % if DOWNf_rxn = 0 => F_rxn < LowMax(i,2)
% for i=1:length(reacsInd)
%     model.constraintNames{numCons+i} = strcat({'DOWNF2_'}, model.rxns(reacsInd(i)));
%     model.rhs(numCons+i) = LowMax(i,2);
%     model.A(numCons+i, indFvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indUSEF(i)) = LowMax(i,2);
%     model.constraintType(numCons+i) = {'>'};
% end
% 
% [numCons,~] = size(model.A);
% % R_rxn + abs(LowMax(i,1))*DOWNR_rxn > abs(LowMax(i,1))
% % if DOWNR_rxn = 1 => R_rxn < abs(LowMax(i,1))
% % if DOWNR_rxn = 0 => R_rxn < 1000 + abs(LowMax(i,1))
% for i=1:length(reacsInd)
%     model.constraintNames(numCons+i) = strcat({'DOWNR2_'}, model.rxns(reacsInd(i)));
%     model.rhs(numCons+i) = abs(LowMax(i,1));
%     model.A(numCons+i, indRvar(reacsInd(i))) = 1;
%     model.A(numCons+i, indUSER(i)) = abs(LowMax(i,1));
%     model.constraintType(numCons+i) = {'>'};
% end

% redefine sum_ALLUSEDOWNGE
% sum(DOWN_rxn) = sum_ALLUSEDOWNGE
if ismember({'sum_allLow'}, model.constraintNames)
    [~,a] = ismember({'sum_allLow'},model.constraintNames);
    model.A(a,indUSE) = 1;
else
    error('make sure the constraint sum_allLow is added in DownExp.m')
end

newModel = model;