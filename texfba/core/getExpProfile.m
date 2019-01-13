function solution = getExpProfile(model, numAlt, indUPDOWN, indNFvar, pathSave)
% Obtain alternative solutions for gene expression integration problem
%
% USAGE:
%
%       solution = getExpProfile(model, numAlt, listofvars, indNFvar, pathSave)
%
% INPUTS:
%    model:           model with TFA structure and gene expression
%                     constraints (output of addExpCnsts.m)
%
%
% OPTIONAL INPUTS:
%    numAlt:          number of alternatives wanted (default = 1)
%    listofvars:      cell defined as {'UP','DOWN'}, or any of them
%                     (default = {'UP','DOWN'})
%    indNFvar:        indexes of net fluxes (default = empty - to find them)
%    pathSave:        path to save all the results (default = '/Documents')
%
% OUTPUTS:
%    solution:        solution field with:
%       v_matrix      net flux solution
%       z_matrix      which combination of UP and DOWN was applied (1=on)
%       store_obj     consistency score
%       sol_matrix    all solutions saved in a matrix
%       model         model with integer cuts included to generate more sols
%
% Vikash Pandey 2015
% Anush Chiappino-Pepe 2017 - redefinition of function

if (nargin < 2)
    numAlt = 1;
end
if (nargin < 3)
    indUPDOWN = [];
end
if (nargin < 5)
    indNFvar = [];
end
if (nargin < 5)
    pathSave = '/Documents';
end

if isempty(indNFvar)
    indNFvar = getAllVar(model,{'NF'});
end

if isempty(indUPDOWN)
    indUPDOWN = [];
    listofvars = {'UP','DOWN'};
    for i = 1:numel(listofvars)
        ind = getAllVar(model,listofvars(i));
        indUPDOWN = [indUPDOWN; ind];
    end
end

sol=optimizeThermoModel(model);

if numAlt>100
    model.var_lb(model.f==1)=sol.val-0.5;
end

if isempty(sol.x)
    v_matrix = [];
    z_matrix = [];
    sol_matrix = [];
    store_obj = [];
else
    v_matrix = zeros(length(sol.x(indNFvar)), numAlt);
    z_matrix = zeros(length(sol.x(indUPDOWN)), numAlt);
    sol_matrix = zeros(length(sol.x), numAlt);
    store_obj = zeros(1, numAlt);
end

NumSols = 0;

% find alternate solutions
fprintf('getting alternative solutions\n');
while ((NumSols < numAlt) && ~(isempty(sol.x)))
    
    [numCons,numVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        fprintf('Number of alternative:\t%d\n',NumSols);
        v_matrix(:,NumSols) = sol.x(indNFvar);
        z_matrix(:,NumSols) = sol.x(indUPDOWN);
        sol_matrix(:,NumSols) = sol.x;
        store_obj(NumSols) = sol.val;
        % we find all the use vectors and formulate them into a new integer cut
        % constraint
        USEvecSol = ones(numVars,1);
        USEvecSol(indUPDOWN) = sol.x(indUPDOWN);
        actUSEvec = find(USEvecSol<0.1);
        
        NewCons = zeros(1,numVars);
        NewCons(actUSEvec) = 1;
        
        model.A(numCons+1,:) = NewCons;
        model.rhs(numCons+1) = 0.5;
        model.constraintNames{numCons+1} = ['CUT_' num2str(NumSols)];
        model.constraintType{numCons+1} = '>';
        
        sol = optimizeThermoModel(model);

        if isempty(sol.x)
            break
        end
        if rem(numAlt, 50) == 0
            save(pathSave,'v_matrix','z_matrix','store_obj','sol_matrix');
        end
    end
end

solution.v_matrix = v_matrix(:,1:NumSols);
solution.z_matrix = z_matrix(:,1:NumSols);
solution.sol_matrix = sol_matrix(:,1:NumSols);
solution.store_obj = store_obj(1:NumSols);
solution.model = model;
solution.maxThrCS = length(model.ExpInfo.indUP) + length(model.ExpInfo.indDOWN);

end