function model = integrateConsProfile(model, solution, ind, selectAlt)
% Integrate an expression profile in the model after getExpProfile.m
%
% USAGE:
%
%       model = integrateConsProfile(model, solution, ind, selectAlt)
%
% INPUTS:
%    model:           model with TFA structure
%    solution:        solution structure from getExpProfile.m
%
% OPTIONAL INPUTS:
%    ind:             indexes of use variables to fix (default = 'UP' and
%                     'DOWN')
%    selectAlt:       number of alternative solution to integrate (default
%                     = 0, which just required max CS)
%
% OUTPUTS:
%    model:           model with TFA structure and a gene expression
%                     profile integrated
%
%
% .. Author:
% Anush Chiappino-Pepe 2017

if (nargin < 4)
    selectAlt = 0;
end

if isequal(selectAlt,0) % define maximum consistency;
    [~,ind] = ismember({'Total_UpsAndDowns'}, model.varNames);
    model.var_lb(ind) = solution.store_obj(1);
else
    model.var_lb(ind(solution.z_matrix(:,selectAlt)==1)) = 1;
    model.var_ub(ind(solution.z_matrix(:,selectAlt)==0)) = 0;
end
end