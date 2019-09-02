function sol = optimizeThermoModel(tModel,tagParTFA,solver,time)
% This function solves a TFA problem using CPLEX
%
% USAGE:
%
%    sol = optimizeThermoModel(tModel,tagParTFA,solver,time)
%
% INPUT:
%    tModel:          Model with TFA structure
%
% OPTIONAL INPUTS:
%    tagParTFA:       Tag to perform parsimonious TFA (default = 0)
%    solver:          Solver to use (default = 'cplex'). Other solvers have
%                     not been implemented yet
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%                     
% OUTPUTS:
%    sol:             Solution structure from the solver. The optimal value
%                     is stored in sol.val. The sol.stat values are defined
%                     here based on the solution from the solver. If there 
%                     is no solution sol.stat will be zero, else it is one.
%
% .. Author:
% Anush Chiappino-Pepe 2016
%

if (nargin < 2)
    tagParTFA = 0;
end
if (nargin < 3)
    solver = 'cplex'; % to be updated
end
if (nargin < 4)
    time = [];
end

sol = solveTFAmodelCplex(tModel,time);

if tagParTFA
    model4tagParTFA = tModel;
    model4tagParTFA.var_lb(model4tagParTFA.f==1) = sol.val;
    FRidx = getAllVar(model4tagParTFA,{'F','R'});
    f = zeros(numel(model4tagParTFA.f),1);
    f(FRidx) = 1;
    model4tagParTFA.f = f;
    model4tagParTFA.objtype = 1;
    sol = solveTFAmodelCplex(model4tagParTFA);
end

if ~isempty(sol.val) && ~isnan(sol.val)
    sol.stat = 1;
else
    sol.stat = 0;
    sol.val = NaN;
end