function [texfba_directory, thermo_data_directory] = ...
    initTEXFBAPaths(saving_directory, mattfa_directory)
% Initialization of paths for texfba
%
% USAGE:
%
%    [texfba_directory, thermo_data_directory] = initTEXFBAPaths(saving_directory, mattfa_directory)
%
% OPTIONAL INPUTS:
%    saving_directory:      Directory where we want to save intermediate  
%                           and final results (default = tmpresults folder  
%                           in parent texfba folder)
%    mattfa_directory:      Directory to the matTFA repository (default = 
%                           empty / provide manually* see note)
%
% OUTPUTS:
%    texfba_directory:Directory to the parent texfba folder
%    thermo_data_directory: Directory to the thermodynamic data within
%                           matTFA
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

if (nargin < 1)
    saving_directory = strrep(mfilename('fullpath'),...
    'texfba/inittexfbaPaths','tmpresults/');
end
texfba_directory = strrep(mfilename('fullpath'),...
    'texfba/initTEXFBAPaths','');
addpath(genpath(texfba_directory));

% NOTE: it seems like uigetdir has a bug in the latest versions of 
% matlab... this is a temporal solution - make sure you indeed install 
% matTFA in the same paths as texfba. When this bug is
% resolved the directories could be empty by default

if (nargin < 2)
    mattfa_directory = strrep(texfba_directory,'texfba/',...
        'matTFA');
end


if ~isdir(mattfa_directory) || isempty(mattfa_directory)
    mattfa_directory = uigetdir('matTFA','select the matTFA directory');
end
addpath(genpath(mattfa_directory));
thermo_data_directory = strcat(mattfa_directory,...
    '/thermoDatabases/thermo_data.mat');

cd(texfba_directory)

% texfba was developed to work with the solver CPLEX. We hence check 
% that you have CPLEX installed and on the path. In future releases, 
% this repository will work with other solvers like gurobi.
cplex_directory = what('cplex');
if isempty(cplex_directory)
    cplex_directory = [];
else
    cplex_directory = cplex_directory.path;
end
addCplexPath(cplex_directory);

% create temporary results folder to save intermediate and final results: 
% note that texfba will not upload the intermediate results 
% - but you might need to use them if matlab crashes
if ~isdir(saving_directory(1:end-1))
    mkdir(saving_directory(1:end-1))
end