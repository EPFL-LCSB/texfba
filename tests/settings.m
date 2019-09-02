% Initialization of paths and loading and preparing the model for the
% tex-fba test
%
% USAGE:
%
%    run('settings.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

clear
clc

% Define a directory where we want to save final results (.mat) 
% - if the directory does not exist it will be created automatically
% Note: run this script (do not copy and paste in the command window)
saving_directory = strrep(mfilename('fullpath'),...
    'tests/settings','tmpresults/');

% Initialize paths to tex-fba and dependencies (matTFA and
% cplex)
addpath(genpath(strrep(mfilename('fullpath'),...
    'tests/settings','texfba')));
[texfba_directory, thermo_data_directory] = initPhenoMappingPaths(...
    saving_directory);
