%% Setting up paths for TEXFBA
% files to help set up the paths for TEX-FBA have been provided in 
% the texfba/tests folder

% prior to running phenomapping you should install matTFA.

% it is recommended that you remove any other file from the matlab path 
% prior to starting phenomapping

% OPTION A: add the required repositories and solver directories to the 
% matlab path in an automatic fashion, meaning that texfba will 
% automatically search for matTFA in the parent folder and if 
% this is not there you will be asked to provide the path. 
% For this option run the "settings.m" script

% OPTION B: add manually the required repositories and solver paths to the 
% matlab path. For this option you can follow the backbone of the code 
% presented below. 
% Please, note that since option B is manual, you should adapt the paths

clear
clc
close all

% adding phenomapping and dependent repositories to the path
cplex_directory = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271'; % provide path to cplex
mattfa_directory = '/Users/Anush/GIT_Folders/matTFA'; % provide path to matTFA repository
texfba_directory = '/Users/Anush/GIT_Folders/texfba'; % provide path to TEXFBA repository

addpath(genpath(mattfa_directory));
addpath(genpath(texfba_directory));
cd(texfba_directory) % work from first folder of texfba as the starting path
