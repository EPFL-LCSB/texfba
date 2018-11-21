clear
clc
close all

addpath(genpath('/Users/Anush/GIT_Folders/fba_toolbox'));
addpath(genpath('/Users/Anush/GIT_Folders/et'));
addpath(genpath('/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271'));
addpath(genpath('/Users/Anush/GIT_Folders/texfba'));
cd('/Users/Anush/GIT_Folders/texfba')

%% INPUTS
mm = load('/Users/Anush/SWITCHdrive/texFBA/texfba/Data/models/FBAmodel_withC13.mat');
model{1} = mm.fba_model;

mm = load('/Users/Anush/SWITCHdrive/texFBA/texfba/Data/models/xFBAmodel.mat');
model{2} = mm.xFBA;

mm = load('/Users/Anush/SWITCHdrive/texFBA/texfba/Data/models/TFBAmodel_fluxomics_metabolomics_ishii.mat');
model{3} = mm.fba_tmodel_conc;

mm = load('/Users/Anush/SWITCHdrive/texFBA/texfba/Data/models/xTFAmodel.mat');
model{4} = mm.xTFA;

essThr = 0.1;
solt = cell(1,4);
indNF = cell(1,4);
essential_rxns_tfa = cell(1,4);
yesrxntfa = cell(1,4);
essential_genes_tfa = cell(1,4);
yesgenetfa = cell(1,4);

%% 
for i = 2:4
    solt{i} = optimizeThermoModel(model{i});
    if solt{i}.val > 1.5
        solution.store_obj = solt{i}.val;
        % define objective function back to original model
        model{i}.f = zeros(numel(model{i}.f),1);
        model{i}.f(ismember(model{i}.varNames, strcat('NF_',model{i}.rxns(model{i}.c==1)))) = 1;
        model{i}.var_lb(model{i}.f==1) = 0;
        indUPDOWN = [model{i}.ExpInfo.indUP; model{i}.ExpInfo.indDOWN];
        model{i} = integrateConsProfile(model{i}, solution, indUPDOWN, 0);
        solt{i} = optimizeThermoModel(model{i});
    end
    indNF{i} = getAllVar(model{i},{'NF'});
    if not(isempty(indNF{i}))
        [grRatio_rxntfa, grRateKO_rxntfa] = thermoSingleRxnDeletion(model{i}, 'TFA', model{i}.rxns, 0, 0, essThr, indNF{i});
        grRateKO_rxntfa(isnan(grRateKO_rxntfa)) = 0;
        essential_rxns_tfa{i} = model{i}.rxns(grRateKO_rxntfa(:,1)<essThr*solt{i}.val);
        yesrxntfa{i} = ismember(model{i}.rxns,essential_rxns_tfa{i});
        
        [grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(model{i}, 'TFA', model{i}.genes, 0, 0, 0, essThr, indNF{i});
        grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0;
        essential_genes_tfa{i} = model{i}.genes(grRateKO_genetfa(:,1)<essThr*solt{i}.val);
        yesgenetfa{i} = ismember(model{i}.genes,essential_genes_tfa{i});
    end
end

