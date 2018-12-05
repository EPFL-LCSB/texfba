function [expModel, solution, newgrRules, expModelm] = integrateGeneExp(tmodel, ...
    levelGenes, mingrRate, lowPvalue, highPvalue, percLow, percHigh, ...
    numAlt, selectAlt, flagPlot, minmax, geneKO, rxnLevelOUT)
% Integrate gene expression data in the model
%
% USAGE:
%
%       [expModel, solution, newgrRules, expModelm] = integrateGeneExp(tmodel, mingrRate, lowPvalue, highPvalue, percHigh, percLow, numAlt, selectAlt, flagPlot, minmax, geneKO)
%
% INPUTS:
%    tmodel:          model with TFA structure
%    levelGenes:      quantitative gene expression value for all genes in
%                     the model. If lacking data, it is recomended to
%                     assign an average expression value for that genes.
%
% OPTIONAL INPUTS:
%    mingrRate:       min growth to achieve (default = 0.01)
%    lowPvalue:       percentile to assign lowly expressend genes (default
%                     = 20)
%    highPvalue:      percentile to assign highly expressend genes (default
%                     = 80)
%    percLow:         percentage from  (default
%                     = 80)
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    selectAlt:       Indicate what alternative expression profile to
%                     integrate in the model (default = common profile to 
%                     max consistency score (CS))
%    flagPlot:        plot distribution of gene expression (default = 1)
%    minmax:          minmax of net reaction fluxes in two columns 
%                     (default = lb and ub)
%    geneKO:          gene IDs to KO
%
% OUTPUTS:
%    expModel:        model with TFA structure, gene expression 
%                     constraints and information integrated (check field
%                     ExpInfo)
%    solution:        solution structure with information on CS,
%                     alternatives, etc
%    newgrRules:      new GPR rules modified from the original GPR (low 
%                     expressed genes are excluded from complex GPR)
%    expModelm:       model with TFA structure and gene expression 
%                     constraints without a profile integrated
% 
%
% .. Author:
% Daniel F. Hernandez & Vikash Pandey 2015
% Anush Chiappino-Pepe 2017 - integration of expression constraints into 
% the model and organization of the function
%

if (nargin < 3)
    mingrRate = 0.1;
end
if (nargin < 4)
    lowPvalue = 20;
end
if (nargin < 5)
    highPvalue = 80;
end
if (nargin < 6)
    percLow = 2E-5;
end
if (nargin < 7)
    percHigh = 2E-3;
end
if (nargin < 8)
    numAlt = 5;
end
if (nargin < 9)
    selectAlt = 0;
end
if (nargin < 10)
    flagPlot = 1;
end
if (nargin < 11)
    minmax = [];
end
if (nargin < 12)
    geneKO = [];
end
if (nargin < 13)
    rxnLevelOUT = {};
end

if isempty(minmax)
    minmax=[tmodel.lb tmodel.ub];
end

path_save = '/Users/Anush/Documents/gei_altsol.mat';

% Prepare model for gene expression data integration
if flagPlot
    % change gene to rxn expression
    levels = geneToreaction_levels(tmodel, tmodel.genes, levelGenes, @min, @max );
    lowPer1 = prctile(log(levels),lowPvalue);
    higPer1 = prctile(log(levels),highPvalue);
    % plot distribution
    hist(log(levels),10)
    hold on
    plot([lowPer1; lowPer1], [300*ones(1,length(lowPer1)); ...
        zeros(1,length(lowPer1))])
    plot([higPer1; higPer1], [300*ones(1,length(higPer1)); ...
        zeros(1,length(higPer1))])
    
    ylabel('Number of genes','FontWeight','bold','FontSize',12);
    xlabel('LOG(expression level)','FontWeight','bold','FontSize',12);
    set(gca,'FontSize',10,'FontName','Arial','FontWeight','normal')
    set(gcf, 'Units', 'centimeters','Position', [15 15 10 7.5])
end

geneRatio = ones(numel(levelGenes),1);
y = prctile(levelGenes,[lowPvalue highPvalue]);
geneRatio(levelGenes>y(2)) = 2;
geneRatio(levelGenes<y(1)) = 0;

% calculate rxn expression (2 1 0: high medium low)
[rxns, rxnExp, newgrRules, optionOut] = fEvalGPR(tmodel, tmodel.genes, geneRatio, @min, @max);
tmodel.genesnewgrRules = optionOut.genes;
tmodel.rxnsnewgrRules = rxns;
tmodel.newgrRules = newgrRules;

% add expression constraint
% find high and low expressed reaction index
[~, highInd] = ismember(rxns(rxnExp==2),tmodel.rxns);
[~, lowInd] = ismember(rxns(rxnExp==0),tmodel.rxns);

if ~isempty(rxnLevelOUT)
    [y,r] = ismember(rxnLevelOUT,tmodel.rxns);
    highInd = highInd(~ismember(highInd,r));
    lowInd = lowInd(~ismember(lowInd,r));
    if ~all(y)
        warning('rxnIDs in rxnLevelOUT are not in the model!')
    end
end

% F, R and NF vars
indRvar = getAllVar(tmodel,{'R'});
indFvar = getAllVar(tmodel,{'F'});
indNFvar = getAllVar(tmodel,{'NF'});
if isempty(indNFvar)
    model = addNetFluxVariables(model);
    indNFvar = getAllVar(model,{'NF'});
end

expModelm = addExpCnsts(tmodel, highInd, lowInd, minmax, percHigh,...
    percLow, indFvar, indRvar);

% calculate the consistency score
expModelm.var_lb(ismember(expModelm.varNames, ...
    tmodel.varNames(tmodel.f==1))) = mingrRate;
alt = optimizeThermoModel(expModelm);
if isempty(alt.val) || isnan(alt.val)
    error('model is not feasible with defined mingrRate')
end

indUPDOWN = [expModelm.ExpInfo.indUP; expModelm.ExpInfo.indDOWN];
if isempty(indUPDOWN)
    indUPDOWN = [];
    listofvars = {'UP','DOWN'};
    for i = 1:numel(listofvars)
        ind = getAllVar(model,listofvars(i));
        indUPDOWN = [indUPDOWN; ind];
    end
end
if ~isempty(geneKO)
    expModelm = thermoDeleteModelGenes(expModelm,geneKO);
end 
solution = getExpProfile(expModelm, numAlt, indUPDOWN, indNFvar, path_save);
% if pipeline broken upload: solution=load(path_save);

% define objective function back to original model
expModelm.f = zeros(numel(expModelm.f),1);
expModelm.f(ismember(expModelm.varNames, tmodel.varNames(tmodel.f==1))) = 1;
expModelm.var_lb(ismember(expModelm.varNames, tmodel.varNames(tmodel.f==1))) = 0;

% integrate expression constraints into model
if isequal(selectAlt,0) % define maximum consistency;
    expModel = integrateConsProfile(expModelm, solution, find(ismember({'Total_UpsAndDowns'},expModelm.varNames)), selectAlt);
else
    expModel = integrateConsProfile(expModelm, solution, indUPDOWN, selectAlt);
end