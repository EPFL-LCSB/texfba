function newModel = addExpCnsts(model, UpReacs1, DownReacs1, mm, ...
percHigh, percLow, indFvar, indRvar, visualize)
% Integrate gene expression constraints
%
% USAGE:
%
%       newModel = addExpCnsts(model, UpReacs1, DownReacs1, mm, percHigh, percLow, indFvar, indRvar, visualize)
%
% INPUTS:
%    model:           model with TFA structure
%    UpReacs1:        indexes of highly expressed reactions. Indexes are
%                     equivalent to the position of the rxn in the rxns
%                     field
%    DownReacs1:      indexes of lowly expressed reactions. Indexes are
%                     equivalent to the position of the rxn in the rxns
%                     field
%    mm:              FVA of the model in a two column matrix
%    percHigh:        percentage that constraints will be set at for highly
%                     expressed reactions
%    percLow:         percentage that constraints will be set at for lowly
%                     expressed reactions
%
%
% OPTIONAL INPUTS:
%    indFvar:         indexes of F rxn variables (default = find them)
%    indRvar:         indexes of R rxn variables (default = find them)
%    visualize:       visualize constraints on rxn ranges (default = 0)
%
% OUTPUTS:
%    newModel:        model with gene expression constraints and field
%                     ".ExpInfo" summarizing information and constraints
%                     integrated
%
%
% .. Author:
% Daniel F. Hernandez 2016
% Anush Chiappino-Pepe 2017 - restructuring

% The constraints are explained below. Note that
% all expression variables (other than the adders) start with UP_ or DOWN_
% and are binary. They act directly on the F_ and R_ variables.
% bidirectional reactions have helper vars UPf_ UPr_ DOWNf_ DOWNr_

% Up reactions
% Unidirectional Up reaction:
% general:
% -F_rxn - R_rxn + UpMin*UP_rxn < 0
% guarantees the following:
% if UP_rxn = 1 => F or R > UpMin

% Bidirectional Up reaction:
% general:
% UPr_rxn + UPf_rxn = UP_rxn
% -F_rxn + UpMin*UPf_rxn < 0
% guarantees the following:
% if UPf_rxn = 1 => F > UpMin
% R_rxn + UpMin*UPr_rxn < 0
% guarantees the following:
% if UPr_rxn = 1 => R > UpMin

% Down reaction
% Unidirectional Down reaction:
% general:
% F_rxn + R_rxn + 1000*UP_rxn < 1000 + lowMax
% guarantees the following:
% if UP_rxn = 1 => F or R > lowMax

% Bidirectional Down reaction:
% general:
% DOWNr_rxn + DOWNf_rxn = DOWN_rxn
% F_rxn + 1000*DOWNf_rxn < 1000 + lowMax
% guarantees the following:
% if DOWNf_rxn = 1 => F < lowMax
% R_rxn + 1000*DOWNr_rxn < 1000 + lowMax
% guarantees the following:
% if DOWNr_rxn = 1 => R < lowMax

% All Expression adders(named as "Total_UpsAndDowns" should be your last variable):
%  SumExpCons: sum(UP_rxn) + sum(DOWN_rxn) - Total_UpsAndDowns = 0


if ~exist('visualize', 'var') || isempty(visualize)
    visualize = 0;
end
if ~exist('indRvar', 'var') || isempty(indRvar)
    indRvar = getAllVar(model,{'R'});
end
if ~exist('indFvar', 'var') || isempty(indFvar)
    indFvar = getAllVar(model,{'F'});
end

% identify transports and remove them from analysis
tptRxns = checkTransport(model);
if ~isfield(model,'isTrans')
    model.isTrans=tptRxns.isTrans;
end
UpReacs2 = setdiff(UpReacs1,find(model.isTrans));
DownReacs2 = setdiff(DownReacs1,find(model.isTrans));

%removed blocked and low-variability reacs
UpReacs   = UpReacs2(diff(mm(UpReacs2,:),[],2)>.1);
DownReacs = DownReacs2(diff(mm(DownReacs2,:),[],2)>.1);

% choose flux thresholds based on minmax
% High
% ranges that dont cross 0: unidirectional, only one constraint
posflux_high = mm(UpReacs,1)>=0; % positive flux
a = mm(UpReacs(posflux_high),1);
b = mm(UpReacs(posflux_high),2);
highPerPos = percHigh*b+(1-percHigh)*a;

negflux_high = mm(UpReacs,2)<=0; % negative flux
a = mm(UpReacs(negflux_high),1) ;
b = mm(UpReacs(negflux_high),2);
highPerNeg = (1-percHigh)*b+percHigh*a;

saverHigh = nan(length(UpReacs),2);
saverHigh(posflux_high,:) = [highPerPos highPerPos]; % only the lb will be used in UpExp
saverHigh(negflux_high,:) = [highPerNeg highPerNeg]; % only the lb will be used in UpExp

% ranges that cross 0: will have to add two constraints
crossflux_high = ~posflux_high & ~negflux_high;
highPerCross = mm(UpReacs(crossflux_high),:)*percHigh;
saverHigh(crossflux_high,:) = highPerCross;

% Low
% ranges that dont cross 0: unidirectional, only one constraint
posflux_low = mm(DownReacs,1)>=0 ;% positive flux
a_low = mm(DownReacs(posflux_low),1) ;
b_low = mm(DownReacs(posflux_low),2);
LowPerPos = percLow*b_low+(1-percLow)*a_low;

negflux_low = mm(DownReacs,2)<=0 ;% negative flux
a_low = mm(DownReacs(negflux_low),1) ;
b_low = mm(DownReacs(negflux_low),2);
LowPerNeg = (1-percLow)*b_low+percLow*a_low;

saverLow = nan(length(DownReacs),2);
saverLow(posflux_low,:) = [LowPerPos LowPerPos];
saverLow(negflux_low,:) = [LowPerNeg LowPerNeg];

% ranges that cross 0 : will have to add two constraints
crossflux_low = ~posflux_low & ~negflux_low;
LowPerCross = mm(DownReacs(crossflux_low),:)*percLow;
saverLow(crossflux_low,:) = LowPerCross;

% plot to visualize
if visualize
    mm4 = mm;
    [~,I] = sort(max(abs(mm4(UpReacs,:)'))');
    herrorbar(zeros(length(UpReacs),1),1:length(UpReacs) , mm4(UpReacs(I),1),...
        mm4(UpReacs(I),2)); clear mm4
    hold on
    scatter([saverHigh(I,1); saverHigh(I,2)],[(1:length(UpReacs))'; (1:length(UpReacs))'])
    
    mm4 = mm;
    [~,I] = sort(max(abs(mm4(DownReacs,:)'))');
    figure;	herrorbar(zeros(length(DownReacs),1),1:length(DownReacs) , mm4(DownReacs(I),1),...
        mm4(DownReacs(I),2)); clear mm4
    hold on
    scatter([saverLow(I,1); saverLow(I,2)],[(1:length(DownReacs))'; (1:length(DownReacs))'])
end

%  add expression constraints
model.f = zeros(size(model.f));

%take transcription data from Gem and map to small model
fprintf('defining MILP problem\n');
[model_UP1, indUP1] = UpExp(model,UpReacs(posflux_high | negflux_high), abs(saverHigh((posflux_high | negflux_high),1)), indFvar, indRvar);
[model_UP, indUP]  = UpExpDiffE(model_UP1, UpReacs(crossflux_high), saverHigh(crossflux_high,:), indFvar, indRvar);
indUP = [indUP1; indUP];

[mBoth1, indDOWN1] = DownExp(model_UP,DownReacs(posflux_low | negflux_low),abs(saverLow((posflux_low | negflux_low),1)), indFvar, indRvar);
[mBoth, indDOWN]  = DownExpDiffE(mBoth1, DownReacs(crossflux_low), saverLow(crossflux_low,:), indFvar, indRvar);
indDOWN = [indDOWN1; indDOWN];

%use this only if we want to fix a total number of the binary expression variables
newModel = UnifiedExpVar(mBoth);

% save info of gene expression constraints
newModel.ExpInfo.mm = mm;
newModel.ExpInfo.percHigh = percHigh;
newModel.ExpInfo.percLow = percLow;
newModel.ExpInfo.UpReacsFinal = UpReacs;
newModel.ExpInfo.DownReacsFinal = DownReacs;
newModel.ExpInfo.UpReacsInput = UpReacs1;
newModel.ExpInfo.DownReacsInput = DownReacs1;
newModel.ExpInfo.indUP = indUP;
newModel.ExpInfo.indDOWN = indDOWN;

