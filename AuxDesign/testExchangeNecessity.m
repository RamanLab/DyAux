function T = testExchangeNecessity(model, options, exchNames, receiver, mode, dfbaFun, minFold)
% Tests whether blocking an exchange stops growth of a species in the community.
% model     : mutant community model
% options   : same options used for the normal dFBA
% exchNames : cell array of exchange names, e.g. {'EX_leu__L_e','EX_ala__D_e'}
% receiver  : species index (1 or 2) that is expected to depend on the metabolite
% mode      : 'uptake' = block only uptake by the receiver
%             'both'   = block uptake and secretion in the receiver
%             'donor'  = block secretion by the other species (knockout of the donor's exchange)
% dfbaFun   : function handle, e.g. @dFBAComFSEOF
% minFold   : final/initial biomass below which a species counts as not growing (default 1.05)
if nargin < 7, minFold = 1.05; end
donor = 3 - receiver;                          % the other species
o0 = options; o0.Products = [];                % skip product tracking/FVA to save time
o0.blockUptakeRxns = []; o0.knockRxns = [];

% baseline: no block
base = dfbaFun(model, o0);
baseFold = base.biomassarr(:,end) ./ options.initBiomass(:);

nMet = numel(exchNames);
foldRecv  = nan(nMet,1);   % receiver final/initial biomass with block
foldDonor = nan(nMet,1);   % donor final/initial biomass with block
nBadSteps = nan(nMet,1);   % infeasible steps with block
found     = false(nMet,1); % exchange exists in both species

for n = 1:nMet
    % locate the exchange reaction index of the receiver and donor
    rxnRecv  = findExchange(model, exchNames{n}, receiver);
    rxnDonor = findExchange(model, exchNames{n}, donor);
    if isempty(rxnRecv) || isempty(rxnDonor), continue; end
    found(n) = true;

    o = o0;
    switch mode
        case 'uptake', o.blockUptakeRxns = rxnRecv;     % receiver cannot import
        case 'both',   o.knockRxns       = rxnRecv;     % receiver has no exchange at all
        case 'donor',  o.knockRxns       = rxnDonor;    % donor cannot secrete
    end

    r = dfbaFun(model, o);
    fold = r.biomassarr(:,end) ./ options.initBiomass(:); %fold change in receiver without blocking exchange of metabolite
    foldRecv(n)  = fold(receiver); %fold change in receiver after blocking uptake
    foldDonor(n) = fold(donor); %fold change in donor
    nBadSteps(n) = sum(r.solnstat ~= 1);
end

% receiver is "dead" if it did not grow, or the simulation became infeasible
stoppedGrowth = foldRecv < minFold | nBadSteps > 0;
T = table(exchNames(:), found, repmat(baseFold(receiver),nMet,1), foldRecv, foldDonor, ...
          'VariableNames', {'exchange','inBothSpecies','foldBaseline','foldReceiver', ...
                      'foldDonor'});
end

function rxnIdx = findExchange(model, exchName, sp)
% index of species sp's exchange reaction whose name ends with exchName
names = model.infoCom.EXsp(:,sp);
k = find(endsWith(names, exchName), 1);        % row in EXsp for this metabolite
if isempty(k) || model.indCom.EXsp(k,sp) == 0
    rxnIdx = [];
else
    rxnIdx = model.indCom.EXsp(k,sp);
end
end
