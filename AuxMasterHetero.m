tic
Soneidensis=readCbModel('iSO783_shewan.mat');Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

% define medium
options.mediumMets = {'glc__D_e[u]';'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'ser__L_e[u]';'so3_e[u]'};

%setup anaerobic conditions
options.mediumMets = setdiff(options.mediumMets, 'o2_e[u]','stable'); %remove oxygen

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;

%define organism list
speciesList = {Soneidensis, Kpneumoniae};
biomassList = {'Biomass','Biomass'};
spNameList = {'Soneidensis','Kpneumoniae'};

orgCount = length(speciesList);
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

% to set Vmax and Km and initial biomass
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;
options.initBiomass = [0.1,0.1];
options.delt = 0.1;
options.maxTime = 10;
options.solver = 'ibm_cplex';
options.carbonSource = 'glc__D_e[u]';

options.abdCutoff = 10;
options.Products={'EX_succ_e'};
Community = Communities{1,1};

% to set medium
Community = setMediumCom(Community, options.mediumMets,options.initMedium);

% If first organism is the donor then Donor=1, similarly Donor=2 if
% organismB is the donor
% To run mutual auxotroph, Donor=0
% Size denotes number of results required for Mutual auxotrophs

[sortedAuxotrophScoreSingle1,sortedAuxotrophResultsSingle1,runInfo1] = FindAuxotroph(Community,pairedModelInfo,options,1,5);
[sortedAuxotrophScoreMutual,sortedAuxotrophResultsMutual,runInfo] = FindAuxotroph(Community,pairedModelInfo,options,0,10);
[sortedAuxotrophScoreSingle2,sortedAuxotrophResultsSingle2,runInfo2] = FindAuxotroph(Community,pairedModelInfo,options,2,5);
toc
