tic
Soneidensis=readCbModel('iSO783_shewan.mat');Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

% define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};

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

[sortedAuxotrophScoreSingle1,sortedAuxotrophResultsSingle1] = FindAuxotroph(Community,pairedModelInfo,options,1,3);
[sortedAuxotrophScoreMutual,sortedAuxotrophResultsMutual] = FindAuxotroph(Community,pairedModelInfo,options,0,3);
[sortedAuxotrophScoreSingle2,sortedAuxotrophResultsSingle2] = FindAuxotroph(Community,pairedModelInfo,options,2,3);
toc
