tic
Soneidensis=readCbModel('iSO783_shewan.mat'); Kpneumoniae = readCbModel('iYL1228_klebsi.mat');
%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;

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
options.ScoreWeights = [1/3, 1/3, 1/3];

options.abdCutoff = 10;
options.Products={'EX_succ_e'};
Community = Communities{1,1};

dels = {'Soneidensis_ALAR','Kpneumoniae_ILETA'};
%to create auxotrophic pair 
[Community] = changeRxnBounds(Community,dels,0,'b');

%to set medium
for i = 1: length(Communities)
    Community = setMediumCom(Community, options.mediumMets,options.initMedium);
end

%Dynamic FBA for all communities
allResults = {};
for i= 1 : size(Communities,2)
    disp (i);
    allResults{i} = dFBAComFSEOF(Community,options);
    CommunityResultsWT{1,i} = compileComFSEOF(Community,pairedModelInfo,options,allResults{i});
end

%perform FSEOF
[FseofAll] = FSEOFall_comm(Community,pairedModelInfo(i,:),options);

j=1;
for i=1:length(FseofAll{1,3})
    if contains(FseofAll{1,3}{j,1},'tex')||contains(FseofAll{1,3}{j,1},'tpp')
        FseofAll{1,3}(j,:)=[];
        j=j-1;
    end
    j=j+1;
end

j=1;
for i=1:length(FseofAll{1,2})
    if contains(FseofAll{1,2}{j,1},'tex') || contains(FseofAll{1,2}{j,1},'tpp')
        FseofAll{1,2}(j,:)=[];
        j=j-1;
    end
    j=j+1;
end

[ Amp, KO, Amp_HO, KO_HO, AmpKO_2, AmpKO_3] = deal(struct('Score', [], 'Intervention', [], 'prdtConc', [],...
'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [], 'ScoreProduct', []));
i=1;

Amp = testresults(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},'amp');
KO = testresults(Community,{pairedModelInfo(i,:)},options,FseofAll{i,3},'ko');

Amp_HO = testresultsSimilarTargets(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},{'amp'},3); %multiple amps (A)
KO_HO = testresultsSimilarTargets(Community,{pairedModelInfo(i,:)},options,FseofAll{i,3},{'ko'},3); %multiple KOs (K)

AmpKO_2 = testresultsHigherOrder(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},FseofAll{i,3},2,3); %Amp+KO (AK)
AmpKO_3 = testresultsHigherOrder(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},FseofAll{i,3},3,3); %Amp+Amp+KO, Amp+KO+KO (AAK,AKK)

Amp.invTypeMax(1,1:length(Amp.Score))= {'A'};
KO.invTypeMax(1,1:length(KO.Score))= {'K'};

Targets = horzcat(Amp.Intervention,KO.Intervention,Amp_HO.Intervention,KO_HO.Intervention,AmpKO_2.Intervention,AmpKO_3.Intervention);

structs = [Amp KO Amp_HO KO_HO AmpKO_2 AmpKO_3];
fields = fieldnames(Amp);

for i = 1:numel(fields)
    TargetsAll.(fields{i}) = horzcat(structs.(fields{i}));
end

if ~isempty(TargetsAll.Score)
    for i= 1:length(TargetsAll.Intervention)
        %to get length of the Intervetion
        if iscell(TargetsAll.Intervention{i})
            TargetLength(i) = numel(TargetsAll.Intervention{i});     % Length of nested cell array
        else
            TargetLength(i) = 1;                  % Treat non-cell content as single unit
        end
        
        if TargetLength(i) == 1
            TargetsScoreTable{i,1} = TargetsAll.Intervention{1,i};
            TargetsScoreTable{i,2} = {};
            TargetsScoreTable{i,3} = {};
            
        elseif TargetLength(i) == 2
            TargetsScoreTable{i,1} = TargetsAll.Intervention{1,i}{1,1};
            TargetsScoreTable{i,2} = TargetsAll.Intervention{1,i}{1,2};
            TargetsScoreTable{i,3} = {};
            
        elseif TargetLength(i) == 3
            TargetsScoreTable{i,1} = TargetsAll.Intervention{1,i}{1,1};
            TargetsScoreTable{i,2} = TargetsAll.Intervention{1,i}{1,2};
            TargetsScoreTable{i,3} = TargetsAll.Intervention{1,i}{1,3};
        end
        TargetsScoreTable{i,4} = TargetsAll.Score(1,i);
        TargetsScoreTable{i,5} = TargetsAll.ScoreProduct(1,i);
        TargetsScoreTable{i,6} = TargetsAll.invTypeMax(1,i);
        TargetsScoreTable{i,7} = TargetsAll.mutantBiomass1{1,i};
        TargetsScoreTable{i,8} = TargetsAll.mutantBiomass2{1,i};
        TargetsScoreTable{i,9} = TargetsAll.mutantAbundance{1,i};
        TargetsScoreTable{i,10} = TargetsAll.prdtConc{1,i};
    end
    TargetsScoreTable = sortrows(TargetsScoreTable,4,'descend');
    header = {'Intervention1','Intervention2','Intervention3','Score','ScoreProduct','Type of intervention','Mutant biomass flux Org1','Mutant biomass flux Org2','Mutant Abundance','Product Concentration'};
    TargetsScoreTable = [header;TargetsScoreTable];
else
    disp ('No targets found');
    TargetsScoreTable = [];
end
toc
