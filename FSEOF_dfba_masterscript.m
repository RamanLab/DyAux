tic
Ecoli = readCbModel('iML1515.mat'); Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');% Bsubtilis = readCbModel('iYO844_subtilis.mat');

% define medium
options.mediumMets = {'glc__D_e[u]';'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'ser__L_e[u]';'so3_e[u]'};

% %define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;


speciesList = {Ecoli,Ecoli};
biomassList = {'BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_Ec_iML1515_core_75p37M'};
spNameList = {'EcoliA','EcoliB'};


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
options.ConcRangeWeights = [0.5, 0.5];
options.ScoreWeights = [1/3, 1/3, 1/3];

options.abdCutoff = 10;
options.Products={'EX_spmd_e'};
Community = Communities{1,1};

dels = {'EcoliA_ARGSL','EcoliB_CHORM'};

% %to create auxotrophic pair
[Community] = changeRxnBounds(Community,dels,0,'b');

%to set medium
Community = setMediumCom(Community, options.mediumMets,options.initMedium);


%Dynamic FBA for all communities

allResults = {};
for i= 1 : size(Communities,2)
    disp (i);
    allResults{i} = dFBAComFSEOF(Community,options);
    CommunityResultsWT{1,i} = compileComFSEOF(Community,pairedModelInfo,options,allResults{i});
end

% perform FSEOF
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

[ Amp, KO, Amp_HO, KO_HO, AmpKO_2, AmpKO_3] = deal(struct('Intervention', [], 'prdtMax', [], 'prdtMaxDiff', [],'prdtMin', [],'prdtMinDiff', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'RawAbd', [],'RawPrdt', [],'RawBm', []));
i=1;

if ~isempty(FseofAll{i,2})
    Amp = testresults(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},'amp');
    Amp_HO = testresultsSimilarTargets(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},{'amp'},3); %multiple amps (A)
elseif ~isempty(FseofAll{i,3})
    KO = testresults(Community,{pairedModelInfo(i,:)},options,FseofAll{i,3},'ko');
    KO_HO = testresultsSimilarTargets(Community,{pairedModelInfo(i,:)},options,FseofAll{i,3},{'ko'},3); %multiple KOs (K)
end

if ~isempty(FseofAll{i,2}) && ~isempty(FseofAll{i,3})
    AmpKO_2 = testresultsHigherOrder(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},FseofAll{i,3},2,3); %Amp+KO (AK)
    AmpKO_3 = testresultsHigherOrder(Community,{pairedModelInfo(i,:)},options,FseofAll{i,2},FseofAll{i,3},3,3); %Amp+Amp+KO, Amp+KO+KO (AAK,AKK)
end

allSets = {Amp.Intervention, KO.Intervention, Amp_HO.Intervention, ...
    KO_HO.Intervention, AmpKO_2.Intervention,  AmpKO_3.Intervention};

Targets = horzcat(allSets{~cellfun(@isempty, allSets)});

structs = [Amp KO Amp_HO KO_HO AmpKO_2 AmpKO_3 ];
fields = fieldnames(Amp);

for i = 1:numel(fields)
    TargetsAll.(fields{i}) = horzcat(structs.(fields{i}));
end

if ~isempty(TargetsAll.Intervention)
    for i= 1:length(TargetsAll.Intervention)
        %         to get length of the Intervetion
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
        
        %         Score the interventions
        TargetsAll.Score(1,i) = options.ScoreWeights(1)*TargetsAll.RawAbd(i) + options.ScoreWeights(2)*TargetsAll.RawPrdt(i) + options.ScoreWeights(3)*TargetsAll.RawBm(i);
        TargetsScoreTable{i,4} = TargetsAll.Score(1,i);
        TargetsScoreTable{i,5} = TargetsAll.invTypeMax(1,i);
        TargetsScoreTable{i,6} = TargetsAll.mutantBiomass1(1,i);
        TargetsScoreTable{i,7} = TargetsAll.mutantBiomass2(1,i);
        TargetsScoreTable{i,8} = TargetsAll.mutantAbundance{1,i};
        TargetsScoreTable{i,9} = TargetsAll.prdtMin(1,i);
        TargetsScoreTable{i,10} = TargetsAll.prdtMax(1,i);
    end
    % sort by score, best first
    TargetsScoreTable = sortrows(TargetsScoreTable, 4, 'descend');
    
    % rank with tolerance: scores within rankTol of a group's top score share a rank
    rankTol     = 1e-3;
    sortedScore = cell2mat(TargetsScoreTable(:,4));    % already descending
    rankVec     = zeros(length(TargetsAll.Intervention), 1);
    groupTop    = sortedScore(1);                      % top score of the current rank group
    currRank    = 1;
    for i = 1:length(TargetsAll.Intervention)
        if groupTop - sortedScore(i) > rankTol         % outside the tolerance band
            currRank = currRank + 1;                   % start a new rank group
            groupTop = sortedScore(i);
        end
        rankVec(i) = currRank;
    end
    TargetsScoreTable(:,11) = num2cell(rankVec);       % rank as last column
    
    header = {'Intervention1','Intervention2','Intervention3','Score','Type of intervention', ...
        'Mutant biomass flux Org1','Mutant biomass flux Org2','Mutant Abundance', ...
        'Minimum Product Concentration','Maximum Product Concentration','Rank'};
    TargetsScoreTable = [header; TargetsScoreTable];
else
    disp('No targets found');
    TargetsScoreTable = [];
end
