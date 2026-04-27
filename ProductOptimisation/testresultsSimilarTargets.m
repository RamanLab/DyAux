function [sortedResults] = testresultsSimilarTargets(model, pairedModelInfo, options, TargetRxns, var,SampleSize)

%%%% Input and Output Parameters
% model: the GSMM with appropriate medium bounds applied
% ProductRxnsList: list of products for which the analysis is done
% TargetRxns: the potential amp/ko target to be tested
% var: nature of intervention - amp/ko
% Result: array of intervention, min and max mutant flux of products, mutant
% biomass flux, wild-type product fluxes

[Intervention, invType, invTypeMax, abundance] = deal({});
[MinMetMut, MaxMetMut, biomass1, biomass2] = deal([]);
% prdtID = find(contains(model.infoCom.EXcom, options.Products) == 1);

%Prdt conc for secreted metabolites in wild type
CommunityResultRaw = dFBAComFSEOF(model,options);
CommunityResultWT = compileComFSEOF(model,pairedModelInfo,options,CommunityResultRaw);
MinMetWT = CommunityResultWT.totPrdtFVAMinConc;
MaxMetWT = CommunityResultWT.totPrdtFVAMaxConc;
MetWT = CommunityResultWT.totPrdtMaxConc;

% Store environment for later restoration
MaxSize=100; temp=0;


%% Amplification of target reactions
if strcmp(var, 'amp')
    for intrvSize = 2:3
        count=0;
        AmpTargetsAll = nchoosek(TargetRxns, intrvSize);
        
        % randomize the target order
        randOrder = randperm(size(AmpTargetsAll,1));
        randOrder = randOrder(1:min(size(AmpTargetsAll,1),MaxSize));
        numTargets = length(randOrder);

        
        % Preallocate arrays
        biomass1 = zeros(1, numTargets);
        biomass2 = zeros(1, numTargets);
        MinMetMut = zeros(1, numTargets);
        MaxMetMut = zeros(1, numTargets);
        abundance = cell(1, numTargets);
        Intervention = cell(1, numTargets);
        invType = cell(1, numTargets);
        
        
        % Loop for amplifications
        for i = 1:double(numTargets)
            % Restore environment and set solver
            solverOK = changeCobraSolver(options.solver, 'all');
            modeltest = model; idx = randOrder(i);
            
            % Create mutant model for amplification
            for j=1: size(AmpTargetsAll,2)
                modelMut = model;
                modelMut.c(:) = 0;
                modelMut.infoCom.spBm(:)=[];
                modelMut.indCom.spBm(:)=[];
                spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargetsAll(idx,j), '_')) == 1);
                modelMut.infoCom.spBm(spIdx) = AmpTargetsAll(idx,j);
                modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargetsAll(idx,j)));
                modelMut.c(find(strcmp(modelMut.rxns, AmpTargetsAll(idx,j)))) = 1;
                
                MutResult = optimizeCbModel(modelMut, 'max', 'one');
                if MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
                    modeltest = changeRxnBounds(modeltest, AmpTargetsAll(idx,j), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
                    modeltest = changeRxnBounds(modeltest, AmpTargetsAll(idx,j), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
                else
                    modeltest = model;
                    continue; %because one of the rxns can't be amplified
                end
            end
            TestResultRaw = dFBAComFSEOF(modeltest, options);
            TestResult = compileComFSEOF(modeltest, pairedModelInfo, options, TestResultRaw);
            
            if ~isempty(TestResult)
                biomass1(i) = TestResult.biomass(1);
                biomass2(i) = TestResult.biomass(2);
                MinMetMut(i) = TestResult.totPrdtMinConc;
                MaxMetMut(i) = TestResult.totPrdtMaxConc;
                abundance{i} = TestResult.abundance;
                Intervention{i} = AmpTargetsAll(idx,:);
                invType{i} = 'A';
            end
            
            if MaxMetMut(i) > 0.01 && MaxMetMut(i) > 1.05 * MaxMetWT(i) && all(abundance{i}) > (options.abdCutoff / 100)
                temp = temp + 1; count= count+1;
                Result.InterventionNew{temp} = Intervention{1,i};
                Result.MinMetConc(temp) = num2cell(MinMetMut(i));
                Result.MaxMetConc(temp) = num2cell(MaxMetMut(i));
                Result.mutantBiomass1(temp) = num2cell(biomass1(i));
                Result.mutantBiomass2(temp) = num2cell(biomass2(i));
                Result.abundance(temp) = abundance(i);
                Result.invTypeMax(temp) = invType(i);
                RawAbd(temp) = 1 - abs(diff(abundance{i}));
                RawPrdt(temp) = MaxMetMut(i) / MetWT;
                RawBm(temp) = (biomass1(i) + biomass2(i)) / sum(CommunityResultWT.biomass);
            end
            if count >= SampleSize
                 break;
            end
        end
    end
    
%% Deletion of target reactions
elseif strcmp(var, 'ko')
    for intrvSize = 2:3
        count = 0;
        KoTargetsAll = nchoosek(TargetRxns, intrvSize);
        % randomize the target order
        randOrder = randperm(size(KoTargetsAll,1));
        randOrder = randOrder(1:min(size(KoTargetsAll,1),MaxSize));
        numTargets = length(randOrder);

        % Preallocate arrays
        biomass1 = zeros(1, numTargets);
        biomass2 = zeros(1, numTargets);
        MinMetMut = zeros(1, numTargets);
        MaxMetMut = zeros(1, numTargets);
        abundance = cell(1, numTargets);
        Intervention = cell(1, numTargets);
        invType = cell(1, numTargets);
        
        
        for i = 1:numTargets
            idx = randOrder(i);
                        
            % Create mutant model for knockout
            modeltest = changeRxnBounds(model, KoTargetsAll(idx,:), 0, 'b');
            TestResultRaw = dFBAComFSEOF(modeltest, options);
            TestResult = compileComFSEOF(modeltest, pairedModelInfo, options, TestResultRaw);
            
            if ~isempty(TestResult)
                biomass1(i) = TestResult.biomass(1);
                biomass2(i) = TestResult.biomass(2);
                MinMetMut(i) = TestResult.totPrdtMinConc;
                MaxMetMut(i) = TestResult.totPrdtMaxConc;
                abundance{i} = TestResult.abundance;
                Intervention{i} = KoTargetsAll(idx,:);
                invType{i} = 'K';
            end
            if MaxMetMut(i) > 0.01 && MaxMetMut(i) > 1.05 * MaxMetWT(i) && all(abundance{i}) > (options.abdCutoff / 100)
                temp = temp + 1; count= count+1;
                Result.InterventionNew{temp} = Intervention{1,i};
                Result.MinMetConc(temp) = num2cell(MinMetMut(i));
                Result.MaxMetConc(temp) = num2cell(MaxMetMut(i));
                Result.mutantBiomass1(temp) = num2cell(biomass1(i));
                Result.mutantBiomass2(temp) = num2cell(biomass2(i));
                Result.abundance(temp) = abundance(i);
                Result.invTypeMax(temp) = invType(i);
                RawAbd(temp) = 1 - abs(diff(abundance{i}));
                RawPrdt(temp) = MaxMetMut(i) / MetWT;
                RawBm(temp) = (biomass1(i) + biomass2(i)) / sum(CommunityResultWT.biomass);
            end
            if count >= SampleSize
                 break;
            end
        end
        
    end
end

%% compile results
if exist('Result','var')
    % Normalize each component to [0, 1]
    normAbd = (RawAbd - min(RawAbd)) / (max(RawAbd) - min(RawAbd));
    normPrdt = (RawPrdt - min(RawPrdt)) / (max(RawPrdt) - min(RawPrdt));
    normBm = (RawBm - min(RawBm)) / (max(RawBm) - min(RawBm));
    
    % Recompute final scores using normalized components
    for k = 1:length(normAbd)
        Result.Score(k) = options.ScoreWeights(1)*normAbd(k) + options.ScoreWeights(2)*normPrdt(k) + options.ScoreWeights(3)*normBm(k);
        Result.ScoreProduct(k) = normAbd(k) * normPrdt(k) * normBm(k);
    end
    for i = 1:numel(Result)
        [~, sortIdx] = sort(Result(i).Score, 'descend');
        
        % Now apply same sorting to every field
        sortedResults(i).Score = Result(i).Score(sortIdx);
        sortedResults(i).Intervention     = Result(i).InterventionNew(sortIdx);
        sortedResults(i).prdtConc     = Result(i).MaxMetConc(sortIdx);
        sortedResults(i).mutantBiomass1     = Result(i).mutantBiomass1(sortIdx);
        sortedResults(i).mutantBiomass2    = Result(i).mutantBiomass2(sortIdx);
        sortedResults(i).mutantAbundance     = Result(i).abundance(sortIdx);
        sortedResults(i).invTypeMax     = Result(i).invTypeMax(sortIdx);
        sortedResults(i).ScoreProduct = Result(i).ScoreProduct(sortIdx);
    end
else
    sortedResults = struct('Score', [], 'Intervention', [], 'prdtConc', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'ScoreProduct', []);
end

clearAllMemoizedCaches;
end
