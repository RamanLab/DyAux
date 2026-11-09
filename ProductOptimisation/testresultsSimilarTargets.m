function [Result] = testresultsSimilarTargets(model, pairedModelInfo, options, TargetRxns, var,SampleSize)

%%%% Input and Output Parameters
% model: the GSMM with appropriate medium bounds applied
% ProductRxnsList: list of products for which the analysis is done
% TargetRxns: the potential amp/ko target to be tested
% var: nature of intervention - amp/ko
% Result: array of intervention, min and max mutant flux of products, mutant
% biomass flux, wild-type product fluxes


%Prdt conc for secreted metabolites in wild type
CommunityResultRaw = dFBAComFSEOF(model,options);
CommunityResultWT = compileComFSEOF(model,pairedModelInfo,options,CommunityResultRaw);
MinMetWT = CommunityResultWT.totPrdtMinConc;
MaxMetWT = CommunityResultWT.totPrdtMaxConc;

MaxSize=100; %Maximum target combinations to be assessed
temp=0;


%% Amplification of target reactions
if strcmp(var, 'amp')
    for intrvSize = 2:3
        count=1;
        AmpTargetsAll = nchoosek(TargetRxns, intrvSize);
        
        % randomize the target order
        randOrder = randperm(size(AmpTargetsAll,1));
        randOrder = randOrder(1:min(size(AmpTargetsAll,1),MaxSize));
        numTargets = length(randOrder);
        
        
        % Preallocate arrays
        [biomass1, biomass2, prdtMin, prdtMax, prdtMinDiff, prdtMaxDiff] = deal(zeros(1, numTargets));
        [abundance, Intervention, invType] = deal(cell(1, numTargets));
        
        % Loop for amplifications
        for i = 1:double(numTargets)
            modeltest = model; idx = randOrder(i);
            count=1; TestResultRaw = []; TestResult = [];
            skipOuter = false; %temp variable
            % Create mutant model for amplification
            for j=1: size(AmpTargetsAll,2)
                modelMut = model;
                modelMut.c(:) = 0;
                modelMut.infoCom.spBm(:)=[];
                modelMut.indCom.spBm(:)=[];
                spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore(AmpTargetsAll(idx,j), '_')) == 1);
                modelMut.infoCom.spBm(spIdx) = AmpTargetsAll(idx,j);
                modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargetsAll(idx,j)));
                modelMut.c(find(strcmp(modelMut.rxns, AmpTargetsAll(idx,j)))) = 1;
                
                MutResult = optimizeCbModel(modelMut, 'max', 'one');
                if MutResult.stat == 1 && MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
                    modeltest = changeRxnBounds(modeltest, AmpTargetsAll(idx,j), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
                    modeltest = changeRxnBounds(modeltest, AmpTargetsAll(idx,j), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
                else
                    count = 0;
                    skipOuter = true; %because one of the rxns can't be amplified
                    break;
                end
            end
            if skipOuter
                continue;
            end
            if count>0
                TestResultRaw = dFBAComFSEOF(modeltest, options);
                TestResult = compileComFSEOF(modeltest, pairedModelInfo, options, TestResultRaw);
            else
                TestResultRaw = []; TestResult = [];
            end
            if ~isempty(TestResult) && count>0
                biomass1(i) = TestResult.biomass(1);
                biomass2(i) = TestResult.biomass(2);
                prdtMin(i) = TestResult.totPrdtMinConc;
                prdtMax(i) = TestResult.totPrdtMaxConc;
                prdtMaxDiff(i) = TestResult.totPrdtMaxConc - MaxMetWT;
                prdtMinDiff(i) = TestResult.totPrdtMinConc - MinMetWT;
                abundance{i} = TestResult.abundance;
                Intervention{i} = AmpTargetsAll(idx,:);
                invType{i} = 'A';
            end
            
            condMax = (prdtMax(i) > 0.01) && (prdtMax(i) > MaxMetWT + 0.05*abs(MaxMetWT));
            condMin = (prdtMin(i) > 0.01) && (prdtMin(i) > MinMetWT + 0.05*abs(MinMetWT));
            condAbd = all(abundance{i} > (options.abdCutoff/100));
            
            if (condMax || condMin) && condAbd
                temp = temp + 1;
                Result.Intervention{temp} = Intervention{1,i};
                
                Result.prdtMax(temp) = prdtMax(i);
                Result.prdtMin(temp) = prdtMin(i);
                
                Result.prdtMaxDiff(temp) = prdtMaxDiff(i);
                Result.prdtMinDiff(temp) = prdtMinDiff(i);
                
                Result.mutantBiomass1(temp) = biomass1(i);
                Result.mutantBiomass2(temp) = biomass2(i);
                
                Result.mutantAbundance(temp) = abundance(i);
                Result.invTypeMax(temp) = invType(i);
                
                Result.RawAbd(temp) = 1 - abs(diff(abundance{i}));
                Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(i)) + (options.ConcRangeWeights(2)*prdtMinDiff(i));
                Result.RawBm(temp) = (biomass1(i) + biomass2(i)) / sum(CommunityResultWT.biomass);
            end
            if temp >= SampleSize
                break;
            end
        end
    end
    
    %% Deletion of target reactions
elseif strcmp(var, 'ko')
    for intrvSize = 2:3
        TestResult=[];
        KoTargetsAll = nchoosek(TargetRxns, intrvSize);
        % randomize the target order
        randOrder = randperm(size(KoTargetsAll,1));
        randOrder = randOrder(1:min(size(KoTargetsAll,1),MaxSize));
        numTargets = length(randOrder);
        
        % Preallocate arrays
        [biomass1, biomass2, prdtMin, prdtMax, prdtMinDiff, prdtMaxDiff] = deal(zeros(1, numTargets));
        [abundance, Intervention, invType] = deal(cell(1, numTargets));
        
        for i = 1:numTargets
            idx = randOrder(i);
            
            % Create mutant model for knockout
            modeltest = changeRxnBounds(model, KoTargetsAll(idx,:), 0, 'b');
            TestResultRaw = dFBAComFSEOF(modeltest, options);
            TestResult = compileComFSEOF(modeltest, pairedModelInfo, options, TestResultRaw);
            
            if ~isempty(TestResult)
                biomass1(i) = TestResult.biomass(1);
                biomass2(i) = TestResult.biomass(2);
                prdtMin(i) = TestResult.totPrdtMinConc;
                prdtMax(i) = TestResult.totPrdtMaxConc;
                prdtMaxDiff(i) = TestResult.totPrdtMaxConc - MaxMetWT;
                prdtMinDiff(i) = TestResult.totPrdtMinConc - MinMetWT;
                abundance{i} = TestResult.abundance;
                Intervention{i} = KoTargetsAll(idx,:);
                invType{i} = 'K';
            end
            condMax = (prdtMax(i) > 0.01) && (prdtMax(i) > MaxMetWT + 0.05*abs(MaxMetWT));
            condMin = (prdtMin(i) > 0.01) && (prdtMin(i) > MinMetWT + 0.05*abs(MinMetWT));
            condAbd = all(abundance{i} > (options.abdCutoff/100));
            
            if (condMax || condMin) && condAbd
                temp = temp + 1;
                Result.Intervention{temp} = Intervention{1,i};
                
                Result.prdtMax(temp) = prdtMax(i);
                Result.prdtMin(temp) = prdtMin(i);
                
                Result.prdtMaxDiff(temp) = prdtMaxDiff(i);
                Result.prdtMinDiff(temp) = prdtMinDiff(i);
                
                Result.mutantBiomass1(temp) = biomass1(i);
                Result.mutantBiomass2(temp) = biomass2(i);
                
                Result.mutantAbundance(temp) = abundance(i);
                Result.invTypeMax(temp) = invType(i);
                
                Result.RawAbd(temp) = 1 - abs(diff(abundance{i}));
                Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(i)) + (options.ConcRangeWeights(2)*prdtMinDiff(i));
                Result.RawBm(temp) = (biomass1(i) + biomass2(i)) / sum(CommunityResultWT.biomass);
            end
            if temp >= SampleSize
                break;
            end
        end
        
    end
end

%% compile results
if exist('Result','var')
else
    Result = struct('Intervention', [], 'prdtMax', [], 'prdtMaxDiff', [],'prdtMin', [],'prdtMinDiff', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'RawAbd', [],'RawPrdt', [],'RawBm', []);
end

clearAllMemoizedCaches;
end
