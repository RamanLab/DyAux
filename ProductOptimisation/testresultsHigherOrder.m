%Uses fastFVA to test the effect of mixture of amplification/deletion of any target on
%the flux of specific products and biomass

function [sortedResults]=testresultsHigherOrder(model,pairedModelInfo,options,AmpTargetsAll,KoTargetsAll,intrvSize,SampleSize)

% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%ProductRxnsList: list of products for which the analysis is done
%AmpTargets: the potential amp target to be tested
%KOTargets: the potential KO target to be tested
%intrvSize: size of intervention
%Result: array of intervention, min and max mutant flux of products,mutant
%biomass flux, wild-type product fluxes

[Intervention,invType]= deal({});
[MinMetMut,MaxMetMut,biomass1,biomass2] = deal([]);
prdtID = find(contains(model.infoCom.EXcom,options.Products)==1);

%wild type FBA solutions
CommunityResultRaw = dFBAComFSEOF(model,options);
CommunityResultWT = compileComFSEOF(model,pairedModelInfo,options,CommunityResultRaw);
MinMetWT = CommunityResultWT.totPrdtFVAMinConc;
MaxMetWT = CommunityResultWT.totPrdtFVAMaxConc;
MetWT = CommunityResultWT.totPrdtMaxConc;

%% Amp + KO
if intrvSize == 2
    inv=0; count=0; temp=0; stopLoop = false;
    %choose random targets
    randomIndices = randperm(size(AmpTargetsAll, 1));
    AmpTargets = AmpTargetsAll(randomIndices, :);
    
    randomIndices = randperm(size(KoTargetsAll, 1));
    KoTargets = KoTargetsAll(randomIndices, :);
    
    for i=1:length(AmpTargets)
        modeltest=model;
        for j= 1:length(KoTargets)
            modelMut = model;
            modelMut.c(:) = 0;
            modelMut.infoCom.spBm(:)=[];
            modelMut.indCom.spBm(:)=[];
            spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i), '_')) == 1);
            modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
            modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
            modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
            
            MutResult = optimizeCbModel(modelMut, 'max', 'one');
            if MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
                modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
                modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
            else
                modeltest = model;
                break; %because one of the rxns can't be amplified
            end
            
            modeltest = changeRxnBounds(model, KoTargets(j), 0, 'b');
            
            TestResultRaw = dFBAComFSEOF(modeltest,options);
            TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
            
            if ~isempty(TestResult)
                inv=inv+1;
                biomass1(inv) = TestResult.biomass(1);
                biomass2(inv) = TestResult.biomass(2);
                abundance{inv} = TestResult.abundance;
                MinMetMut(inv) = TestResult.totPrdtMinConc;
                MaxMetMut(inv) = TestResult.totPrdtMaxConc;
                Intervention{inv} = [AmpTargets(i),KoTargets(j)];
                invType(inv) = {'AK'};
                
                
                if MaxMetMut(inv)>0.01 && MaxMetMut(inv)>1.05*MaxMetWT(inv) && all(abundance{inv})>(options.abdCutoff/100)
                    temp=temp+1; count=count+1;
                    Result.InterventionNew(temp) = Intervention(inv);
                    Result.MinMetConc(temp)=num2cell(MinMetMut(inv));
                    Result.MaxMetConc(temp)=num2cell(MaxMetMut(inv));
                    Result.mutantBiomass1(temp) = num2cell(biomass1(inv));
                    Result.mutantBiomass2(temp) = num2cell(biomass2(inv));
                    Result.abundance(temp) = abundance(inv);
                    Result.invTypeMax(temp) = invType(inv);
                    RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                    RawPrdt(temp) = MaxMetMut(inv) / MetWT;
                    RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
                    clearAllMemoizedCaches
                end
            end
            if temp>= SampleSize
                stopLoop = true;
                break;
            end
        end
        if stopLoop
            break;
        end
    end
    %% choose targets with non-zero flux improvement
    
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
        sortedResults = [];
    end
    
elseif intrvSize == 3
    %% AMP + KO + KO
    inv=0;stopLoop = false; temp=0; count=0;
    
    randomIndices = randperm(size(AmpTargetsAll, 1));
    AmpTargets = AmpTargetsAll(randomIndices, :);
    
    randomIndices = randperm(size(KoTargetsAll, 1));
    KoTargets = KoTargetsAll(randomIndices, :);
    
    for i=1:length(AmpTargets)
        modeltest=model;
        for j= 1:length(KoTargets)-1
            modelMut = model;
            modelMut.c(:) = 0;
            modelMut.infoCom.spBm(:)=[];
            modelMut.indCom.spBm(:)=[];
            spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore(AmpTargets(i), '_')) == 1);
            modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
            modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
            modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
            
            MutResult = optimizeCbModel(modelMut, 'max', 'one');
            if MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
                modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
                modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
            else
                modeltest = model;
                break; %because one of the rxns can't be amplified
            end
            
            modeltest = changeRxnBounds(modeltest,KoTargets(j:j+1),0,'b');
            
            TestResultRaw = dFBAComFSEOF(modeltest,options);
            TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
            
            if ~isempty(TestResult)
                inv=inv+1;
                biomass1(inv) = TestResult.biomass(1);
                biomass2(inv) = TestResult.biomass(2);
                abundance{inv} = TestResult.abundance;
                MinMetMut(inv) = TestResult.totPrdtMinConc;
                MaxMetMut(inv) = TestResult.totPrdtMaxConc;
                Intervention{inv} = [AmpTargets(i),KoTargets(j),KoTargets(j+1)];
                invType(inv) = {'AKK'};
                if MaxMetMut(inv)>0.01 && MaxMetMut(inv)>1.05*MaxMetWT(inv) && all(abundance{inv})>(options.abdCutoff/100)
                    temp=temp+1; count=count+1;
                    Result.InterventionNew(temp) = Intervention(inv);
                    Result.MinMetConc(temp)=num2cell(MinMetMut(inv));
                    Result.MaxMetConc(temp)=num2cell(MaxMetMut(inv));
                    Result.mutantBiomass1(temp) = num2cell(biomass1(inv));
                    Result.mutantBiomass2(temp) = num2cell(biomass2(inv));
                    Result.abundance(temp) = abundance(inv);
                    Result.invTypeMax(temp) = invType(inv);
                    RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                    RawPrdt(temp) = MaxMetMut(inv) / MetWT;
                    RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
                end
                clearAllMemoizedCaches
            end
            if temp>= SampleSize
                stopLoop = true;
                break;
            end
        end
        if stopLoop
            break;
        end
    end

%% AMP + AMP + KO
count=0;stopLoop = false;
for i=1:length(AmpTargets)-1
    modeltest=model;
    for j= 1:length(KoTargets)
        modelMut = model;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
            modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modeltest = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            modeltest = model;
            break; %because one of the rxns can't be amplified
        end
        
        modelMut = model;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i+1), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i+1);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i+1)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i+1)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
            modeltest = changeRxnBounds(modeltest, AmpTargets(i+1), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modeltest = changeRxnBounds(modeltest, AmpTargets(i+1), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            modeltest = model;
            break; %because one of the rxns can't be amplified
        end
        
        modeltest = changeRxnBounds(modeltest,KoTargets(j),0,'b');
        
        TestResultRaw = dFBAComFSEOF(modeltest,options);
        TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
        
        if ~isempty(TestResult)
            inv=inv+1;
            biomass1(inv) = TestResult.biomass(1);
            biomass2(inv) = TestResult.biomass(2);
            abundance{inv} = TestResult.abundance;
            MinMetMut(inv) = TestResult.totPrdtMinConc;
            MaxMetMut(inv) = TestResult.totPrdtMaxConc;
            
            Intervention{inv} = [AmpTargets(i),AmpTargets(i+1),KoTargets(j)];
            invType(inv) = {'AAK'};
            
            if MaxMetMut(inv)>0.01 && MaxMetMut(inv)>1.05*MaxMetWT(inv) && all(abundance{inv})>(options.abdCutoff/100)
                temp=temp+1; count=count+1;
                Result.InterventionNew(temp) = Intervention(inv);
                Result.MinMetConc(temp)=num2cell(MinMetMut(inv));
                Result.MaxMetConc(temp)=num2cell(MaxMetMut(inv));
                Result.mutantBiomass1(temp) = num2cell(biomass1(inv));
                Result.mutantBiomass2(temp) = num2cell(biomass2(inv));
                Result.abundance(temp) = abundance(inv);
                Result.invTypeMax(temp) = invType(inv);
                RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                RawPrdt(temp) = MaxMetMut(inv) / MetWT;
                RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
            end
            clearAllMemoizedCaches
        end
        if temp >= (SampleSize*2)
            stopLoop = true;
            break;
        end
    end
    if stopLoop
        break;
    end
end

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
end
end

