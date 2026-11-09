%Uses fastFVA to test the effect of mixture of amplification/deletion of any target on
%the flux of specific products and biomass

function [Result]=testresultsHigherOrder(model,pairedModelInfo,options,AmpTargetsAll,KoTargetsAll,intrvSize,SampleSize)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%ProductRxnsList: list of products for which the analysis is done
%AmpTargets: the potential amp target to be tested
%KOTargets: the potential KO target to be tested
%intrvSize: size of intervention
%Result: array of intervention, min and max mutant flux of products,mutant
%biomass flux, wild-type product fluxes

[Intervention,invType]= deal({});
[prdtMin,prdtMax,prdtMaxDiff,prdtMinDiff,biomass1,biomass2,abundance] = deal([]);
MaxSize = 100; %Maximum target combinations to be assessed

%wild type FBA solutions
CommunityResultRaw = dFBAComFSEOF(model,options);
CommunityResultWT = compileComFSEOF(model,pairedModelInfo,options,CommunityResultRaw);
MinMetWT = CommunityResultWT.totPrdtMinConc;
MaxMetWT = CommunityResultWT.totPrdtMaxConc;

%% Amp + KO
if intrvSize == 2
    inv=0; temp=0; stopLoop = false;
    %choose random targets
    randomIndices = randperm(size(AmpTargetsAll, 1));
    AmpTargets = AmpTargetsAll(randomIndices, :);
    
    randomIndices = randperm(size(KoTargetsAll, 1));
    KoTargets = KoTargetsAll(randomIndices, :);
    
    for i=1:length(AmpTargets)
        modeltest = model;
        modelMut = model; count=1;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.stat == 1 && MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
            modelAmp = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modelAmp = changeRxnBounds(modelAmp, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            count = 0;
            continue; %because one of the rxns can't be amplified
        end
        for j= 1:length(KoTargets)
            
            if count>0
                modeltest = modelAmp;
                modeltest = changeRxnBounds(modeltest, KoTargets(j), 0, 'b');
                
                TestResultRaw = dFBAComFSEOF(modeltest,options);
                TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
            else
                TestResultRaw = []; TestResult = [];
            end
            
            if ~isempty(TestResult) && count > 0
                inv=inv+1;
                biomass1(inv) = TestResult.biomass(1);
                biomass2(inv) = TestResult.biomass(2);
                abundance{inv} = TestResult.abundance;
                prdtMin(inv) = TestResult.totPrdtMinConc;
                prdtMax(inv) = TestResult.totPrdtMaxConc;
                prdtMaxDiff(inv) = TestResult.totPrdtMaxConc - MaxMetWT;
                prdtMinDiff(inv) = TestResult.totPrdtMinConc - MinMetWT;
                Intervention{inv} = [AmpTargets(i),KoTargets(j)];
                invType(inv) = {'AK'};
                
                
                condMax = (prdtMax(inv) > 0.01) && (prdtMax(inv) > MaxMetWT + 0.05*abs(MaxMetWT));
                condMin = (prdtMin(inv) > 0.01) && (prdtMin(inv) > MinMetWT + 0.05*abs(MinMetWT));
                condAbd = all(abundance{inv} > (options.abdCutoff/100));
                
                if (condMax || condMin) && condAbd
                    temp=temp+1;
                    Result.Intervention(temp) = Intervention(inv);
                    
                    Result.prdtMin(temp)=prdtMin(inv);
                    Result.prdtMax(temp)=prdtMax(inv);
                    
                    Result.prdtMaxDiff(temp) = prdtMaxDiff(inv);
                    Result.prdtMinDiff(temp) = prdtMinDiff(inv);
                    
                    Result.mutantBiomass1(temp) = biomass1(inv);
                    Result.mutantBiomass2(temp) = biomass2(inv);
                    
                    Result.mutantAbundance(temp) = abundance(inv);
                    Result.invTypeMax(temp) = invType(inv);
                    
                    Result.RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                    Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(inv)) + (options.ConcRangeWeights(2)*prdtMinDiff(inv));
                    Result.RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
                    clearAllMemoizedCaches
                end
            end
            if temp >= SampleSize
                stopLoop = true;
                break;
            elseif inv >= MaxSize
                stopLoop = true;
                break;
            end
        end
        if stopLoop
            break;
        end
    end
    
    if exist('Result','var')
    else
        Result = struct('Intervention', [], 'prdtMax', [], 'prdtMaxDiff', [],'prdtMin', [],'prdtMinDiff', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'RawAbd', [],'RawPrdt', [],'RawBm', []);
    end
    
elseif intrvSize == 3
    %% AMP + KO + KO
    inv=0;stopLoop = false; temp=0;
    
    randomIndices = randperm(size(AmpTargetsAll, 1));
    AmpTargets = AmpTargetsAll(randomIndices, :);
    
    randomIndices = randperm(size(KoTargetsAll, 1));
    KoTargets = KoTargetsAll(randomIndices, :);
    
    for i=1:length(AmpTargets)
        modeltest=model;
        modelMut = model; count=1;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore(AmpTargets(i), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.stat == 1 && MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
            modelAmp = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modelAmp = changeRxnBounds(modelAmp, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            count =0;
            continue; %because one of the rxns can't be amplified
        end
        
        for j= 1:length(KoTargets)-1
            
            if count>0
                modeltest = modelAmp;
                modeltest = changeRxnBounds(modeltest,KoTargets(j:j+1),0,'b');
                
                TestResultRaw = dFBAComFSEOF(modeltest,options);
                TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
            else
                TestResultRaw = []; TestResult = [];
            end
            if ~isempty(TestResult)
                inv=inv+1;
                biomass1(inv) = TestResult.biomass(1);
                biomass2(inv) = TestResult.biomass(2);
                abundance{inv} = TestResult.abundance;
                prdtMin(inv) = TestResult.totPrdtMinConc;
                prdtMax(inv) = TestResult.totPrdtMaxConc;
                prdtMaxDiff(inv) = TestResult.totPrdtMaxConc - MaxMetWT;
                prdtMinDiff(inv) = TestResult.totPrdtMinConc - MinMetWT;
                Intervention{inv} = [AmpTargets(i),KoTargets(j),KoTargets(j+1)];
                invType(inv) = {'AKK'};
                
                condMax = (prdtMax(inv) > 0.01) && (prdtMax(inv) > MaxMetWT + 0.05*abs(MaxMetWT));
                condMin = (prdtMin(inv) > 0.01) && (prdtMin(inv) > MinMetWT + 0.05*abs(MinMetWT));
                condAbd = all(abundance{inv} > (options.abdCutoff/100));
                
                if (condMax || condMin) && condAbd
                    temp=temp+1;
                    Result.Intervention(temp) = Intervention(inv);
                    Result.prdtMin(temp)=prdtMin(inv);
                    Result.prdtMax(temp)=prdtMax(inv);
                    Result.prdtMaxDiff(temp) = prdtMaxDiff(inv);
                    Result.prdtMinDiff(temp) = prdtMinDiff(inv);
                    Result.mutantBiomass1(temp) = biomass1(inv);
                    Result.mutantBiomass2(temp) = biomass2(inv);
                    Result.mutantAbundance(temp) = abundance(inv);
                    Result.invTypeMax(temp) = invType(inv);
                    Result.RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                    Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(inv)) + (options.ConcRangeWeights(2)*prdtMinDiff(inv));
                    Result.RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
                end
                clearAllMemoizedCaches
            end
            if temp >= SampleSize
                stopLoop = true;
                break;
            elseif inv >= MaxSize
                stopLoop = true;
                break;
            end
        end
        if stopLoop
            break;
        end
    end
    
    %% AMP + AMP + KO
    stopLoop = false;
    for i=1:length(AmpTargets)-1
        modeltest=model; count =1;
        
        %First amplification
        modelMut = model;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.stat == 1 && MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3
            modelAmp = changeRxnBounds(modeltest, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modelAmp = changeRxnBounds(modelAmp, AmpTargets(i), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            count=0;
            continue; %because one of the rxns can't be amplified
        end
        %second amplification
        modelMut = model;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        spIdx = find(contains(modelMut.infoCom.spAbbr, extractBefore( AmpTargets(i+1), '_')) == 1);
        modelMut.infoCom.spBm(spIdx) = AmpTargets(i+1);
        modelMut.indCom.spBm(spIdx) = find(strcmp(modelMut.rxns, AmpTargets(i+1)));
        modelMut.c(find(strcmp(modelMut.rxns, AmpTargets(i+1)))) = 1;
        
        MutResult = optimizeCbModel(modelMut, 'max', 'one');
        if MutResult.stat == 1 && MutResult.x(modelMut.indCom.spBm(spIdx))>1e-3 && count>0
            modelAmp = changeRxnBounds(modelAmp, AmpTargets(i+1), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) - 1e-3, 'l');
            modelAmp = changeRxnBounds(modelAmp, AmpTargets(i+1), 0.5 * MutResult.x(modelMut.indCom.spBm(spIdx)) + 1e-3, 'u');
        else
            count = 0;
            continue; %because one of the rxns can't be amplified
        end
        
        for j= 1:length(KoTargets)
            
            if count>0
                modeltest = modelAmp;
                modeltest = changeRxnBounds(modeltest,KoTargets(j),0,'b');
                
                TestResultRaw = dFBAComFSEOF(modeltest,options);
                TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
            else
                TestResult = [];
            end
            
            if ~isempty(TestResult) && count > 0
                inv=inv+1;
                biomass1(inv) = TestResult.biomass(1);
                biomass2(inv) = TestResult.biomass(2);
                abundance{inv} = TestResult.abundance;
                prdtMin(inv) = TestResult.totPrdtMinConc;
                prdtMax(inv) = TestResult.totPrdtMaxConc;
                prdtMaxDiff(inv) = TestResult.totPrdtMaxConc - MaxMetWT;
                prdtMinDiff(inv) = TestResult.totPrdtMinConc - MinMetWT;
                
                Intervention{inv} = [AmpTargets(i),AmpTargets(i+1),KoTargets(j)];
                invType(inv) = {'AAK'};
                
                condMax = (prdtMax(inv) > 0.01) && (prdtMax(inv) > MaxMetWT + 0.05*abs(MaxMetWT));
                condMin = (prdtMin(inv) > 0.01) && (prdtMin(inv) > MinMetWT + 0.05*abs(MinMetWT));
                condAbd = all(abundance{inv} > (options.abdCutoff/100));
                
                if (condMax || condMin) && condAbd
                    temp=temp+1;
                    Result.Intervention(temp) = Intervention(inv);
                    Result.prdtMin(temp)=prdtMin(inv);
                    Result.prdtMax(temp)=prdtMax(inv);
                    Result.prdtMaxDiff(temp) = prdtMaxDiff(inv);
                    Result.prdtMinDiff(temp) = prdtMinDiff(inv);
                    Result.mutantBiomass1(temp) = biomass1(inv);
                    Result.mutantBiomass2(temp) = biomass2(inv);
                    Result.mutantAbundance(temp) = abundance(inv);
                    Result.invTypeMax(temp) = invType(inv);
                    Result.RawAbd(temp) = 1 - abs(diff(abundance{inv}));
                    Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(inv)) + (options.ConcRangeWeights(2)*prdtMinDiff(inv));
                    Result.RawBm(temp) = (biomass1(inv) + biomass2(inv)) / sum(CommunityResultWT.biomass);
                end
                clearAllMemoizedCaches
            end
            if temp >= (SampleSize)
                stopLoop = true;
                break;
            elseif inv >= MaxSize
                stopLoop = true;
                break;
            end
        end
        if stopLoop
            break;
        end
    end
    
    if exist('Result','var')
    else
        Result = struct('Intervention', [], 'prdtMax', [], 'prdtMaxDiff', [],'prdtMin', [],'prdtMinDiff', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'RawAbd', [],'RawPrdt', [],'RawBm', []);
    end
    
end

end
