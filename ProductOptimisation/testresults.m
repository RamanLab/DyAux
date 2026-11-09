function [Result]=testresults(model,pairedModelInfo,options,testRxns,var)

%wild type FBA solutions
CommunityResultRaw = dFBAComFSEOF(model,options);
CommunityResultWT = compileComFSEOF(model,pairedModelInfo,options,CommunityResultRaw);
[prdtMax, prdtMin, prdtMaxDiff, prdtMinDiff] = deal(zeros(length(testRxns),1));
[biomass1,biomass2,abundance] = deal(cell(1,length(testRxns)));
MinMetWT = CommunityResultWT.totPrdtMinConc;
MaxMetWT = CommunityResultWT.totPrdtMaxConc;

for i=1:length(testRxns)
    
    %for amplification targets set value to maxflux
    if strcmp(var,'amp')
        index = find(strcmp(model.rxns,testRxns(i)));
        %to change objective of model
        modelMut=model;
        modelMut.c(:) = 0;
        modelMut.infoCom.spBm(:)=[];
        modelMut.indCom.spBm(:)=[];
        modelMut.infoCom.spBm(find(contains(modelMut.infoCom.spAbbr,extractBefore(testRxns(i),'_'))==1))=testRxns(i);
        modelMut.indCom.spBm(find(contains(modelMut.infoCom.spAbbr,extractBefore(testRxns(i),'_'))==1))=find(strcmp(modelMut.rxns,testRxns(i)));
        modelMut.c(modelMut.indCom.spBm(modelMut.indCom.spBm(:,:)~=0))=1;
        MutResult = optimizeCbModel(modelMut,'max','one');
        if MutResult.stat == 1 && MutResult.x(index)>1e-3
            modeltest = changeRxnBounds(model,testRxns(i),0.5*MutResult.x(index)-1e-3,'l');
            modeltest = changeRxnBounds(modeltest,testRxns(i),0.5*MutResult.x(index)+1e-3,'u');
            
            TestResultRaw = dFBAComFSEOF(modeltest,options);
            TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
        else
            TestResultRaw = []; TestResult = [];
        end
        if ~isempty(TestResult)
            prdtMax(i) = TestResult.totPrdtMaxConc;
            prdtMin(i) = TestResult.totPrdtMinConc;
            prdtMaxDiff(i) = TestResult.totPrdtMaxConc - CommunityResultWT.totPrdtMaxConc;
            prdtMinDiff(i) = TestResult.totPrdtMinConc - CommunityResultWT.totPrdtMinConc;
            biomass1{i} = TestResult.biomass(1);
            biomass2{i} = TestResult.biomass(2);
            abundance{i} = TestResult.abundance;
        end
        %for knockout targets set value to zero
    elseif strcmp(var,'ko')
        modeltest = changeRxnBounds(model,testRxns(i),0,'b');
        TestResultRaw = dFBAComFSEOF(modeltest,options);
        TestResult = compileComFSEOF(modeltest,pairedModelInfo,options,TestResultRaw);
        if ~isempty(TestResult)
            prdtMax(i) = TestResult.totPrdtMaxConc;
            prdtMin(i) = TestResult.totPrdtMinConc;
            prdtMaxDiff(i) = TestResult.totPrdtMaxConc - CommunityResultWT.totPrdtMaxConc;
            prdtMinDiff(i) = TestResult.totPrdtMinConc - CommunityResultWT.totPrdtMinConc;
            biomass1{i} = TestResult.biomass(1);
            biomass2{i} = TestResult.biomass(2);
            abundance{i} = TestResult.abundance;
        end
    end
end

%choose targets with non-zero flux improvement
temp=0;
for i=1:length(testRxns)
    condMax = (prdtMax(i) > 0.01) && (prdtMax(i) > MaxMetWT + 0.05*abs(MaxMetWT));
    condMin = (prdtMin(i) > 0.01) && (prdtMin(i) > MinMetWT + 0.05*abs(MinMetWT));
    condAbd = all(abundance{i} > (options.abdCutoff/100));
    
    if (condMax || condMin) && condAbd
        temp=temp+1;
        Result.Intervention{temp} = testRxns(i);
        
        Result.prdtMax(temp) = prdtMax(i);
        Result.prdtMin(temp) = prdtMin(i);
        
        Result.prdtMaxDiff(temp) = prdtMaxDiff(i);
        Result.prdtMinDiff(temp) = prdtMinDiff(i);
        
        Result.mutantBiomass1(temp) = biomass1{i};
        Result.mutantBiomass2(temp) = biomass2{i};
        Result.mutantAbundance{temp} = abundance{i};
        
        Result.RawAbd(temp) = 1 - abs(diff(abundance{i}));
        Result.RawPrdt(temp) = (options.ConcRangeWeights(1)*prdtMaxDiff(i)) + (options.ConcRangeWeights(2)*prdtMinDiff(i));
        Result.RawBm(temp) = (biomass1{i} + biomass2{i}) / sum(CommunityResultWT.biomass);
    end
    
    if exist('Result','var')
    else
        Result = struct('Intervention', [], 'prdtMax', [], 'prdtMaxDiff', [],'prdtMin', [],'prdtMinDiff', [],'mutantBiomass1', [], 'mutantBiomass2', [], 'mutantAbundance', [],'invTypeMax', [],'RawAbd', [],'RawPrdt', [],'RawBm', []);
    end
    
end


if strcmp(var,'amp') && exist('Result','var')
    Result.invTypeMax = repmat({'A'}, 1, length(Result.Intervention));
elseif exist('Result','var')
    Result.invTypeMax = repmat({'K'}, 1, length(Result.Intervention));
end

end
