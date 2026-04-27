function [sortedAuxotrophScore,sortedAuxotrophResults] = FindAuxotroph(Community,pairedModelInfo,options,Donor,Size)

options.Products=[]; %to omit running FVA for product
sortedAuxotrophScore =[]; sortedAuxotrophResults=[];

WTResults = dFBAComFSEOF(Community, options);
CommunityWT = compileComFSEOF(Community, pairedModelInfo, options, WTResults);


if Donor>0
    if Donor==1
        singleModelInfo{1,1} = pairedModelInfo{1,4};
        Mono = changeRxnBounds(Community,Community.rxns(strmatch(strcat(pairedModelInfo{1,2},'_'),Community.rxns)), 0,'b');
    elseif Donor==2
        singleModelInfo{1,1} = pairedModelInfo{1,2};
        Mono = changeRxnBounds(Community,Community.rxns(strmatch(strcat(pairedModelInfo{1,4},'_'),Community.rxns)), 0,'b');
    end
    
    
    %to find lethal genes in recepient
    
    grRatio = singleRxnDeletion(Mono);
    singleLethals = Mono.rxns(find(grRatio < 0.01 ));
    
    %remove exchange reactions
    excRxns = Mono.rxns(find(findExcRxns(Mono)==1));
    transRxns = findTransRxns(Mono);
    singleLethals = setdiff(singleLethals,[excRxns;transRxns]);
    
    %to find the genes rescued
    randIndices = randperm(length(singleLethals));
    temp = 0;
    
    for k = 1:length(randIndices)
        i = randIndices(k);
        
        KoRecepient = changeRxnBounds(Community, singleLethals(i),0,'b');
        
        soln= optimizeCbModel(KoRecepient);
        if soln.f>0.001
            allResults = dFBAComFSEOF(KoRecepient, options);
            KoResult = compileComFSEOF(Community, pairedModelInfo, options, allResults);
            
            if min(KoResult.abundance) > 0.105 && min(KoResult.biomass) > 0.105
                temp = temp + 1;
                AuxotrophScore(temp, 1) = singleLethals(i);
                AuxotrophResults(temp) = KoResult;
                AuxotrophScore{temp, 2} = (1 - abs(diff(KoResult.abundance)));
                AuxotrophScore{temp, 3} = AuxotrophResults(temp).biomass;
                AuxotrophScore{temp, 4} = AuxotrophResults(temp).abundance;
            end
        end
        
        if temp >= Size || k> 0.7*length(randIndices)
            break;
        end
    end
    
    %sort the Auxotrophs acc to score
    if exist('AuxotrophScore','var') && size(AuxotrophScore,1)>0
        for i=1:size(AuxotrophScore,1)
            AuxotrophScore{i,2} = AuxotrophScore{i,2}*(sum(AuxotrophResults(i).biomass)/max(cellfun(@sum,{AuxotrophResults(:).biomass})));
        end
        [~, sortIdx] = sort(cell2mat(AuxotrophScore(:,2)), 'descend');
        sortedAuxotrophScore = AuxotrophScore(sortIdx, :);
        AuxotrophResults = rmfield(AuxotrophResults,{'prdtFVAMinConc1','prdtFVAMinConc2','totPrdtFVAMinConc','prdtFVAMaxConc1','prdtFVAMaxConc2','totPrdtFVAMaxConc','prdtID','prdtRxns','prdtMinConc1','prdtMinConc2','totPrdtMinConc','prdtMaxConc1','prdtMaxConc2','totPrdtMaxConc'});
        sortedAuxotrophResults = AuxotrophResults(:,sortIdx);
    end
    
elseif Donor==0
    
    singleModelInfoA{1,1} = pairedModelInfo{1,2};
    MonoA = changeRxnBounds(Community,Community.rxns(strmatch(strcat(pairedModelInfo{1,4},'_'),Community.rxns)), 0,'b');
    
    singleModelInfoB{1,1} = pairedModelInfo{1,4};
    MonoB = changeRxnBounds(Community,Community.rxns(strmatch(strcat(pairedModelInfo{1,2},'_'),Community.rxns)), 0,'b');
    
    %to find lethal genes in recepient
    grRatioA = singleRxnDeletion(MonoA);
    singleLethalsA = MonoA.rxns(find(grRatioA < 0.01 ));
    grRatioB = singleRxnDeletion(MonoB);
    singleLethalsB = MonoB.rxns(find(grRatioB < 0.01 ));
    
    %remove exchange reactions
    excRxnsA = MonoA.rxns(find(findExcRxns(MonoA)==1));
    transRxnsA = findTransRxns(MonoA);
    singleLethalsA = setdiff(singleLethalsA,[excRxnsA;transRxnsA]);
    
    excRxnsB = MonoB.rxns(find(findExcRxns(MonoB)==1));
    transRxnsB = findTransRxns(MonoB);
    singleLethalsB = setdiff(singleLethalsB,[excRxnsB;transRxnsB]);
    
    j=1;
    for i=1:length(singleLethalsA)
        if contains(singleLethalsA{j},'tex')||contains(singleLethalsA{j},'tpp')
            singleLethalsA(j)=[];
            j=j-1;
        end
        j=j+1;
    end

    j=1;
    for i=1:length(singleLethalsB)
        if contains(singleLethalsB{j},'tex') || contains(singleLethalsB{j},'tpp')
            singleLethalsB(j)=[];
            j=j-1;
        end
        j=j+1;
    end
    
    
    %to find the genes rescued
    temp=0; AuxotrophScore={};
    % Generate random index combinations
    pairs = combvec(1:length(singleLethalsA), 1:length(singleLethalsB))';
    numPairs = size(pairs,1);
    randOrder = randperm(numPairs);
    
    for k = 1:numPairs
        i = pairs(randOrder(k), 1);
        j = pairs(randOrder(k), 2);
        
        KoRecepient = changeRxnBounds(Community, {singleLethalsA{i}, singleLethalsB{j}},0,'b');
        
        soln= optimizeCbModel(KoRecepient);
        if soln.f>0.001
            allResults = dFBAComFSEOF(KoRecepient, options);
            KoResult = compileComFSEOF(Community, pairedModelInfo, options, allResults);
            
            if min(KoResult.abundance) > 0.105 && min(KoResult.biomass) > 0.105
                temp = temp + 1;
                AuxotrophScore{temp, 1} = {singleLethalsA{i}, singleLethalsB{j}};
                AuxotrophResults(temp) = KoResult;
                AuxotrophScore{temp, 2} = (1 - abs(diff(KoResult.abundance)));
                AuxotrophScore{temp, 3} = AuxotrophResults(temp).biomass;
                AuxotrophScore{temp, 4} = AuxotrophResults(temp).abundance;
            end
        end
        
        if temp >= Size || k> 0.9*numPairs
            break;
        end
    end
    
    % Sort the Auxotrophs according to score
    if exist('AuxotrophScore','var') && size(AuxotrophScore,1)>0
        for i=1:size(AuxotrophScore,1)
            AuxotrophScore{i,2} = AuxotrophScore{i,2}*(sum(AuxotrophResults(i).biomass)/max(cellfun(@sum,{AuxotrophResults(:).biomass})));
        end
        [~, sortIdx] = sort(cell2mat(AuxotrophScore(:,2)), 'descend');
        sortedAuxotrophScore = AuxotrophScore(sortIdx, :);
        AuxotrophResults = rmfield(AuxotrophResults, {'prdtFVAMinConc1', 'prdtFVAMinConc2', 'totPrdtFVAMinConc', ...
            'prdtFVAMaxConc1', 'prdtFVAMaxConc2', 'totPrdtFVAMaxConc', 'prdtID', 'prdtRxns', ...
            'prdtMinConc1', 'prdtMinConc2', 'totPrdtMinConc', 'prdtMaxConc1', 'prdtMaxConc2', 'totPrdtMaxConc'});
        sortedAuxotrophResults = AuxotrophResults(:, sortIdx);
    end
end
end
