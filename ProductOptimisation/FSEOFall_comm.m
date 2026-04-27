function[FseofAll]= FSEOFall_comm(model,pairedModelInfo,options)

    %% set up elimination lists
    %finding exchange and transport rxns
    excRxns = model.rxns(findExcRxns(model)==1);
    transRxns = union(excRxns,findTransRxns(model));
%     FseofAll = [];

    %elimination list for amplification targets
    eliListInc = union(transRxns,'ATPM');
    eliListInc = union(eliListInc, model.rxns(model.c==1));

    %elimination list for knockout targets
    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model);
    SL = delRxn(grRateKO < 1e-5);
    eliListDec = union(eliListInc,SL);
%     eliListDec = eliListInc;

    
    
    %% running FSEOF for all exchange metabolites
    for i = 1:length(options.Products)
        [FluxMatrix, IncFlux, DecFlux] = FSEOF_comm (model,pairedModelInfo,options.Products(i),eliListInc,eliListDec,options);
        if ~isempty(IncFlux) 
            IncAll{i,1} = IncFlux(:,4); 
        else
            IncAll{i,1} = [];
        end
        if ~isempty(DecFlux)
            DecAll{i,1} = DecFlux(:,4); 
        else 
            DecAll{i,1} = [];
        end
    end
    Fseof_temp = horzcat(options.Products,IncAll,DecAll);

    %removing empty entries
    temp=0;
    for j=1:length(options.Products)    
        if ~isempty(Fseof_temp{j,2}) || ~isempty(Fseof_temp{j,3})
            temp=temp+1;
            FseofAll(temp,:) = Fseof_temp(j,:);
        else
            FseofAll = cell(0,3);
        end
    end


end
