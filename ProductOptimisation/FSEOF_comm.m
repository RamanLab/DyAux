%FSEOF
function [FluxMatrix, IncFlux, DecFlux] = FSEOF_comm (model,pairedModelInfo,TargetRxn,eliListInc,eliListDec,options)

%Maximum biomass
allResults = dFBAComFSEOF(model,options); %%%
CommunityResults = compileComFSEOF(model,pairedModelInfo,options,allResults);
Vbiomass = CommunityResults.biomass;

%% Initial product and maximum product
Viniprdt = CommunityResults.totPrdtMaxConc;

%to change objective of model
modelPrdt=model;
TargetRxn1 = strcat(pairedModelInfo{1,2},'_',TargetRxn);
TargetRxn2 = strcat(pairedModelInfo{1,4},'_',TargetRxn);
modelPrdt.infoCom.spBm(:,:) = [TargetRxn1,TargetRxn2];
modelPrdt.indCom.spBm(:,:) = [find(strcmp(model.rxns,TargetRxn1)),find(strcmp(model.rxns,TargetRxn2))];
modelPrdt.c(:)=0;
modelPrdt.c(modelPrdt.indCom.spBm(:,:))=1;
FBAsolnPrdt = optimizeCbModel(modelPrdt,'max','one');

TargetRxnComm = modelPrdt.infoCom.EXcom(find(contains(modelPrdt.infoCom.EXcom,TargetRxn)==1));
TargetRxnCommID = modelPrdt.indCom.EXcom(find(contains(modelPrdt.infoCom.EXcom,TargetRxn)==1));
if FBAsolnPrdt.stat == 1
    if FBAsolnPrdt.x(TargetRxnCommID)>1e-4
        Vmaxprdt = FBAsolnPrdt.x(TargetRxnCommID);
    else
        Vmaxprdt = 0;
        fprintf('metabolite cannot be produced, only gets consumed in');
        disp(TargetRxn);
    end
else
    disp('Infeasible solution')
    Vmaxprdt = 0;
end


%% Setup of flux matrix
steps = randperm(10);
FluxMatrix = zeros(length(model.rxns),length(steps)+1);
if size(CommunityResults.solnarr,1) == size(FluxMatrix,1)
    FluxMatrix(:,1) = CommunityResults.soln_arr;
else
    FluxMatrix(:,1) = NaN;
end
IncFlux = [];
DecFlux = [];

%% Enforcing flux
if Vmaxprdt>=1e-4
    for j =1:length(steps)
        Venfprdt = (steps(j)/10)*Vmaxprdt;
        model_mut = changeRxnBounds(model, TargetRxnComm, Venfprdt-1e-3, 'l');
        model_mut = changeRxnBounds(model_mut, TargetRxnComm, Venfprdt+1e-3, 'u');
        MutResultRaw = dFBAComFSEOF(model_mut,options);
        MutResults = compileComFSEOF(model_mut,pairedModelInfo,options,MutResultRaw);
        if ~isempty(MutResults) && size(MutResults.soln_arr,1) == size(FluxMatrix,1)
            FluxMatrix(:,steps(j)+1) = MutResults.soln_arr;
        else
            FluxMatrix(:,steps(j)+1) = NaN;
            fprintf('%02d: infeasible solution at Venfprdt = %f \n\n',steps(j),Venfprdt);
        end
    end
    
    %% identification of amplification and knockout targets
    % for RxnIndex = 1:length(model.rxns)
    %% identification of amplification and knockout targets
    for RxnIndex = 1:length(model.rxns)
       if FluxMatrix (RxnIndex, 1) < FluxMatrix (RxnIndex, 3) && FluxMatrix (RxnIndex, 4) <= FluxMatrix (RxnIndex, 6) && ~ismember(model.rxns(RxnIndex),eliListInc)
            NewRowI = FluxMatrix(RxnIndex,:);
            Score = 1/((Vmaxprdt-Viniprdt)/abs(abs(FluxMatrix (RxnIndex, 9)) - abs(FluxMatrix (RxnIndex, 2))));
            NewRowI = [Score, NewRowI , RxnIndex, model.rxns(RxnIndex)];
            IncFlux = [IncFlux; NewRowI];
       elseif FluxMatrix (RxnIndex, 1) > FluxMatrix (RxnIndex, 3) && FluxMatrix (RxnIndex, 4) >= FluxMatrix (RxnIndex, 6) && ~ismember(model.rxns(RxnIndex),eliListDec)
            NewRowII = FluxMatrix(RxnIndex,:);
            Score = 1/((Vmaxprdt-Viniprdt)/abs(abs(FluxMatrix (RxnIndex, 9)) - abs(FluxMatrix (RxnIndex, 2))));
            NewRowII = [Score, NewRowII , RxnIndex, model.rxns(RxnIndex)];
            DecFlux = [DecFlux ; NewRowII];         
       end
    end
    
    if ~isempty(IncFlux)
        IncFlux = sortrows(IncFlux, 1,'descend');
    end
    if ~isempty(DecFlux)
        DecFlux = sortrows(DecFlux, 1,'descend');
    end
    %Score in first column ; rows arranged according to inc order of scores
    %Flux ranging across steps 0(biomass obj) to 10(product obj) in columns 2-12
    %RxnNo in column 13
end
end
