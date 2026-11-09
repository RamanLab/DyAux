%% runs dynamic flux balance analysis for a community
% this outputs the concentrations of metabolites estimated by FBA and not FVA

function [result] = dFBAComCross(model,options)

% medium_nutrient: change in medium at each time step
% biomassarr: change in biomass at each time step
% mu_est:growth rate at each time step
% solnstat: array of solution stat (returns 1 if feasible and 0 otherwise)
% bioAdded: biomass added due to growth at each time step
% timearr: array of timeSteps
% prdtFBAConc1: product concentration for organism A in community at each time step
% prdtFBAConc2: product concentration for organism B in community at each time step
% totPrdtFBAConc: total product concentration of community at each time step
% prdtID: Product reaction numbers
% prdtRxns: Product reaction names
% modelName: Community model name

%% to initiate necessary arrays
result = struct;
minConc = 0.00001; %the minimal detectable concentration in medium
temp = 0; %iteration number through timeSteps
threshold = 0.0000001;

%define product IDs
prdtID = zeros(length(options.Products),length(options.initBiomass));
prdtRxns= cell(length(options.Products),length(options.initBiomass));

%to monitor how the medium changes over time
%array of nutrient concentrations across timeSteps
medium_nutrient = zeros(length(options.initMedium),round(options.maxTime/options.delt,0));
medium = options.initMedium; %to set initial state of medium
medium_nutrient(:,1) = medium;

%array of biomass, added biomass and growth rate across timeSteps
biomass = zeros(length(model.infoCom.spName),1);
[biomassarr,biomass_added,mu_est,CrossMets] = deal(zeros(length(options.initBiomass),round(options.maxTime/options.delt,0)));
biomassarr(:,1) = biomass;

%time and solution.stat arrays
[timearr,statarr] = deal(zeros(1,round(options.maxTime/options.delt,0)));
infeasibleCount=0; %to keep track of infeasible solutions

%initiate biomass
biomass(:) = options.initBiomass;
%to estimate concentration of products
[prdtFBAConc1,prdtFBAConc2]= deal(zeros(length(options.Products),round(options.maxTime/options.delt,0)));

for i=1:length(model.infoCom.spName)
    model.lb(nonzeros(model.indCom.EXsp(:,i)))=-1000;
    model.ub(nonzeros(model.indCom.EXsp(:,i)))=1000;
end

%% dynamic FBA
for timeStep = 1: round(options.maxTime/options.delt,0)
    tcurr = timeStep*options.delt;
    temp = temp+1; %iteration number
    timearr(temp) = tcurr;
    %exchange rxns of each community
    for i=1:length(biomass)
        mediumExcRxns{:,i} = intersect(model.indCom.EXsp(:,i),find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets))));
    end
    %estimating medium consumption
    for i = 1:length(options.mediumMets)
        uptake = computeUptakeNew(options.mediumMets(i), medium_nutrient(i), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        if temp>1
            uptake = computeUptakeNew(options.mediumMets(i), medium_nutrient(i,temp-1), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        end
        for j = 1:length(biomass)
            mediumExchange = intersect(model.indCom.EXsp(:,j),find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(i)))));
            %update model bounds with computed uptake
            if ~isempty(mediumExchange) && ~isempty(uptake) && length(uptake)==length(biomass)
                if biomass(j) > 0
                    % update fluxes are negative
                    model.lb(mediumExchange) = - uptake(j)/(biomass(j)*options.delt);
                end
            end
        end
    end
    % optional: block uptake (lb = 0) and/or all flux (lb = ub = 0) of chosen exchanges
    if isfield(options,'blockUptakeRxns') && ~isempty(options.blockUptakeRxns)
        model.lb(options.blockUptakeRxns) = 0;      % species can still secrete, not take up
    end
    if isfield(options,'knockRxns') && ~isempty(options.knockRxns)
        model.lb(options.knockRxns) = 0;            % species can neither take up
        model.ub(options.knockRxns) = 0;            % nor secrete
    end
    %standard FBA of model
    solution = optimizeCbModel(model,'max');
    %% update biomass and growth rates
    for i=1:length(biomass)
        statarr(temp) = solution.stat;
        if biomass(i)==0 && solution.stat==1
            biomass_added(i,temp) = 0;
        else
            if solution.stat ==1 && infeasibleCount<=5
                biomass_added(i,temp) = solution.x(model.indCom.spBm(i))*biomass(i)*options.delt;
                mu_est(i,temp) = biomass_added(i,temp)/(biomass(i)*options.delt);
                
                %to perform FBA for products
                for k=1:length(options.Products)
                    Flux1=[]; Flux2=[];
                    if ~isempty(find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1,1))
                        prdtID(k,i) = (find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1));
                    end
                    if i==1 && prdtID(k,i)~=0
                        prdtRxns(k,1) = model.infoCom.EXsp(prdtID(k,1),1);
                        Flux1 = solution.x(model.indCom.EXsp(prdtID(k,1),1));
                        if ~isempty(Flux1) && ~isnan(Flux1)
                            prdtFBAConc1(k,temp) = Flux1*biomass(i)*options.delt;
                        end
                    end
                    if i==2 && prdtID(k,i)~=0
                        prdtRxns(k,2) = model.infoCom.EXsp(prdtID(k,2),2);
                        Flux2 = solution.x(model.indCom.EXsp(prdtID(k,2),2));
                        if ~isempty(Flux2) && ~isnan(Flux2)
                            prdtFBAConc2(k,temp) = Flux2*biomass(i)*options.delt;
                        end
                    end
                    if ~isempty(prdtFBAConc1(k,temp)) && ~isempty(prdtFBAConc2(k,temp))
                        if (prdtFBAConc1(k,temp)>=threshold) && (prdtFBAConc2(k,temp)<=-(threshold))
                            CrossMets(k,temp) =  min(abs(prdtFBAConc1(k,temp)),abs(prdtFBAConc2(k,temp)));
                            
                        elseif ((prdtFBAConc1(k,temp)<=-(threshold)) && (prdtFBAConc2(k,temp)>=threshold))
                            CrossMets(k,temp) = -min(abs(prdtFBAConc1(k,temp)),abs(prdtFBAConc2(k,temp)));
                        else
                            CrossMets(k,temp)=0;
                        end
                    end
                end
                
            else
                biomass_added(i,temp)=0;
                infeasibleCount = infeasibleCount + 1;
                mu_est(i,temp) = NaN;
            end
        end
        %update all arrays
        biomass(i) = biomass(i) + biomass_added(i,temp);
        biomassarr(i,temp) = biomass(i);
    end
    
    %% to update the concentration of all media components
    for m =1:length(options.mediumMets)
        %only update media components for species that are actually present
        if sum(biomass)>0  && solution.stat ==1
            %in case a specific model does not have metabolite m
            metIndex = intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(m)))));
            if solution.stat == 1 && ~isempty(metIndex)
                medium(m) = medium(m) + solution.x(metIndex)*sum(biomass)*options.delt;
            end
        end
    end
    
    %% set all medium components that have fallen below a threshold value to zero, to avoid solver problems with very small values
    for m = 1: length(options.mediumMets)
        if medium(m) <= minConc
            medium(m)=0;
        end
    end
    %change in media components and biomass with time
    medium_nutrient(:,temp) = medium;
    
end

%% update results
result.medium_nutrient = medium_nutrient;
result.biomassarr = biomassarr;
result.mu_est = mu_est;  result.solnstat = statarr;
result.bioAdded = biomass_added;
result.timearr = timearr; result.CrossMets = CrossMets;
result.prdtFBAConc1 = prdtFBAConc1; result.prdtFBAConc2 = prdtFBAConc2;
result.totPrdtFBAConc = prdtFBAConc1 + prdtFBAConc2;
result.prdtID = prdtID; result.prdtRxns = prdtRxns;
result.modelName = model.description;
end
