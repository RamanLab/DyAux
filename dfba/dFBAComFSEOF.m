%% runs dynamic flux balance analysis for a community
%this gives the raw concentrations and are compiled based on abundance
%cutoff by the compileCom function
%% uses fastFVA function - use fluxVariability for other solvers
function [result] = dFBAComFSEOF(model,options)

% medium_nutrient: change in medium at each time step
% biomassarr: change in biomass at each time step
% mu_est:growth rate at each time step
% solnstat: array of solution stat (returns 1 if feasible and 0 otherwise)
% solnarr: array of all fluxes
% bioAdded: biomass added due to growth at each time step
% timearr: array of timeSteps
% prdtFVAConc1: product concentration for organism A in community at each time step
% prdtFVAConc2: product concentration for organism B in community at each time step
% totPrdtFVAConc: total product concentration of community at each time step
% prdtID: Product reaction numbers
% prdtRxns: Product reaction names
% modelName: Community model name

%% to initiate necessary arrays
result = struct; 
minConc = 0.001; %the minimal detectable concentration in medium
temp = 0; %iteration number through timeSteps

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
[biomassarr,biomass_added,mu_est] = deal(zeros(length(options.initBiomass),round(options.maxTime/options.delt,0)));
biomassarr(:,1) = biomass;  solnarr=[];

%time and solution.stat arrays
[timearr,statarr] = deal(zeros(1,round(options.maxTime/options.delt,0)));
infeasibleCount=0; %to keep track of infeasible solutions

%initiate biomass
biomass(:) = options.initBiomass;
%to estimate concentration of products
[prdtFVAMinConc1,prdtFVAMinConc2,prdtFVAMaxConc1,prdtFVAMaxConc2]= deal(zeros(length(options.Products),round(options.maxTime/options.delt,0)));

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
        uptake = computeUptake(options.mediumMets(i), medium_nutrient(i), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        if temp>1
            uptake = computeUptake(options.mediumMets(i), medium_nutrient(i,temp-1), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        end
        for j = 1:length(biomass)
            mediumExchange = intersect(model.indCom.EXsp(:,j),find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(i)))));
            %update model bounds with computed uptake
            if ~isempty(mediumExchange) && ~isempty(uptake) && length(uptake)==length(biomass)
                model.lb(mediumExchange) = - uptake(j)/(options.delt*biomass(j)); %uptake fluxes are negative
            end
        end
    end
    
    %standard FBA of model
    solution = optimizeCbModel(model,'max','one');
    %% update biomass and growth rates
    for i=1:length(biomass)
        statarr(temp) = solution.stat;
        
        if biomass(i)==0 && solution.stat==1
            biomass_added(i,temp) = 0;
            solnarr(:,temp) = 0; 
        else
            if solution.stat ==1 && infeasibleCount<=5
                biomass_added(i,temp) = solution.x(model.indCom.spBm(i))*biomass(i)*options.delt;
                mu_est(i,temp) = biomass_added(i,temp)/(biomass(i)*options.delt);
                solnarr(:,temp) = solution.x;
                
                %to perform FVA for products
                for k=1:length(options.Products)
                    maxFlux1=[]; maxFlux2=[];
                    if ~isempty(find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1,1))
                        prdtID(k,i) = (find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1));
                    end
                    if i==1 && prdtID(k,i)~=0
                        prdtRxns(k,1) = model.infoCom.EXsp(prdtID(k,1),1);
                        if options.solver == 'ibm_cplex'
                            [minFlux1, maxFlux1] = fastFVA(model, 90, 'max','ibm_cplex',prdtRxns(k,1));
                        else
                            [minFlux1, maxFlux1] = fluxVariability(model, 90, 'max',prdtRxns(k,1));
                        end
                        if ~isempty(minFlux1) && ~isnan(minFlux1)
                            prdtFVAMinConc1(k,temp) = minFlux1*biomass(i)*options.delt;               
                        end
                        if ~isempty(maxFlux1) && ~isnan(maxFlux1)
                            prdtFVAMaxConc1(k,temp) = maxFlux1*biomass(i)*options.delt;
                        end
                    end
                    if i==2 && prdtID(k,i)~=0
                        prdtRxns(k,2) = model.infoCom.EXsp(prdtID(k,2),2);
                        if options.solver == 'ibm_cplex'
                            [minFlux2, maxFlux2] = fastFVA(model, 90, 'max','ibm_cplex',prdtRxns(k,2));
                        else
                            [minFlux2, maxFlux2] = fluxVariability(model, 90, 'max',prdtRxns(k,2));
                        end
                        if ~isempty(minFlux2) && ~isnan(minFlux2)
                            prdtFVAMinConc2(k,temp) = minFlux2*biomass(i)*options.delt;
                        end
                        if ~isempty(maxFlux2) && ~isnan(maxFlux2)
                            prdtFVAMaxConc2(k,temp) = maxFlux2*biomass(i)*options.delt;
                        end
                    end
                end
                
            else
                biomass_added(i,temp)=0;
                infeasibleCount = infeasibleCount + 1;
                mu_est(i,temp) = NaN;
                solnarr(:,temp)= NaN;
            end
        end
        
        biomass(i) = biomass(i) + biomass_added(i,temp);
        biomassarr(i,temp) = biomass(i);
%         solnarr(1:size(model.rxns,1),temp) = soln_arr';
    end
    
    %% to update the concentration of all media components
    for m =1:length(options.mediumMets)
        for i=1:length(biomass)
            %only update media components for species that are actually present
            if biomass(i)>0  && solution.stat ==1
                %in case a specific model does not have metabolite m
                metIndex = intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(m)))));
                if solution.stat == 1 && ~isempty(metIndex)
                    medium(m) = medium(m) + solution.x(metIndex)*biomass(i)*options.delt;
                end
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
result.timearr = timearr; result.solnarr = solnarr;
result.prdtFVAMinConc1 = prdtFVAMinConc1; result.prdtFVAMinConc2 = prdtFVAMinConc2;
result.totPrdtFVAMinConc = prdtFVAMinConc1 + prdtFVAMinConc2;
result.prdtFVAMaxConc1 = prdtFVAMaxConc1; result.prdtFVAMaxConc2 = prdtFVAMaxConc2;
result.totPrdtFVAMaxConc = prdtFVAMaxConc1 + prdtFVAMaxConc2;
result.prdtID = prdtID; result.prdtRxns = prdtRxns;
result.modelName = model.description;
end
