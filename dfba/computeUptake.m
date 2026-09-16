%% computes the substrate uptake bound based on kinetic parameters and
%availability limits

function [actualImportBound] = computeUptake(metid, conc, excRxns, model, biomass, Vmax, Km, timeSteps)

% Initialize arrays
[modelValidity,max_transport,max_availability] = deal(zeros(1,length(biomass)));

% Check if organism has minimal biomass and valid exchange reaction
for i = 1:length(biomass)
    if any(contains(model.rxns(excRxns{:,i}), erase(metid,'[u]'))) && biomass(i)>0
        modelValidity(i)=1;
    else
        modelValidity(i)=0;
    end
end

% Compute transport limits 
for i = 1:length(biomass)
    if modelValidity(i)==1
        max_transport(i)= Vmax*(conc/(Km + conc));
        max_transport(i)= max_transport(i)*biomass(i)*timeSteps;
    end
end

% Compute total biomass 
sumbiomass = sum(biomass(modelValidity==1));

% Compute maximal availability 
for i = 1:length(biomass)
    if modelValidity(i)==1
        % maximal availability per gram of total biomass and unit time
        max_availability(i) = conc/(sumbiomass*timeSteps);
        % convert to absolute mmoles for this strain
        max_availability(i) = max_availability(i)*biomass(i)*timeSteps;
    end
end

% Set actual upper bound for FBA
actualImportBound = zeros(1,length(biomass));
for i = 1:length(biomass)
    if modelValidity(i)==1
        actualImportBound(i) = min(max_transport(i), max_availability(i));
    end
end

%(the actual amount imported is determined by FBA)
%if the sum of the actual import bounds exceeds the total nutrient concentration, then
%the total amount of nutrient will be imported

% Normalize if sum exceeds concentration
max_sum_import = sum(actualImportBound(modelValidity==1));
if max_sum_import > conc
    for i = 1:length(biomass)
        if modelValidity(i)==1
            actualImportBound(i) = actualImportBound(i)*(conc/max_sum_import);
        end
    end
end
end
