%% computes the substrate uptake bound based on kinetic parameters and
%availability limits

function [actualImportBound] = computeUptake(metid, conc, excRxns, model, biomass, Vmax, Km, timeSteps)
[modelValidity,max_transport,max_availability] = deal(zeros(1,length(biomass)));

%checking if organism has minimal biomass
for i = 1:length(biomass)
    if any(contains(model.rxns(excRxns{:,i}),erase(metid,'[u]'))==1) && biomass(i)>0
        modelValidity(i)=1;
    else
        modelValidity(i) =0;
    end
end

for i = 1:length(biomass)
    if modelValidity(i)==1
        max_transport(i)= Vmax*(conc/(Km + conc));
        %total mmoles of nutrient that can be imported by the amount of present biomass
        max_transport(i)= max_transport(i)*biomass(i)*timeSteps;
    end
end

% maximal amount of nutrient available to each strain
sumbiomass=0;   % weighted sum of biomasses for all models that consumed the nutrient

for i = 1:length(biomass)
    if modelValidity(i)==1
        sumbiomass=sumbiomass+biomass(i);
        max_availability(i)=conc/(sumbiomass*timeSteps);   %this is the maximal availability per gram of total biomass and unit time
    end
end

%calculate the maximal availability in terms of absolute mmoles of nutrient
for i = 1:length(biomass)
    if modelValidity(i)==1
        max_availability(i)=max_availability(i)*biomass(i)*timeSteps;
    end
end

%set the actual upper bound for FBA
actualImportBound=[];
%calculate the minimum of the transport limit and the availability limit,
for i = 1:length(biomass)
    if modelValidity(i)==1
        actualImportBound(i)= min(max_transport(i), max_availability(i));
    end
end

%(the actual amount imported is determined by FBA)
%if the sum of the actual import bounds exceeds the total nutrient concentration, then
%the total amount of nutrient will be imported

max_sum_import=0;

for i = 1:length(biomass)
    if modelValidity(i)==1
        max_sum_import= max_sum_import + actualImportBound(i);
    end
end
if max_sum_import>conc
    for i = 1:length(biomass)
        if modelValidity(i)==1
            actualImportBound(i) = actualImportBound(i)*(conc/max_sum_import);
        end
    end
end
end
