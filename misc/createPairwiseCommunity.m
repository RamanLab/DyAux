%% creates pairwise community models from a list of organisms
function [Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList)

%to alter model bounds from +/- 99999 to +/- 1000
for i = 1: length(speciesList)
    speciesList{i} = alterbounds(speciesList{i});
end

Communities ={};
pairedModelInfo ={}; %contains organism and biomass rxn lists
count=0;

for i = 1:length(speciesList)
    for j = i+1 :length(speciesList)
        model1 = speciesList{i}; model2 = speciesList{j};
        models = {model1,model2};
        options.spBm = {biomassList{i}, biomassList{j}};
        options.spAbbr = {spNameList{i},spNameList{j}};
        options.metExId = ('_e');
        options.sepUtEx = false;
        count = count+1;
        
        ModelCom = createCommModels(models,options);
        Communities{count} = ModelCom;
        
        % information on joined models
        pairedModelInfo{count, 1} = strcat(spNameList{i}, '_', spNameList{j});
        pairedModelInfo{count, 2} = spNameList{i};
        pairedModelInfo{count, 3} = biomassList{i};
        pairedModelInfo{count, 4} = spNameList{j};
        pairedModelInfo{count, 5} = biomassList{j};
    end
end
end
