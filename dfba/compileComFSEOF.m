%% Compile the community results for a single community
%% gives biomass and abundance without using the abundance cutoff
function [result] = compileComFSEOF(Community,pairedModelInfo,options,result)

% mean of solutions across a timeframe
noOfPoints = 5; %time interval across which growth and production rates are averaged

%Compile each community result
biomass = [0;0]; abundance =[0;0]; cs_cons =[];
%Products to be analyzed
[prdtMinConc1, prdtMinConc2, totPrdtMinConc] = deal(zeros(length(options.Products),1));
[prdtMaxConc1, prdtMaxConc2, totPrdtMaxConc, soln_arr] = deal(zeros(length(options.Products),1));

%Analyze data for a set of datapoints at stationary phase
startIndex = max(find(result.solnstat==1))-noOfPoints;
if startIndex > 0 %to ensure enough datapoints
    index = startIndex + find(result.solnstat(startIndex:startIndex-1+noOfPoints)==1);
    
    %mean biomass and abundance across datapoints
    biomass = mean(result.biomassarr(:,index),2);
    abundance = biomass/sum(biomass,1);
    
    csID = find(ismember(options.mediumMets,options.carbonSource));
    cs_cons = mean(options.initMedium(find(ismember(options.mediumMets,options.carbonSource)),1)-result.medium_nutrient(csID,index));
    
    %Product concentration
    prdtMinConc1 = mean(result.prdtFVAMinConc1(:,index),2);
    prdtMinConc2 = mean(result.prdtFVAMinConc2(:,index),2);
    prdtMaxConc1 = mean(result.prdtFVAMaxConc1(:,index),2);
    prdtMaxConc2 = mean(result.prdtFVAMaxConc2(:,index),2);
    totPrdtMinConc = mean(result.totPrdtFVAMinConc(:,index),2);
    totPrdtMaxConc = mean(result.totPrdtFVAMaxConc(:,index),2);
    soln_arr = mean(result.solnarr(:,index),2);
end

result.biomass = biomass;
result.abundance = abundance;
result.cs_cons = cs_cons;
result.soln_arr = soln_arr;

result.prdtMinConc1 = prdtMinConc1;
result.prdtMinConc2 = prdtMinConc2;
result.totPrdtMinConc = totPrdtMinConc;

result.prdtMaxConc1 = prdtMaxConc1;
result.prdtMaxConc2 = prdtMaxConc2;
result.totPrdtMaxConc = totPrdtMaxConc;
end
