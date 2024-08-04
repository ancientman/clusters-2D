function [newGlobalNucLabelMat] = globalHotSpotLabelCorrect(globalNucLabelMat, nucLabelStats, nucPropStruct)
nNuc = max(nucPropStruct{1}.labelMat, [], 'all');
nucRadius = sqrt(nucPropStruct{1}.area./pi);
globalNucLabel = nonzeros(unique(globalNucLabelMat));
newGlobalNucLabelMat = zeros(size(globalNucLabelMat));
for i=1:length(globalNucLabel)
    for k=1:nNuc
        if(pdist([nucPropStruct{1}.centroid(k,1), nucPropStruct{1}.centroid(k,2);...
                nucLabelStats(i).Centroid(1), nucLabelStats(i).Centroid(2)], 'euclidean')...
                <nucRadius)
            newGlobalNucLabelMat(globalNucLabelMat==i) = k;
        end
    end
end
end