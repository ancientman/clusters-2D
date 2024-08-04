function [alignedNucProp] = alignedNucProp(CCcur, structPre, I, t, stopMark)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compares the centroid of the labels found in the prior iteration with the
%centroid in the current iteration.
%
%Takes the label corrected structure from the previous timestep and updates
%the current timestep, if there is a label mismatch. 
%
%Only works if the total labels in the previous timestep is equal to the
%current. 
%
%Note that the data structure of the two timepoints are different. The
%Current frame data structure is diretly obtained from regionprops,
%while the previous step data structure is constructed (corrected).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nNucPre = max(structPre.labelMat,[],'all');
nNucCur = CCcur.NumObjects;
if(nNucCur~=nNucPre)
    error("unequal number of nuclei");
end
labelCur = labelmatrix(CCcur);
meanCur = structPre.meanRelIntensity;
Scur = regionprops(CCcur,'centroid');%current iteration centroid
centroidCur = cat(1,Scur.Centroid);
Acur = regionprops(CCcur,'area');
areaCur = cat(1,Acur.Area);
labelPre = structPre.labelMat;
meanPre = structPre.meanRelIntensity;
centroidPre = structPre.centroid;
areaPre = structPre.area;
newCentroidCur = zeros(nNucPre,2);
newLabelMat = zeros(size(labelPre));% New label martix
newLabel = zeros(nNucPre,1);% One dimensional array storing the indices of the new labels
distCent = zeros(nNucPre,1);
newPixelIdxList = cell(nNucPre,1);% Cell array of pixel info for he labels
meanRelNucIntensity = zeros(nNucPre,1);
meanAbsNucIntensity = zeros(nNucPre,1);
stDevNucIntensity = zeros(nNucPre,1);
maxNucIntensity = zeros(nNucPre, 1);
newArea = zeros(nNucPre,1);
for i = 1:nNucPre
    for j = 1:nNucPre% Calculate distance of the current centroid with all centroids from the provious frame
        centCombine = vertcat(centroidCur(i,:),centroidPre(j,:));
        distCent(j) = pdist(centCombine,'Euclidean');% get the distances of the current centroid (i) from all centroids 
    end
    [~, idx] = min(distCent); %idx gives the centroid of which j is closest to the centroid of i
    newLabel(i) = idx;% assign the index to the new array
    newPixelIdxList(i) = CCcur.PixelIdxList(idx);
    newCentroidCur(i,1) = centroidCur(idx,1);
    newCentroidCur(i,2) = centroidCur(idx,2);
    newLabelMat((labelCur(:)==i)) = idx;
    meanRelNucIntensity(i) = mean(rescale(I(newPixelIdxList{i})));
    stDevNucIntensity(i) = std(rescale(I(newPixelIdxList{i})));
    meanAbsNucIntensity(i) = mean((I(newPixelIdxList{i})));
    maxNucIntensity(i) = max(I(newPixelIdxList{i}));
    newArea(i) = areaCur(idx); 
end
alignedNucProp.labelMat = newLabelMat;
alignedNucProp.meanRelIntensity = meanRelNucIntensity;
alignedNucProp.meanAbsIntensity = meanAbsNucIntensity;
alignedNucProp.stdIntensity = stDevNucIntensity;
alignedNucProp.maxIntensity = maxNucIntensity;
alignedNucProp.centroid = newCentroidCur;
alignedNucProp.area = newArea;
if t == stopMark
    fprintf("in nucProp");
end
end