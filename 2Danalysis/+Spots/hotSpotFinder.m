function [hotSpotStruct] = hotSpotFinder(im, nucPropStruct, thres, nucEl, spotEl, nucleusDetectMode, minNucSize, minSpotSize, clearBorder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This technique detects spots based on the fraction of the maximum
%   brightness per nucleus label. So be cautioned that in the nuclei where
%   the spots are not exceptionally bright there might be false positives,
%   as there is no cutoff for deviation from the mean background (nuclear)
%   intensity mean.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imScale = rescale(im);
% sx = nucEl; sy = nucEl;
% nucleusDetectMode = 1;
BW = Nucleus.nucDetect(im, thres, nucleusDetectMode, nucEl, minNucSize, clearBorder, 0, 1);
labelMat = bwlabel(BW);
nucLabelStats = regionprops(labelMat, 'Area', 'Centroid');
[newGlobalNucLabelMat] = Spots.globalHotSpotLabelCorrect(labelMat, nucLabelStats, nucPropStruct);
labelMat = newGlobalNucLabelMat; % reassign the labels after checking if consistent with the raw frames
nucLabel = nonzeros(unique(labelMat, 'sorted'));
nucLabelStats = regionprops(labelMat, 'Area', 'Centroid');
spotDetect = 4; % Mode #4 is dedicated for hotspot detection
spotFilter = 1; % Gaussian
minSpotSize = 2*minSpotSize; % Double the size of individual spots
% Segment global hotspots
[bwHotSpot, hotSpotLabelMatNL] = Spots.spotDetect(imScale, labelMat, spotEl, spotDetect, spotFilter, minSpotSize);
% sort nuclei with detected spots
uniqSpotLabelNL = nonzeros(unique(hotSpotLabelMatNL, 'sorted'));
noSpotNuc = setdiff(nucLabel, uniqSpotLabelNL);
if ~isempty(noSpotNuc)
    uniqSpotLabelNL = HelperFunctions.vectorInsertAfter(uniqSpotLabelNL, 0, noSpotNuc);
    uniqSpotLabelNL = sort(uniqSpotLabelNL);
end
% compute the properties of individual spots, and asign them to cell array,
% with each array element representing one nucleus
spotCentroidUniq = cell(1, length(uniqSpotLabelNL));
spotAreaUniq = cell(1, length(uniqSpotLabelNL));
spotIdxListUniq = cell(1, length(uniqSpotLabelNL));
spotNLtemp = hotSpotLabelMatNL;
for i=1:length(uniqSpotLabelNL)
    spotNLtemp(spotNLtemp~=uniqSpotLabelNL(i)) = 0;
    spotNLtemp(spotNLtemp==uniqSpotLabelNL(i)) = 1;
    CC = bwconncomp(spotNLtemp);
    stats = regionprops(CC, 'Centroid', 'Area', 'PixelIdxList');
    spotCentroidUniq{i} = cat(1,stats.Centroid);
    spotAreaUniq{i} = cat(1,stats.Area);
    spotIdxListUniq{i} = cat(1, stats.PixelIdxList);
    spotNLtemp = hotSpotLabelMatNL;
end
hotspotAreaNL = regionprops(hotSpotLabelMatNL,'Area');
hotSpotStruct.globalNucLabel = labelMat;
hotSpotStruct.globalNucAreaNL = cat(1, nucLabelStats.Area);
hotSpotStruct.globalNucCentroidNL = cat(1, nucLabelStats.Centroid);
hotSpotStruct.hotspotBW = bwHotSpot;
hotSpotStruct.hotspotNL = hotSpotLabelMatNL;
hotSpotStruct.hotspotAreaNL = cell2mat(struct2cell(hotspotAreaNL));
hotSpotStruct.hotspotCentroidUniq = spotCentroidUniq;
hotSpotStruct.hotspotAreaUniq = spotAreaUniq;
hotSpotStruct.hotspotIdxListUniq = spotIdxListUniq;
hotSpotStruct.noSpotNuc = noSpotNuc;
end
