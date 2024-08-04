function [maskSpot, labelSpot] = noWatershedSpot(labelMat, bwSpot, minSpotSize)
maskSpot = bwareaopen(bwSpot,minSpotSize);
labelSpot = labelMat.*maskSpot;
end
