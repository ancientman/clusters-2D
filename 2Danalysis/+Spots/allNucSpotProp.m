function [allSpots] = allNucSpotProp(hotSpotStruct, spotHotSpotStruct, XpixelSize, timeCount)
areaAllSpotIn = (cellfun(@(s)s.areaInHotspotNL,spotHotSpotStruct,'uni',0));
[numAllSpotIn, ~] = cellfun(@size,areaAllSpotIn,'uni',false);
numAllSpotIn = vertcat(numAllSpotIn{:});
hotSpotAreaNL = (cellfun(@sum,hotSpotStruct.hotspotAreaUniq,'uni',0));
hotSpotAreaNL = vertcat(hotSpotAreaNL{:});
hotSpotAreaNL = hotSpotAreaNL.*(XpixelSize)^2;
% spotsPerHotSpotAreaIn = zeros(size(hotSpotAreaNL));
if ~isempty(nonzeros(numAllSpotIn))
    spotsPerHotSpotAreaIn = nonzeros(numAllSpotIn)./nonzeros(hotSpotAreaNL);
else
    spotsPerHotSpotAreaIn = 0;
end
spotsPerHotSpotAreaIn(hotSpotAreaNL==0) = 0;
spotRateIn = spotsPerHotSpotAreaIn./timeCount;
areaAllSpotIn = vertcat(areaAllSpotIn{:});
areaAllSpotOut = (cellfun(@(s)s.areaOutHotspotNL,spotHotSpotStruct,'uni',0));
[numAllSpotOut, ~] = cellfun(@size,areaAllSpotOut,'uni',false);
numAllSpotOut = vertcat(numAllSpotOut{:});
notHotSpotAreaNL = nonzeros(hotSpotStruct.globalNucAreaNL) - (hotSpotAreaNL);
notHotSpotAreaNL = notHotSpotAreaNL.*(XpixelSize)^2;
% spotsPerHotSpotAreaOut = zeros(size(hotSpotAreaNL));
try
    spotsPerHotSpotAreaOut = nonzeros(numAllSpotOut)./nonzeros(notHotSpotAreaNL);
catch
    aaa = 0;
end
    
spotsPerHotSpotAreaOut(notHotSpotAreaNL==0) = 0;
spotRateOut = spotsPerHotSpotAreaOut./timeCount;
areaAllSpotOut = vertcat(areaAllSpotOut{:});
allSpots.numAllSpotIn = numAllSpotIn;
allSpots.areaAllSpotIn = areaAllSpotIn;
allSpots.spotRateIn = spotRateIn;
allSpots.numAllSpotOut = numAllSpotOut;
allSpots.areaAllSpotOut = areaAllSpotOut;
allSpots.spotRateOut = spotRateOut;
end