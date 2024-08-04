function [hotSpotPropStruct] = hotSpotProp(hotSpotStruct)
nNuc = size(nonzeros(unique(hotSpotStruct.globalNucLabel, 'sorted')));
[spotsPerNuc, ~] = cellfun(@size,hotSpotStruct.hotspotAreaUniq,'uni',false);
hotSpotsPerNuc = vertcat(spotsPerNuc);
meanHotSpotArea = cellfun(@mean,hotSpotStruct.hotspotAreaUniq,'uni',false);
meanHotSpotArea = vertcat(meanHotSpotArea);
allHotSpotAreas = vertcat(hotSpotStruct.hotspotAreaUniq{:});
hotSpotPropStruct.hotSpotsPerNuc = hotSpotsPerNuc;
hotSpotPropStruct.meanHotSpotArea = meanHotSpotArea;
hotSpotPropStruct.allHotSpotAreas = allHotSpotAreas;
end