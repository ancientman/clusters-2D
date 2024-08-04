function [spotPropStruct] = spotProp(I, nucSpotLabel, nNuc, t, stopMark)
frameTemp = zeros(size(nucSpotLabel));
CC = struct([]);
SC = struct([]);
for i = 1:nNuc
    if ~ismember(i, nucSpotLabel)
        SC{i}.meanSpotRelIntensity = 0;
        SC{i}.sumSpotAbsIntensity = 0;
        SC{i}.meanSpotAbsIntensity = 0;
        SC{i}.stDevSpotIntensity = 0;
        SC{i}.spotMols = 0;
        SC{i}.spotArea = 0;
        SC{i}.spotCentroid = [0 0];        
    end
    frameTemp(nucSpotLabel == i) = 1;
    CC{i} = bwconncomp(frameTemp);
    if CC{i}.NumObjects ~=0
        for k = 1:CC{i}.NumObjects
            SC{i}.meanSpotRelIntensity(k) = mean(rescale(I(CC{i}.PixelIdxList{k})));
            SC{i}.sumSpotAbsIntensity(k) = sum((I(CC{i}.PixelIdxList{k})));
            SC{i}.meanSpotAbsIntensity(k) = mean((I(CC{i}.PixelIdxList{k})));
            SC{i}.stDevSpotIntensity(k) = std(rescale(I(CC{i}.PixelIdxList{k})));       
        end
        SC{i}.spotMols = round(SC{i}.sumSpotAbsIntensity./min(SC{i}.sumSpotAbsIntensity));
        A = regionprops(CC{i}, 'area');
        SC{i}.spotArea = cat(1, A.Area);
        S = regionprops(CC{i}, 'centroid');
        SC{i}.spotCentroid = cat(1, S.Centroid);
        Idx = regionprops(CC{i}, 'PixelIdxList');
        SC{i}.spotIdx = struct2cell(Idx);%cat(1, Idx.PixelIdxList);
    else
        SC{i}.meanSpotRelIntensity(1) = 0;
        SC{i}.sumSpotAbsIntensity(1) = 0;
        SC{i}.meanSpotAbsIntensity(1) = 0;
        SC{i}.stDevSpotIntensity(1) = 0;
        SC{i}.spotMols = 0;
        SC{i}.spotArea = 0;        
        SC{i}.spotCentroid = [];
        SC{i}.spotIdx = {[]};
    end
    frameTemp = zeros(size(nucSpotLabel));
end
spotPropStruct = SC;
if t == stopMark
    fprintf("in nucProp");
end
end