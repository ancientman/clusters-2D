function [spotMask, spotLabel] = spotDetect(im, labelMat, el, spotDetect, spotFilter, minSpotSize)
% mask = ones(size(im));
% mask(labelMat==0) = 0;
% If = uint16(medfilt2(im, [el, el]));
% If(mask==0) = 0;
if spotDetect == 1 % deviationFromMean
    [bwSpot, ~] = Spots.spotPixCrude(im, labelMat, el, spotFilter); %based on deviations from base
elseif spotDetect == 2 % deviationFromMax
    [bwSpot, ~] = Spots.spotPixCrude2(im, labelMat, el, spotFilter); %based on threshold from spot max
elseif spotDetect == 3 % deviationofGradient
    [bwSpot, ~] = Spots.spotPixCrude3(im, labelMat, el, spotFilter); %based on image gradient
elseif spotDetect == 4 % only for global hotspots
    [bwSpot, ~] = Spots.globalSpotSpixCrude(im, labelMat, el, spotFilter); %based on image gradient
elseif spotDetect == 5 % for adaptive threshold
    [bwSpot, labelSpot] = Spots.spotPixAdapThres(im, labelMat); %based on image gradient
else %deviationFromMean
    [bwSpot, ~] = Spots.spotPixCrude(labelMat, im, el, spotFilter); %based on deviations from base
end

if spotDetect == 5
    spotMask = bwareaopen(bwSpot, 2);
    spotLabel = labelSpot;
else

    if spotDetect == 4
        minSpotSize = 3*minSpotSize;
    end
    [maskSpot, labelSpot] = Spots.watershedSpot(im, labelMat, bwSpot, minSpotSize, spotFilter);
    % [maskSpot, labelSpot] = Spots.noWatershedSpot(labelMat, bwSpot, minSpotSize);
    spotMask = maskSpot;
    spotLabel = labelSpot;
end
end
