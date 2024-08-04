function [bwSpot, nucLabelSpot] = spotPixCrude2(im, labelMat, el, spotFilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Gets called if "deviationFromMax" is used
%   This technique detects spots based on the fraction of the maximum
%   brightness per nucleus label. So be cautioned that in the nuclei where
%   the spots are not exceptionally bright there might be falso positives,
%   as there is no cutoff for deviation from the mean background (nuclear)
%   intensity mean.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nNuc = max(labelMat,[],'all');
thresNuc = zeros(nNuc, 1);
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
switch spotFilter
    case 1% gaussian
        thresholdFraction = 0.5;
        imBlur = imgaussfilt(im, 5);
    case 2% bilateral
        thresholdFraction = 0.30;
        imBlur = imbilatfilt(im);
    case 3% tophat
        thresholdFraction = 0.9;
        se = strel('disk', el, 8);
        imBlur = imtophat(im,se); 
    case 4% median
        thresholdFraction = 0.30;
        imBlur = medfilt2(im);
    otherwise
        warning('No spot detect type defined, using median');
        thresholdFraction = 0.30;
          imBlur = medfilt2(im);
end

for i=1:nNuc
    thresNuc(i) = thresholdFraction*prctile(imBlur(labelMat==i), (99),'all'); 
    imSpotTemp = imBlur>thresNuc(i);
    imSpotTemp(labelMat ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end  

bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end