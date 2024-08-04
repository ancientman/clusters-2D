function [bwSpot, nucLabelSpot] = spotPixCrude3(im, labelMat, el, spotFilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called if "deviationofGradient" is used
%   Use tophat with extreme caution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imScale = rescale(im);
nNuc = max(labelMat,[],'all');
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
meanNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc, 1);

thresholdFraction = 0.95;

switch spotFilter
    case 1% gaussian
        imBlur = imgaussfilt(im);
    case 2% bilateral
        imBlur = imbilatfilt(im);
    case 3% tophat
        se = strel('disk', el, 8);
        imBlur = imtophat(imScale,se);
%         imBlur = imadjust(imBlur);
    case 4% median
        warning('No spot detect type defined, using median');
        imBlur = medfilt2(im);
    otherwise
        warning('No spot detect type defined, using median');
        imBlur = im;
end
imAdd = imgradient(imScale, 'sobel') + medfilt2((imBlur));
for i = 1:nNuc
    meanNuc(i) = mean2(imAdd(labelMat==i));
    thresNuc(i) = prctile(imAdd(labelMat==i & imAdd>meanNuc(i)),...
        (thresholdFraction*100),'all');    
    imSpotTemp = imAdd>thresNuc(i);    
    imSpotTemp(labelMat ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end
bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end


