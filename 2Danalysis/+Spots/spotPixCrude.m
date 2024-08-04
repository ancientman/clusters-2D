function [bwSpot, nucLabelSpot] = spotPixCrude(im, labelMat, el, spotFilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called if "deviationFromMean" is used
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = rescale(im);
nNuc = max(labelMat, [], 'all');
meanNuc = zeros(nNuc,1);
stdNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc,1);
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
switch spotFilter
    case 1% gaussian
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = imgaussfilt(im);
    case 2% bilateral
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = imbilatfilt(im);
    case 3% tophat
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        se = strel('disk', el, 8);
        imBlur = imtophat(im,se); 
%         imBlur = imadjust(imBlur);
    case 4% median
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = medfilt2(im);
    otherwise
        warning('No spot detect type defined, using median');
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = im;
end

for i = 1:nNuc
    meanNuc(i) = mean2(imBlur(labelMat==i));
    stdNuc(i) = std(imBlur(labelMat==i), 1, 'all');    
    thresTemp1 = prctile(imBlur(labelMat==i & imBlur>meanNuc(i)),...
        (thresholdFraction*100),'all');
    thresTemp2 = meanNuc(i)+stdNuc(i)*howManySigma;
    thresNuc(i)  = max([thresTemp1, thresTemp2]);
%     thresNuc(i) = thresTemp1;
    imSpotTemp = imBlur>thresNuc(i);    
    imSpotTemp(labelMat ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end
bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end