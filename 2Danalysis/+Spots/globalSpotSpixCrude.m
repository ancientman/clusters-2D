function [bwSpot, nucLabelSpot] = globalSpotSpixCrude(im, labelMat, el, spotFilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called for hotspot detection only
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
        thresholdFraction = 0.9;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3; %%2.0; 
        imBlur = imgaussfilt(im, howManySigma);
    case 2% bilateral
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = imbilatfilt(im);
    case 3% tophat
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        se = strel('disk', el, 8);
        imBlur = imtophat(im,se); 
        imBlur = imadjust(imBlur);
    case 4% median
        thresholdFraction = 0.90;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = medfilt2(im);
    otherwise
        warning('No spot detect type defined, using no smoothing');
        thresholdFraction = 0.95;% quantile of data with max = absolute max and min = nuclear mean
        howManySigma = 3;
        imBlur = im;
end

se = strel('disk',2,8);


thresholdFraction = 0.85;%0.5;%0.95;
howManySigma = 1.5;%2.5;
for i = 1:nNuc
    meanNuc(i) = mean2(im(labelMat==i));
    stdNuc(i) = std(im(labelMat==i), 1, 'all');    
    thresTemp1 = prctile(im(labelMat==i & im>meanNuc(i)),...
        (thresholdFraction*100),'all');
    thresTemp2 = meanNuc(i)+stdNuc(i)*howManySigma;
    thresNuc(i)  = max([thresTemp1, thresTemp2]);
%     thresNuc(i) = thresTemp1;
    imSpotTemp = im>thresNuc(i);    
    imSpotTemp(labelMat ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');    
    imSpotTemp = imdilate(imSpotTemp,se);
    imSpotTemp = bwareaopen(imSpotTemp, 16);
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end

% for i = 1:nNuc
%     meanNuc(i) = mean2(imBlur(labelMat==i));
%     stdNuc(i) = std(imBlur(labelMat==i), 1, 'all');    
%     thresTemp1 = prctile(imBlur(labelMat==i & imBlur>meanNuc(i)),...
%         (thresholdFraction*100),'all');
%     thresTemp2 = meanNuc(i)+stdNuc(i)*howManySigma;
%     thresNuc(i)  = max([thresTemp1, thresTemp2]);
%     thresNuc(i) = thresTemp1;
%     imSpotTemp = imBlur>thresNuc(i);    
%     imSpotTemp(labelMat ~= i) = 0;
%     imSpotTemp = imfill(imSpotTemp, 'holes');    
%     imSpotTemp = imdilate(imSpotTemp,se);
%     imSpotTemp = bwareaopen(imSpotTemp, 16);
%     imSpot = imSpot + imSpotTemp;
%     imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
% end
bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end