function [bwSpot, nucLabelSpot] = spotPixAdapThres(im, labelMat)
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

for i = 1:nNuc
    imSpotTemp = im;    
    imSpotTemp(labelMat ~= i) = 0;
%     imSpotTemp = imbinarize(imSpotTemp);
    imSpotTemp(imSpotTemp>0) = 1;
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end
bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end