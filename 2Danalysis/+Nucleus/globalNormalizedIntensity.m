function [imNormAdd, imSpotAdd, imSpotMaskAdd] = globalNormalizedIntensity(nucPropStruct, im, imSpot, tFinal)
nNuc = max(nucPropStruct{1}.labelMat, [], 'all');
im = double(im);
imSpotAdd = double(zeros(size(im, 1),size(im, 2)));
imSpotMaskAdd = double(zeros(size(im, 1),size(im, 2)));
for i=1:tFinal
    imSpot(:,:,i) = rescale(imSpot(:,:,i));
    imSpotAdd = imSpotAdd+imSpot(:,:,i);
    imSpotMaskAdd = imSpotMaskAdd+imbinarize(imSpot(:,:,i));
end
imAdd = double(zeros(size(im, 1),size(im, 2),tFinal));
bgMean = double(zeros(1,tFinal));
imNormAdd = double(zeros(size(im, 1),size(im, 2)));
for t = 1:tFinal
    imTemp = double(im(:,:,t));
    bgMean(t) = (mean(imTemp(nucPropStruct{t}.labelMat==0)));
    imTemp = double(im(:,:,t))-double(bgMean(t));
    for i = 1:nNuc
        imTemp(nucPropStruct{t}.labelMat~=i) = 0;
        imTemp = rescale(imTemp);
        imAdd(:,:,t) = imadd(double(imAdd(:,:,t)),imTemp);
        imTemp = double(im(:,:,t))-double(bgMean(t));
    end
    imNormAdd = double(imNormAdd)+ double(imAdd(:,:,t));    
end
imNormAdd = double(imNormAdd)/tFinal;
end