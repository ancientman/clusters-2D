function [fim] = imNormAv(nucPropStruct, im, tStart, tEnd)
nNuc = max(nucPropStruct{tStart}.labelMat, [], 'all');
im = double(im);
imNorm = double(zeros(size(im, 1),size(im, 2),(tEnd-tStart+1)));
bgMean = double(zeros(1,(tEnd-tStart+1)));
imNormAv = double(zeros(size(im, 1),size(im, 2)));
for t = tStart:tEnd
    tTemp = t-tStart+1;
    imTemp = double(im(:,:,t));
    bgMean(t) = (mean(imTemp(nucPropStruct{t}.labelMat==0)));
    imTemp = double(im(:,:,t))-double(bgMean(t));
    for i = 1:nNuc
        imTemp(nucPropStruct{t}.labelMat~=i) = 0;
        imTemp = rescale(imTemp);
        imNorm(:,:,tTemp) = imadd(double(imNorm(:,:,tTemp)),imTemp);
        imTemp = double(im(:,:,t))-double(bgMean(t)); % beware of the background values
    end
    if(t>1)
%         imNorm(:,:,t) = imregister(imNorm(:,:,1), imNorm(:,:,t), 'rigid', 'RegularStepGradientDescent', 'MeanSquares');
    end
    imNormAv = double(imNormAv)+ double(imNorm(:,:,tTemp));    
end
imNormAv = double(imNormAv)/(tEnd-tStart+1);
[fim] = imNormAv;
end