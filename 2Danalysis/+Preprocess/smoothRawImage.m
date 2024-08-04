function [fim] = smoothRawImage(im, filterType, filterParam)
imTemp = im;
sx = filterParam;
sy = filterParam;
switch filterType
    case 1  %   gaussian
    imSmooth = HelperFunctions.gaussianFilter2D(imTemp,sx,sy);    
    case 2  %   wiener
        imSmooth = HelperFunctions.wienerFilter2D(imTemp, sx, sy);
    case 3  %   mean
        imSmooth = HelperFunctions.meanFilter2D(imTemp, sx, sy);
    case 4  %   median
        imSmooth = HelperFunctions.medianFilter2D(imTemp, sx, sy);
    case 5  %   bilateral
        imSmooth = HelperFunctions.bilateralFilter2D(imTemp, sx, sy);
    case 6  %   nonlocal
        imSmooth = HelperFunctions.nonlocalFilter2D(imTemp, sx, sy);
    case 7  %   dog
        %   using abandoned filter type: DoG
        Sxy = [19, 25];
        imSmooth = HelperFunctions.dogFilter2D(imTemp, Sxy);
    otherwise
        warning ('no nucleus smooth type mentioned, using non local')
        imSmooth = HelperFunctions.wienerFilter2D(imTemp, sx, sy);
end
% imDiff = im - imSmoothClean;
% figure; imshowpair(im, imSmoothClean, 'montage'); title (strcat('smooth as ', {' '}, 'ass'));
fim = imSmooth;
end