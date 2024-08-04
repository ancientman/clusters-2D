function [fim] = bilateralFilter2D(im,sx,sy)
BW1 = imgaussfilt(im, sx);
se = strel('rectangle',[sx, sy]);
BW2 = imdilate (BW1, se);
BW3 = imclose(BW2, se);
BW4 = imbinarize(uint8(BW3));
BWbg = single(im).*(~BW4);
bgVar = std2(single(im).*(~BWbg))^2;
DegreeofSmooth = sx*sy*bgVar;
spatialSigma = 2; %controls smoothing
J = imbilatfilt(im,DegreeofSmooth, spatialSigma);
fim = J;
end