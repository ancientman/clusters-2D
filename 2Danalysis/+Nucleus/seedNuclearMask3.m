function [fim] = seedNuclearMask3(im, el, nucSize, blurControl, t, stopMark)
imScale = rescale(im);
imBlur = HelperFunctions.gaussianFilter2D(imScale, el, el);
imThresh = adaptthresh(imScale, blurControl);
bw1 = imbinarize(imBlur, imThresh);
bw1 = imfill(bw1, 'holes');
SE = strel('disk', el);
bw2 = imerode(bw1, SE);
bw2 = imerode(bw2, SE);
% bw2 = imerode(bw2, SE); %option for super erode
% bw3 = bwconvhull(bw2, 'objects');
fim = bw2;
if t == stopMark
    fprintf("in seedNuclearMask at t = %d\n", stopMark);
end
end