function [fim] = seedNuclearMask2(im, el, nucSize, blurControl, t, stopMark)
imScale = rescale(im);
sx = el; sy = el;
se = strel('disk',el);
imBlur = HelperFunctions.nonlocalFilter2D(imScale, sx, sy);%nonlocalFilter2D(im, el, el);
bw1 = imbinarize(imBlur, 'adaptive', 'Sensitivity', blurControl);
bw1 = imfill(bw1, 'holes');
bw2 = imclose(bw1, se);
bw2 = imfill(bw2, 'holes');
windowSize = 15; %hard coded (roughly 1/5 the nucleus diameter)
expFactor = 1.95;
kernel = ones(windowSize) / (windowSize^expFactor);
imBlur2 = conv2(single(bw2), kernel, 'same');
bw3 = imBlur2 > blurControl;
bw3 = imfill(bw3, 'holes');
bw4 = imfill(bw3, 'holes');
fim = bw4;
if t == stopMark
    fprintf("in seedNuclearMask at t = %d\n", t);
end
end