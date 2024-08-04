function [fim] = seedNuclearMask(im, el, nucSize, blurControl, t, stopMark)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a crude segmentation technique that produces smoothass boundaries to
% segmented nuclei. 
% This corresponds to option "watershed" of nucleusDetectMethod
% This function segments the nuclei specifically by watershed. 
% This can be used for most applications.
% It first blurs the input image then performs watershed operation. 
% Then uses a convolution of the blurred watershed labels to create the
% final segmented image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imScale = rescale(im);
sx = 5*el; sy = 5*el;
imBlur = HelperFunctions.wienerFilter2D(imScale, sx, sy);
se = strel('disk', el);
bw1 = imerode(imBlur, se);
bw2 = imreconstruct(bw1, imBlur);
bw3 = imdilate(bw2, se);
bw4 = imreconstruct(imcomplement(bw3),imcomplement(bw2));
bw4 = imcomplement(bw4);
bw5 = imbinarize(imgaussfilt(bw4, 2)); 
se2 = strel(ones(el,el));
bw6 = imclose(bw5,se2);
bw6 = imerode(bw6,se2);
windowSize = 25; %hardcoded
kernel = ones(windowSize) / windowSize ^1.8; % original = 1.8
imBlur2 = conv2(single(bw6), kernel, 'same');
bw7 = rescale(imBlur2) > blurControl;
% bw7 = bwareaopen(bw7, nucSize);
fim = bw7;
if t == stopMark
    fprintf("in seedNuclearMask at t = %d\n", t);
end
end