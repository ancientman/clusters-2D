function [fim, t1] = noisyNuclearMask(im, mask, blurControl, nucSize, clearBorder, t, stopMark)
% mask2 = activecontour(im,mask,200, 'Chan-Vese', 'SmoothFactor',5, 'ContractionBias',-1);
mask2 = activecontour(im, mask, 100, 'edge', 'SmoothFactor', 1, 'ContractionBias', 1);
imScale = rescale(im);
imFilt = imgaussfilt(imadjust(imScale));
t1 = adaptthresh(imFilt, blurControl, 'ForegroundPolarity', 'bright', 'Statistic', 'gaussian'); 
bw1 = imbinarize(imFilt, t1);
bw2 = mask & mask2;%mask & bw1;%
bw3 = imfill(bw2, 'holes');
bw3 = bwareaopen(bw3, ceil(nucSize/3));
D = -bwdist(~bw3);
D(bw3)= -inf;
DL = watershed(D);
bgm = DL == 0;
gmag = imimposemin(imgradient(imScale), mask2 | bgm);
DL2 = watershed(gmag);
S = regionprops(DL2,'Area');
Area = vertcat(S.Area);
[~,bgValue] = max(Area);
bw4 = any((DL2~=bgValue & DL2~=0),3);
bw5 = imfill(bw4, 'holes');
bw5 = bwareaopen(bw5, nucSize);
if clearBorder
    bw5 = imclearborder(bw5);
end
fim = bw5;
if t == stopMark
    fprintf("in noisyNuclearMask");
end
end