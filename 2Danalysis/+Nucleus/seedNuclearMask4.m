function [fim, bgm] = seedNuclearMask4(im, el, nucSize, blurControl, t, stopMark)
imScale = rescale(im);
imBlur = HelperFunctions.gaussianFilter2D((imScale), el, el);
imBlur = rescale(imBlur);
% imBlur = HelperFunctions.dogFilter2D((imScale), [7, 9]);
% imBlur = imadjust(imBlur);
% bw0 = imbinarize(adapthisteq(imBlur));
% 
% bw1 = imextendedmax(imBlur, blurControl*max(imBlur,[], 'all')/5);
% bw1 = imfill(bw1, 'holes');
% bw2 = false(size(im));
% bw2(bw0|bw1) = 1;

% gaussian1 = fspecial('Gaussian', 101, 91);
% gaussian2 = fspecial('Gaussian', 101, 81);
% dog = gaussian1 - gaussian2;
% dogFilterImage = conv2(double(imBlur), dog, 'same');


bw1 = imbinarize(rescale(imBlur));

%-------------------------------------------------------
% imAd = imadjust(imBlur, [0, graythresh(imBlur)]);
imAd = imBlur.*(imadjust(imBlur, [0, graythresh(imBlur)]).^10);
imAd = HelperFunctions.gaussianFilter2D(rescale(imAd), el, el);
% bwa = imextendedmax(imAd, graythresh(imAd));
bwa = imbinarize(rescale(imAd), 1*graythresh(imBlur));
bwa = imfill(bwa, 'holes');
bwa = imopen(bwa,strel([el el]));
bwa = imclose(bwa,strel([el el]));
bwa = imfill(bwa, 'holes');
Da = -bwdist(~bwa);
Da(bwa) = -inf;
La = watershed(Da);
bgma = La == 0;
bgma = imdilate(bgma,strel([el el]));
bwa(imbinarize(imAd) & bgma) = 0;
bwa = imfill(bwa, 'holes');
fim = bwa;
bgm = bgma;
%-------------------------------------------------------

bw2 = bwa;
bw0 = bw2;

bw3 = imfill(bw2, 'holes');
bw4 = imclose(bw3,strel([el el]));
bw4 = imerode(bw4, strel([el el]));
bw4 = imfill(bw4, 'holes');
D0 = -bwdist(~bw4);
D0(bw4) = -inf;
DL0 = watershed(D0);
bgm0 = DL0==0;
D = -bwdist(~bw1);
D(bw1)= -inf;
DL = watershed(D);
bgm = DL == 0;
bgm = imdilate(bgm,strel([el el]));
gmag0 = imimposemin(imgradient(imScale), bw1 | (bgm0 | bgm));
DL2 = watershed(gmag0);
S = regionprops(DL2,'Area');
Area = vertcat(S.Area);
[~,bgValue] = max(Area);
bw5 = any((DL2~=bgValue & DL2~=0),3);
fim = bw5;
% if t == stopMark
%     fprintf("in seedNuclearMask4 at t = %d\n", stopMark);
% end
end