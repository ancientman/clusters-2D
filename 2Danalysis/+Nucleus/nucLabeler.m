function [CC] = nucLabeler(mask, I)
% mask2 = imfill(imdilate(mask,strel('disk',2)), 'holes');
mask2 = mask;
CC = bwconncomp(mask2);
% meanNucIntensity = zeros(CC.NumObjects, 1);
% for k = 1:CC.NumObjects
%     meanNucIntensity(k) = mean(rescale(I(CC.PixelIdxList{k})));
% end
% L = labelmatrix(CC);
% S = regionprops(CC, 'centroid');
% centroids = cat(1, S.Centroid);
% % Lrgb = label2rgb(L,'jet','w','shuffle');
% iOverlay = labeloverlay(rescale(I), L);
% figure; imshow(iOverlay); title ('labels');
% hold on;
% plot (centroids (:, 1), centroids (:, 2), 'b*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L = bwlabel(nucMask);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Fast Marching%%%%%%%%%%%%%%%%%%%%%%%%
% W = graydiffweight(I,mask);
% BW = imsegfmm(W,mask,0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Redundant Module Dont Use%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imThresh1 = adaptthresh(I, 1);
% I = imgaussfilt(imThresh1);
% gmag = imgradient(I);
% D = bwdist(mask);
% DL = watershed(D);
% bgm = DL == 0;
% gmag2 = imimposemin(gmag, bgm | mask);
% L = watershed(gmag2);%(gmag2);
% labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*mask;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lrgb = label2rgb(L,'jet','w','shuffle');
% figure; imshow(Lrgb);
% I4 = labeloverlay(I,labels);
% figure; imshowpair(rescale(I), I4, 'montage'); 
% title('watershed segmented');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nucLabel = L;
% nucMean = meanNucIntensity;
% nucCentroid = centroids;

end