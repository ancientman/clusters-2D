function [nucProp] = nucProp(CC, I, t, stopMark)
meanNucIntensity = zeros(CC.NumObjects, 1);
meanAbsNucIntensity = zeros(CC.NumObjects, 1);
stDevNucIntensity = zeros(CC.NumObjects, 1);
maxNucIntensity = zeros(CC.NumObjects, 1);
for k = 1:CC.NumObjects
    meanNucIntensity(k) = mean(rescale(I(CC.PixelIdxList{k})));
    meanAbsNucIntensity(k) = mean(I(CC.PixelIdxList{k}));
    stDevNucIntensity(k) = std(rescale(I(CC.PixelIdxList{k})));
    maxNucIntensity(k) = max(I(CC.PixelIdxList{k}));
end
L = labelmatrix(CC);
S = regionprops(CC, 'centroid');
centroids = cat(1, S.Centroid);
A = regionprops(CC, 'area');
nucArea = cat(1, A.Area);
% Lrgb = label2rgb(L,'jet','w','shuffle');
% iOverlay = labeloverlay(rescale(I), L);
% figure; imshow(iOverlay); title ('labels');
% hold on;
% plot (centroids (:, 1), centroids (:, 2), 'b*');
nucProp.labelMat = L;
nucProp.meanIntensity = meanNucIntensity;
nucProp.meanAbsIntensity = meanAbsNucIntensity;
nucProp.stdIntensity = stDevNucIntensity;
nucProp.maxIntensity = maxNucIntensity;
nucProp.area = nucArea;
nucProp.centroid = centroids;
if t == stopMark
    fprintf("in nucProp");
end
end