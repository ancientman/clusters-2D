function [localPeaks] = localPeakDetect(I, labelMat, globalIteration, outputFolder)
% imBgSub = medfilt2(I);
imBgSub = I;
mask = ones(size(labelMat));
mask(~labelMat) = 0;
imBgSub(~mask) = 0;
% imBgSub(~mask) = NaN; %gate the input image by the nuclear mask 
%Try difference of top hat filter for internal spot capturing
se = strel('disk', 8);
imPeaksTemp1 = imtophat(imtophat(imBgSub, se), se);
% imPeaksTemp1 = imPeaksTemp1./max(imPeaksTemp1);
se = strel('disk', 20);
imPeaksTemp2 = imtophat(imtophat(imBgSub, se), se);
% imPeaksTemp2 = imPeaksTemp2./max(imPeaksTemp2);
%Try difference of gaussian or difference of "Laplacian of Gaussian"%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imFilt1 = fspecial('log', [5, 5], 0.3); %log
% % imfilt1 = fspecial('gaussian', [5, 5], 0.3); %dog
% imPeaksTemp1 = imfilter(imBgSub, imFilt1, 'replicate');
% % imPeaksTemp1 = imPeaksTemp1./max(imPeaksTemp1);% 
% imFilt2 = fspecial('log', [5, 5], 0.5); %log
% % imfilt2 = fspecial('gaussian', [5, 5], 0.5); %dog
% imPeaksTemp2 = imfilter(imBgSub, imFilt2, 'replicate');
% % imPeaksTemp2 = imPeaksTemp2./max(imPeaksTemp2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% imPeaksTemp2 = imtophat(imtophat(imBgSub, se), se);
% imPeaksTemp2 = imPeaksTemp2./max(imPeaksTemp2);
% imPeaks2 = imPeaksTemp2;%imPeaksTemp1 .* imPeaksTemp2;
% imPeaks2(imPeaksTemp1<imPeaksTemp2) = 0;

% imPeaksTemp1(imPeaksTemp1== 0) = NaN;
% imPeaksTemp2(imPeaksTemp2 == 0) = NaN;
% imPeaks2 = imsubtract(imPeaksTemp2, imPeaksTemp1);
% imPeaks2(imPeaks2 == 0) = NaN;

% imUse = imgaussfilt(imPeaksTemp2, 0.5);
% mean0 = mean2(imUse);
% std0 = std2(imUse);
% imUse(imUse < (mean0 + 4*std0)) = 0;
% imPeaks5 = imfill(imUse, 'holes');
% imPeaks6 = bwareaopen(imPeaks5, 20);
imPeaks3 = imgradient(imPeaksTemp2, 'sobel');
imPeaks4 = imgaussfilt(imPeaks3, 0.3);
% imPeaks4 = imtophatfilt(imPeaks3, 0.7);
se = strel('disk', 4);
% imPeaks4 = imtophat(imPeaks3, se);
% imPeaks4(imPeaks4 == 0) = NaN;
mean1 = mean2(imPeaks4(imPeaks4 ~=0));
std1 = std2(imPeaks4);
howManySigma = 3; %determine the cutoff of peaks
peakThreshold = (mean1 + howManySigma*std1);
bwPeaks4 = imbinarize(imPeaks4, peakThreshold);
% imPeaks4(imPeaks4 < (mean1 + 2.5*std1)) = 0;
% imPeaks4(isnan(imPeaks4)) = 0;
bwPeaks5 = imfill(bwPeaks4, 'holes');
bwPeaks6 = imclose(bwPeaks5, se);
bwPeaks7 = bwareaopen(bwPeaks6, 30);
spotMask = bwPeaks7;
Bspot = bwboundaries(spotMask);
BspotSmooth = Bspot;
windowWidth = 5;
polynomialOrder = 2;
for ii = 1:length(Bspot)
    boundary = Bspot{ii};
    xx = boundary(:, 1);
    yy = boundary(:, 2);
    smoothX = sgolayfilt(xx, polynomialOrder, windowWidth);
    smoothY = sgolayfilt(yy, polynomialOrder, windowWidth);
    BspotSmooth{ii} = horzcat(smoothX, smoothY);
end
% figure;
% imshowpair(imBgSub, imPeaks7, 'montage');
% imhist(imPeaks2);
% bar(binLocation, count1);%, 'filled', 'LineWidth', 2);
%random experiments%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% features = extractLBPFeatures(imBgSub);
% peakThresh1 = imextendedmax(imPeaks,700, 4);
% peakThresh2 = imregionalmax(imPeaks, 4);
% [count, binLocation] = imhist(imBgSub, 256);
% cumDist = cdf('Normal', binLocation, count);
% figure; 
% yyaxis left;
% stem(binLocation, count, 'filled', 'LineWidth', 2); %bar(count); %
% ylim([0, 15000]);
% set(gca, 'YScale', 'log');
% hold on;
% yyaxis right;
% % plot(binLocation, cumDist,'r-', 'LineWidth', 2);
% cdfplot(count);
% xlim([0, 30000]);%binLocation(end)]);
% ylim([0.95, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; stem(binLocation, count, 'filled');
% xlim([0, 4000]);
% ylim();
% figure; 
% hold on;
% plot(features);
%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bwFuse1 = imfuse(imBgSub, imPeaks2, 'blend');
% B = bwboundaries(Mask);
% figure;
% imshow(imBgSub, []);
% hold on;
% for ii = 1:length(B)
%     boundary = B{ii};
%     plot(boundary(:,2), boundary(:,1), 'g', 'Linewidth', 2);
% end
% hold on;
% for ii = 1:length(BspotSmooth)
%     boundary = BspotSmooth{ii};
%     plot(boundary(:,2), boundary(:,1), 'r', 'Linewidth', 1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write stacks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movieFrame = imoverlay(uint16(im), (max(uint16(im)).*uint16(bwperim(spotMask))), 'r');
% movieFrame = imoverlay(movieFrame, (max(uint16(im)).*uint16(bwperim(Mask))), 'g');
%%%%%%%use this%%%%%%
modIm4Display = repmat(I, [1, 1, 1]);
movieFrame1 = imoverlay(uint16(modIm4Display(:,:,1)), 256.*uint16(bwperim(spotMask)), [1, 0, 0]);
movieFrame = imoverlay((movieFrame1), 256.*uint16(bwperim(mask)), [0, 1, 0]);
% outFileName = strcat(outputFolder , filesep, 'out', num2str(globalIteration), '.tif');
outFileName = strcat(outputFolder , filesep, 'Bcd2Xpeaks','.tif');
imwrite(movieFrame, outFileName, 'WriteMode','append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
localPeaks = spotMask;
mean = mean1;
std = std1;
end