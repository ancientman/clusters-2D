function imageCorrAnalysis()
maxIter = 2; % total iterations for recursive thresholding
% filePath = 'E:\Dropbox (Princeton)\Data\Rawdata_RM\beads\TS bead_4um_2.tif';
% filePath = 'E:\Dropbox (Princeton)\Data\Rawdata_RM\beads\TS bead_0.1um_3_zStack-1.tif';
filePath = 'E:\Dropbox (Princeton)\Data_from_Rahul\bcd_integration_time\f3_bcd2x_gfp_532ms.tif';
% filePath = 'E:\Dropbox (Princeton)\Data_from_Rahul\gfp_nls\f2_nls_gfp_38ms.tif';
stack = HelperFunctions.loadTiff(filePath); 
stack1 = rescale(stack(:,:,1));
%------------------------------------------------------------
%   Use fixed ROI
xAxis = 60:120;
yAxis = 60:120;
stack1Mod = stack1(xAxis, yAxis); 
stackMod = stack(xAxis, yAxis, :); 
[X, Y] = meshgrid(xAxis, yAxis);
roi(X,Y) = 1;
stack1Roi = zeros(size(stack1));
stack1Roi(X,Y) = 1; 
%------------------------------------------------------------
% Draw manual ROI
% imagesc(stack(:,:,1));
% axis('square');
% axis('off');
% roi = roipoly;
%------------------------------------------------------------
figure('Color', 'w')
imagesc(stack1+ max(stack1, [], 'all').*(bwperim(stack1Roi)));
axis('square');
hold off;
%-----------------------------------------------------------
stackMod = stack.*uint16(stack1Roi);
% stackMod = stack.*uint16(roi);
for i=1:size(stackMod, 3)
    stack(:,:,i) = medfilt2(stackMod(:,:,i));
%     fprintf('i = %d\n', i);
end
%   Initialize for multi threshold operation
%___________________________________________________
imThres = cell(1, maxIter);
for j = 1:maxIter
    imThres{j} = zeros(size(stack));
end
thres = zeros(maxIter, size(stack, 3));
thres(1,:) = 0;
%-----------------------------------------------------------

%___________________________________________________
%   Use single threshold
% howManySigmas = 1; % Change here
% imThres = cell(1, 1);
% for j = 1:1
%     imThres{j} = zeros(size(stack));
% end
% thres = zeros(1, size(stack, 3));
% for i =1:size(stack, 3)
%     imTemp = rescale(stack(:,:,i));
%     imTemp = imTemp - graythresh(imTemp);
%     imTemp(imTemp<0) = 0;
%     thres(i) = howManySigmas*std2(nonzeros(imTemp));
%     imTemp = imTemp - ...
%         howManySigmas*std2(nonzeros(imTemp));    
%     imTemp(imTemp<0) = 0;
%     imThres{1}(:,:,i) = imTemp;
% end
%___________________________________________________


%___________________________________________________
%   Calculate  multi-thresholds. technique#1
% if maxIter>1 && maxIter<=21
%     maxMultiLevel = maxIter;
%     for i= 1:size(stack, 3)
%         imTemp = rescale(stack(:,:,i));
%         thres(2:maxMultiLevel, i) = ...
%         multithresh(imTemp, maxMultiLevel-1);
%     end    
% elseif maxIter>21
%     maxMultiLevel = 21;
%     for i = 1:size(stack, 3)
%         imTemp = rescale(stack(:,:,i));
%         thresTemp(2:maxMultiLevel) = ...
%             multithresh(imTemp, maxMultiLevel-1);
%         thres(2:maxIter, i) = interp1(linspace(1, 2, maxMultiLevel-1), ...
%             thresTemp(2:maxMultiLevel), linspace(1, 2, maxIter-1));
%     end
% end
%-----------------------------------------------------------
%   Apply multi threshold to stacks
% for i = 1:size(stack, 3)    
%     imTemp = rescale(stack(:,:,i));
%     for j = 1:maxIter
%         imTemp = imTemp - thres(j, i);
%         imTemp(imTemp<0) = 0;
%         imThres{j}(:,:,i) = imTemp;
%     end
% end
%___________________________________________________

%___________________________________________________
%   Calculate  multi-thresholds. technique#2
for i = 1:size(stack, 3)    
    imTemp = rescale(stack(:,:,i));
    for j = 1:maxIter
        if j == 1
            thres(j, i) = 0;
        elseif j == 2       
            thres(j, i) = graythresh(imTemp); %  image threshold (otsu) 
        else
            thres(j, i) = thres(2, i)+std2(nonzeros(imTemp))*0.2*j;
        end
        imTemp = imTemp - thres(j, i);
        imTemp(imTemp<0) = 0;
        imThres{j}(:,:,i) = imTemp;
        
    end
end
%___________________________________________________


%___________________________________________________
%   Temporal autocorrelation
aCorr = cell(1, length(imThres));
G0 = cell(1, length(imThres));
for j = 1:length(imThres)
    aCorr{j} = zeros(size(imThres{j}));
    G0{j} = zeros(size(imThres{j}, 3), 1);
end
for j = 1:length(imThres)
    for i=1:size(imThres{j}, 3)
        im = (imThres{j}(:,:, i));
        imSub=im-mean2(nonzeros(im)); %    subtract mean
        imSub(imSub<0) = 0;
        imNorm=imSub./sqrt(sum(nonzeros(imSub(:)).^2)); %  normalize magnitude
        fft_I=fft2(imNorm); %   compute fft2
        if i==1
            imFFT1 = fft_I;
        end    
        aCorr{j}(:,:,i)=real(fftshift(ifft2(imFFT1.*conj(fft_I)))); %  compute autocorrelation
        G0{j}(i) = max(aCorr{j}(:,:,i), [], 'all');
    end
end
%-----------------------------------------------------------
%   Axis transformation

for j = 1:length(imThres)
    figure('Color', 'w'); 
    plot((1:length(G0{j})).*0.532, G0{j});
    xlabel('{\Delta}{\tau}');
    ylabel('g_{(0,0)}{\tau}');
    title('G(0) bgSub thresh')
end
%___________________________________________________

%___________________________________________________

%___________________________________________________
plotBreak = min(maxIter, 5); % use 5 (5 rows per subplot figure)
sRows = maxIter;
sCols = 4;
g0 = zeros(maxIter, 1);
sigma = zeros(maxIter, 1);
fwhm = zeros(maxIter, 1);

i = 1;
while i <= maxIter     
    im(im<thres(i)) = 0;
    imSub=im-mean2(nonzeros(im)); %    subtract mean
    imSub(imSub<0) = 0;
    imNorm=imSub/sqrt(sum(nonzeros(imSub(:)).^2)); %  normalize magnitude
    fft_I=fft2(imNorm); %   compute fft2
    sCorr=real(fftshift(ifft2(fft_I.*conj(fft_I)))); %  compute autocorrelation
    if thres(i)<max(im,[], 'all')
        maxRows = plotBreak;
        sRows = maxRows;
        j = floor(i/maxRows - 0.01);
        plotIter = i-j*plotBreak;
        [g0(i), sigma(i), fwhm(i)] = subplotFun1(im, stackMod(:,:,1), sCorr, thres(i), sRows, sCols, plotIter, maxRows);
    else
        break;
    end
        i = i+1;
end
thres(g0==0) = 0;
% plot1(g0, sigma, fwhm, thres);

aCorr= zeros(size(stackMod));

for i=1:size(stack, 3)
    im = rescale(stackMod(:,:, i));
    im = im - 1.5*graythresh(im);
    im(im<0) = 0;
    imSub=im-mean2(nonzeros(im)); %    subtract mean
    imSub(imSub<0) = 0;
    imNorm=imSub./sqrt(sum(nonzeros(imSub(:)).^2)); %  normalize magnitude
    fft_I=fft2(imNorm); %   compute fft2
    if i==1
        imFFT1 = fft_I;
    end    
    aCorr(:,:,i)=real(fftshift(ifft2(imFFT1.*conj(fft_I)))); %  compute autocorrelation
    G0t(i) = max(aCorr(:,:,i), [], 'all');
end
figure('Color', 'w'); 
plot(G0t);
xlabel('{\Delta}{\tau}');
ylabel('g_{(0,0)}{\tau}');
title('G(0) fixed threshold')
% figure; 
% subplotFun2(aCorr);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stack1deNo = medfilt2(stack1, [23, 23]);
% corr1 = xcorr2(stack1);
% corrDeNo = xcorr2(stack1deNo);

% X = 0.042*(1:size(stack, 2));
% Y = 0.042*(1:size(stack, 1));
% sRows = 2;
% sCols = 3;
% figure('Color','w');
% subplot(sRows, sCols, 1)
% imagesc(X, Y, stack1);
% title('Original image');
% xlabel('{\mu}m');
% ylabel('{\mu}m');
% axis('square');
% subplot(sRows, sCols, 4)
% imagesc(X, Y, stack1deNo);
% title('Median filtering');
% xlabel('{\mu}m');
% ylabel('{\mu}m');
% axis('square');
% 
% subplot(sRows, sCols, 2)
% [p1, e1] = histcounts(stack1, 'Normalization','probability');
% % figure('Color', 'w');
% binC1 = e1 + (e1(2) - e1(1))/2;
% % plot(binC1(1:end-1), p1, 'r-');
% hold on;
% H1=area(binC1(1:end-1), p1);
% set(H1(1),'FaceColor',[0.9 0.9 0.7]);
% hold on;
% t1 = graythresh(stack1);
% idx=binC1<t1;
% Ht1=area(binC1(idx),p1(idx));
% set(Ht1(1),'FaceColor',[0.7 0.7 0.5]);
% xlabel('Normalized counts')
% ylabel('Probability density')
% title('thresholding')
% axis('square')
% xlim([0 1]);
% % legend('nucleus', 'background');
% hold on;
% subplot(sRows, sCols, 3)
% bin1 = imbinarize(stack1, t1);
% imshow(bin1);
% title('Nucleus segmentation')
% axis('square');
% subplot(sRows, sCols, 5)
% [p2, e2] = histcounts(stack1deNo, 'Normalization','probability');
% binC2 = e2 + (e2(2) - e2(1))/2;
% % plot(binC2(1:end-1), p2, 'k-');
% hold on;
% H2=area(binC2(1:end-1), p2);
% set(H2(1),'FaceColor',[0.9 0.9 0.7]);
% hold on;
% t2 = graythresh(stack1deNo);
% idx=binC2<t2;
% Ht2=area(binC2(idx),p2(idx));
% set(Ht2(1),'FaceColor',[0.7 0.7 0.5]);
% xlabel('Normalized counts')
% ylabel('Probability density')
% title('thresholding')
% axis('square')
% xlim([0 1]);
% subplot(sRows, sCols, 6)
% bin1deNo = imbinarize(stack1deNo, t2);
% imshow(bin1deNo);
% title('Nucleus segmentation')
% axis('square');
% hold off;
% 
% figure('Color', 'w');
% % frameMin = (min(stack, [], [1 2]));
% % frameMax = (max(stack, [], [1 2]));
% % reStackAll = rescale(stack, 'InputMin', frameMin,'InputMax', frameMax);
% % meanAll = mean(reStackAll,3);
% 
% imAdd = zeros(size(stack1));
% for i = 1:size(stack, 3)
%     reStack = rescale(stack(:,:,i));
%     imAdd = imAdd + rescale(stack(:,:,i));    
% end
% meanAll = imAdd/size(stack, 3);
% imagesc(X, Y, meanAll);
% title('Mean projection');
% xlabel('{\mu}m');
% ylabel('{\mu}m');
% axis('square')
% hold off;

% stack1=double(stack1); %convert to double
% t1 = graythresh(stack1); %  image threshold (otsu)
% 
% 
% im = stack1-t1;
% im=im-mean(im(:)); %subtract mean
% im=im/sqrt(sum(im(:).^2)); %normalize magnitude
% fft_I=fft2(im); %compute fft2
% aCorr=real(fftshift(ifft2(fft_I.*conj(fft_I)))); %compute autocorrelation

% if strcmp(bypass, 'yes')
%     maxIter = 1;
% end
% plotBreak = min(maxIter, 5); % use 5 (5 rows per subplot figure)
% sRows = maxIter;
% sCols = 4;
% i = 1;
% g0 = zeros(maxIter, 1);
% sigma = zeros(maxIter, 1);
% fwhm = zeros(maxIter, 1);
% thres = zeros(maxIter, 1);
% thres(1) = 0;
% if maxIter>1 && maxIter<=21
%     maxMultiLevel = maxIter;
%     thres(2:maxMultiLevel) = multithresh(stack1Mod, maxMultiLevel-1);
% elseif maxIter>21
%     maxMultiLevel = 21;
%     thresTemp(2:maxMultiLevel) = multithresh(stack1Mod, maxMultiLevel-1);
%     thres(2:maxIter) = interp1(linspace(1, 2, maxMultiLevel-1), thresTemp(2:maxMultiLevel), linspace(1, 2, maxIter-1));
% end
% 
% while i <= maxIter   
%     im = stack1Mod;     
% %      if i == 1
% %          im = stackMod;
% %          thres(i) = 0;
% %      elseif i == 2        
% %         thres(i) = graythresh(im); %  image threshold (otsu)        
% %      else
% %         thres(i) = thres(2)+std2(nonzeros(stackMod))*0.2*i;       
% %     end
% 	im(im<thres(i)) = 0;
%     imSub=im-mean2(nonzeros(im)); %    subtract mean
%     imSub(imSub<0) = 0;
%     imNorm=imSub/sqrt(sum(nonzeros(imSub(:)).^2)); %  normalize magnitude
%     fft_I=fft2(imNorm); %   compute fft2
%     sCorr=real(fftshift(ifft2(fft_I.*conj(fft_I)))); %  compute autocorrelation
%     if thres(i)<max(im,[], 'all')
%         maxRows = plotBreak;
%         sRows = maxRows;
%         j = floor(i/maxRows - 0.01);
%         plotIter = i-j*plotBreak;
%         [g0(i), sigma(i), fwhm(i)] = subplotFun1(im, stack1Mod, sCorr, thres(i), sRows, sCols, plotIter, maxRows);
%     else
%         break;
%     end
%     i = i+1;
% end
% thres(g0==0) = 0;
% % plot1(g0, sigma, fwhm, thres);
% 
% aCorr= zeros(size(stackMod));
% 
% for i=1:size(stack, 3)
%     im = rescale(stackMod(:,:, i));
%     im = im - 1.5*graythresh(im);
%     im(im<0) = 0;
%     imSub=im-mean2(nonzeros(im)); %    subtract mean
%     imSub(imSub<0) = 0;
%     imNorm=imSub./sqrt(sum(nonzeros(imSub(:)).^2)); %  normalize magnitude
%     fft_I=fft2(imNorm); %   compute fft2
%     if i==1
%         imFFT1 = fft_I;
%     end    
%     aCorr(:,:,i)=real(fftshift(ifft2(imFFT1.*conj(fft_I)))); %  compute autocorrelation
%     G0(i) = max(aCorr(:,:,i), [], 'all');
% end
% figure('Color', 'w'); 
% plot(G0);
% xlabel('{\Delta}{\tau}');
% ylabel('g_{(0,0)}{\tau}');
% title('100 nm bead')
% % figure; 
% % subplotFun2(aCorr);
end

function subplotFun2(aCorr)
for i = 1:size(aCorr, 3)
    subplot(size(aCorr, 3), 1, i);
    surf(aCorr(:,:,i), 'EdgeColor', 'none');
    ylabel(['{\Delta}{\tau} = ', num2str(i)]);
end
end

function [g0, sigma, fwhm] = subplotFun1(im, stackMod, sCorr, t1, sRows, sCols, i, maxIter)
if i==1
    figure('Color','w');
end
xScale = 0.042;
yScale = 0.042;

X = yScale*(1:size(stackMod, 2));
Y = xScale*(1:size(stackMod, 1));
s1 = subplot(sRows, sCols, (i-1)*sCols+1);
[p1, e1] = histcounts(stackMod, 'Normalization','probability');
binC1 = e1 + (e1(2) - e1(1))/2;
% plot(binC1(1:end-1), p1, 'r-');
hold on;
H1=area(binC1(1:end-1), p1);
set(H1(1),'FaceColor',[0.9 0.9 0.7]);
hold on;
if t1~=0
    idx=binC1<t1;
    Ht1=area(binC1(idx),p1(idx));
    set(Ht1(1),'FaceColor',[0.7 0.7 0.5]);
end
if i==maxIter
    xlabel('Counts')
end
ylabel('{\rho}')
if i==1
    title('Threshold')
end
axis('square')
xlim([0 max(binC1)]);
hold on;
s2 = subplot(sRows, sCols, (i-1)*sCols+2);
imagesc(X, Y, im);
if i==1 && t1 == 0
        title('raw image')
elseif i==2 && t1 ~= 0
    title('Subtracted');
end
if i==maxIter
    xlabel('{\mu}m');
end
ylabel('{\mu}m');
axis('square');
s3 = subplot(sRows, sCols, (i-1)*sCols+3);
xLag = (-ceil((size(sCorr, 1)-1)/2): floor((size(sCorr, 1)-1)/2)).*xScale;
yLag = (-ceil((size(sCorr, 2)-1)/2): floor((size(sCorr, 2)-1)/2)).*yScale;
surf(xLag, yLag, sCorr, 'EdgeColor', 'none');
if i==1
    title('Autocorrelation');
end
if i==maxIter
    xlabel('{\xi} {\mu}m');
    ylabel('{\eta} {\mu}m');
end
axis('square');
s4 = subplot(sRows, sCols, (i-1)*sCols+4);
% p1 = plot((1:size(aCorr, 2)).*xScale, aCorr(floor((size(aCorr, 1)-1)/2), :), '-k');

p1 = plot(xLag, sCorr(ceil((size(sCorr, 1)-1)/2), :), '-k');
hold on;
% f1 = fit(xLag', sCorr(ceil((size(sCorr, 1)-1)/2), :)', 'gauss1');
%------------------------------------------------
%   Turn on if gauss1 can't fit profile
%------------------------------------------------
options = fitoptions('gauss2', 'Lower', [0 -Inf 0 0 -Inf 0]);
f1 = fit(xLag', sCorr(ceil((size(sCorr, 1)-1)/2), :)', 'gauss2', options);
%------------------------------------------------

g0 = f1.a1;
sigma = f1.c1/(sqrt(2));
fwhm = 2*sqrt(log(2)) * f1.c1;
plot(f1);
legend('hide')
if i==1
    title('Fit');
end
if i==maxIter
    xlabel('{\xi} {\mu}m');
end
ylabel('g');
axis('square');
hold off;
pos1 = get(s1, 'position');
annotation('textbox', pos1.*[1, 1, 1, 1], 'EdgeColor', 'none', 'FitBoxToText', 'on','verticalalignment', 'top', 'horizontalalignment', 'right', 'String', ['t = ', num2str(t1)])
pos4 = get(s4, 'position');
annotation('textbox', pos4.*[1, 1, 1, 1], 'EdgeColor', 'none', 'FitBoxToText', 'on','verticalalignment', 'top', 'horizontalalignment', 'right', 'String', ['fwhm = ', num2str(fwhm)])
end

function plot1(g0, sigma, fwhm, thres)
f1 = figure('Color', 'w');
colorLeft = [0.3, 0.8, 0.6];
colorRight = [0.4, 0.2, 0.9];
set(f1,'defaultAxesColorOrder',[colorLeft; colorRight]);
yyaxis left;
% p1 = plot(nonzeros(thres), nonzeros(g0), '--o');
p1 = plot((thres), (g0), '--o');
p1.LineWidth = 2;
p1.Color = [colorLeft];
ylabel('g(0,0)');
yyaxis right;
% p2 = plot(nonzeros(thres), nonzeros(fwhm), '-.^');
p2 = plot((thres), (fwhm), '-.^');
ylabel('fwhm ({\mu}m)');
% p2 = plot((thres), (sigma), '-.^');
% ylabel('{\sigma} ({\mu}m)');
p2.LineWidth = 2;
p2.Color = [colorRight];
xlabel('Threshold (a. u.)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
axis('square');
end