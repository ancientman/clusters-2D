function [fim, totalSpots] = adapThresVis(im, nucLabel)
maxIter = 4; %6 use for bcd % nuclear threshold
technique = 2;
imF= im;
imF(nucLabel == 0) = 0;
imF= medfilt2(imF);

thres = cell(1, max(nucLabel, [], 'all'));
imThresNuc = cell(1, max(nucLabel, [], 'all'));
% imThresLoc = cell(1, max(nucLabel, [], 'all'));
for i = 1:max(nucLabel, [], 'all') 
    thres{i} = zeros(maxIter, 1);
    imThresNuc{i} = zeros([size(im), maxIter]);
%     imThresLoc{i} = zeros([size(im), maxIter]);
end

if technique == 1
%___________________________________________________
%   Calculate  multi-thresholds. technique#1: auto-otsu
imThresNuc = nucThres(nucLabel, im, maxIter, thres, imThresNuc);       
%___________________________________________________
elseif technique == 2
%___________________________________________________
%   Calculate  multi-thresholds. technique#2: deviation from mean
imThresNuc = nucThres2(nucLabel, im, maxIter, thres, imThresNuc);       
end
 %___________________________________________________
%   Calculate  local thresholds
[imThresLoc, totalSpots] = locThres(nucLabel, imThresNuc); 
%___________________________________________________
fim = imThresLoc;
end

%___________________________________________________
function imThres = nucThres(nucLabel, im, maxIter, thres, imThres)
for i= 1:max(nucLabel, [], 'all')
    imTemp = im;
    imTemp(nucLabel~=i) = 0;
    imTemp = rescale(imTemp);
    if maxIter>1 && maxIter<=21
        maxMultiLevel = maxIter;
        thres{i}(2:maxIter) = ...
        multithresh(imTemp, maxMultiLevel-1);
    elseif maxIter>21
        maxMultiLevel = 21;
        thresTemp(2:maxMultiLevel) = ...
            multithresh(imTemp, maxMultiLevel-1);
        thres{i}(2:maxIter) = interp1(linspace(1, 2, maxMultiLevel-1), ...
            thresTemp(2:maxMultiLevel), linspace(1, 2, maxIter-1));
    end
%     figure('Color', 'w');
    [row, col] = find(nucLabel==i);
    nucCrop = imcrop(imTemp, [min(col), min(row), ...
                (max(col) - min(col)), (max(row) - min(row))]);
            
    %   Plot histogram
%     plotImgHist(nucCrop, thres{i});
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     figure('color', 'w');
    for j=1:maxIter
        imTemp = imTemp - thres{i}(j);
        imTemp(imTemp<0) = 0;
        imTemp = rescale(imTemp);
        imThres{i}(:,:,j) = imTemp;        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %   Visualization
%         if j<=9
%             thresCrop = imcrop(imTemp, [min(col), min(row), ...
%                     (max(col) - min(col)), (max(row) - min(row))]);
%             subplot(3, 3, j)
%             imshow(thresCrop,[]);
%         end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end    
%     ax = gca;
%     ax = gca;
%     ax.FontSize = 12;
%     ax.LineWidth = 1.5;
%     box(ax,'on');
%     grid off;
%     x0 = 100;
%     y0= 100;
%     plotWidth=500;
%     plotHeight=500;
%     set(gcf,'position',[x0,y0,plotWidth,plotHeight])
%     hold off;
end
end


function imThres = nucThres2(nucLabel, im, maxIter, thres, imThres)       
thresholdFraction = 0.9;
howManySigma = 1;
for i= 1:max(nucLabel, [], 'all')    
    imTemp = im;
    imTemp(nucLabel~=i) = 0;
    for j = 1:maxIter
        imTemp = rescale(imTemp);
        meanNuc = mean2(nonzeros(imTemp));
        stdNuc = std(nonzeros(imTemp), 1, 'all');    
%         thres{i}(j) = thresholdFraction*prctile(nonzeros(imTemp), 0.5,'all'); 
        thres{i}(j) = meanNuc+stdNuc*0.5*howManySigma;
        imTemp = imTemp - thres{i}(j);
        imTemp(imTemp<0) = 0;
        imTemp(nucLabel ~= i) = 0;
        imThres{i}(:,:,j) = rescale(imTemp);
    end
end   
end

function    [imThresLoc, totalSpots] = locThres(nucLabel, im)
imThresLoc = zeros(size(nucLabel));
filtWinSize = 15;
maxIter = 9; %%%%%% (4 for bcd)change here. use 2, 3 or 4 
maxIter = maxIter^2; % local threshold
nNucThres = size(im{1}, 3);
imTempStack = cell(1,maxIter);
spotBin = cell(1,(max(nucLabel, [], 'all')));
bgMean = cell(1, (max(nucLabel, [], 'all')));
totalSpots = zeros(1, max(nucLabel, [], 'all'));
for i= 1:max(nucLabel, [], 'all') % Nuclei
    cropSpotValStack = cell(1, maxIter);
    cropSpotBinStack = cell(1, maxIter);
    [row, col] = find(nucLabel==i);
    kk = 1;
%     cropSpotStack = zeros([max(row)-min(row)+1, max(col)-min(col)+1, maxIter]);
    for j = 1:3%1:nNucThres %nNucThres:nNucThres%     Nuclear Thresholds
        imTemp = im{i}(:,:,j);
        imTemp(nucLabel~=i) = 0;
        imTemp = rescale(imTemp);
%         figure('Color', 'w');
        for p= 1:maxIter  % Local threshold iterations
            imAvTemp = conv2(imTemp, ones(filtWinSize)/filtWinSize^2, 'same');
            imStdTemp = stdfilt(imTemp, ones(filtWinSize));
            %~~~~~~~~~~~~~~~~~~~~~~~
%             subplot(sqrt(maxIter), sqrt(maxIter), p)
%             noisePlotter(imAvTemp, imStdTemp);
            %~~~~~~~~~~~~~~~~~~~~~~~
            thres = imAvTemp + 0.25.*imStdTemp;
            imTemp = imTemp - thres;
            imTemp(imTemp<0) = 0;
            imTemp = rescale(imTemp);
            imTempStack{j}(:,:,p) = imTemp;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Cut one nucleus out            
            imFiltCrop = imcrop(imTempStack{j}(:,:,p), [min(col), min(row), ...
                (max(col) - min(col)), (max(row) - min(row))]);
            cropSpotValStack{j}(:,:,p) = imFiltCrop;

%             cropSpotBinStack{j}(:,:,p) = imbinarize(imFiltCrop);
            imFiltCropBin = imFiltCrop;
            imFiltCropBin(imFiltCrop>0) = 1;
            cropSpotBinStack{j}(:,:,p) = imFiltCropBin;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Calculate background pixels
            if p==1 && j==1
%                 spotBin{i}(:,:,p) = imbinarize(imTemp);
                imTempBin = imTemp;
                imTempBin(imTemp>0) = 1;
                spotBin{i}(:,:,p) = imTempBin;
                imTemp2 = im{i}(:,:,j);
                bgMean{i} = mean2(imTemp2(nucLabel==i & spotBin{i}(:,:,p)==0));
            end
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if p==1 && j==1
                cropSpotVal1 = cropSpotValStack{j}(:,:,p);
                cropSpotBin1 = cropSpotBinStack{j}(:,:,p);
            end

            CC = bwconncomp(cropSpotBinStack{j}(:,:,p), 8);
            numSpots(kk) = CC.NumObjects; 
       
            [ssimVal(kk),~] = ssim(cropSpotVal1,  cropSpotValStack{j}(:,:,p));
            corrVal(kk) = corr2(cropSpotVal1,  cropSpotValStack{j}(:,:,p));
            [ssimBin(kk),~] = ssim(single(cropSpotBin1),  single(cropSpotBinStack{j}(:,:,p)));
            corrBin(kk) = corr2(single(cropSpotBin1),  single(cropSpotBinStack{j}(:,:,p)));
                      
            A = regionprops(CC, 'area');
            spotAreaAll = cat(1, A.Area);
            spotAreaAv(kk) = mean(spotAreaAll);
            spotAreaStd(kk) = std(spotAreaAll);
            spotAreaSem(kk) = spotAreaStd(kk)./length(spotAreaAll);

            nucBinCrop = imcrop(nucLabel, [min(col), min(row), ...
                (max(col) - min(col)), (max(row) - min(row))]);

            kk = kk+1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%             %   Visualization   
            subplot(sqrt(maxIter), sqrt(maxIter), p)
            imshow(rescale(cropSpotValStack{j}(:,:,p)) + bwperim(nucBinCrop));
%             histogram(cropSpotValStack{j}(:,:,p),0.05:0.05:1);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
    end
    
%------------------------------------------------------------
% %   Visualization   
    figure('color', 'w'); 
    plot (numSpots, 'k');
    ylabel ('# Spots per nucleus');
    yyaxis right;
    ee = errorbar(spotAreaAv, spotAreaStd);
    ee.CapSize = 0;
    ylabel('Mean spot area (pixels)');            
    xlabel('Threshold iteration #');
    
    ax = gca;
    ax = gca;
    ax.FontSize = 10;
    ax.LineWidth = 1;
    box(ax,'on');
    grid off;
    x0 = 100;
    y0= 100;
    plotWidth=250;
    plotHeight=250;
    set(gcf,'position',[x0,y0,plotWidth,plotHeight])
    hold off;    

%------------------------------------------------------------
% %   Visualization
    figure; plot(ssimVal); title('ssim val');
    figure; plot(ssimBin); title('ssim bin');
    figure; plot(corrVal); title('corr val');
    figure; plot(corrBin); title('corr bin');
    figure; plot(numSpots); title('total spots');

%------------------------------------------------------------
% %   Visualization
    figure('color', 'w'); 
    plot(numSpots, 'k'); 
    hold on;
    findpeaks(numSpots); 
    [pkt,lct] = findpeaks(numSpots);
    plot(lct,pkt,'o','MarkerSize',5);
    ylabel('#spots per nucleus');
    xlabel('Threshold iterations');
    hold on; 
    yyaxis right;
%     plot(ssimBin, 'k'); findpeaks(ssimBin); ylabel('ssim spot pix');
    plot(corrBin, 'k'); findpeaks(corrBin); ylabel('Corr coeff. spot pix');
    ylim([0 1]);    
    ax = gca;
    ax = gca;
    ax.FontSize = 10;
    ax.LineWidth = 1;
    box(ax,'on');
    grid off;
    x0 = 100;
    y0= 100;
    plotWidth=250;
    plotHeight=250;
    set(gcf,'position',[x0,y0,plotWidth,plotHeight])
    hold off;
%------------------------------------------------------------
    
%     [~, thresInd] = max(corrVal(2:end));    

    thresInd = maxIter;    %%%%%%%%%%%%% change here %%%%%%%%%
    
%     jj = fix(thresInd/maxIter)+1;
%     if jj==0
%         jj=1;
%     end
%     if rem(thresInd,  maxIter) >0
%         pp = rem(thresInd,  maxIter);
%     else 
%         pp = 2;
%     end
% imThresLoc = imThresLoc + imTempStack{jj}(:,:,pp);
    imThresLoc = imThresLoc + imTempStack{j}(:,:,thresInd);
    figure; imshow(imcrop(imThresLoc, [min(col), min(row), (max(col) - min(col)), (max(row) - min(row))]),[]);
    totalSpots(1, i) = numSpots(thresInd);
end
end

%_____________________________________________________________
%   Plots pixel distribution and marks the bins of threshold
function plotImgHist(im, thres)
im = rescale(im);
[p1, e1] = histcounts(im, 'Normalization','probability');
binC = e1 + (e1(2) - e1(1))/2;
plot(binC(1:end-1), p1, 'r-');
hold on;
H1=area(binC(1:end-1), p1);
set(H1(1),'FaceColor',[0.9 0.9 0.7]);
hold on;
% t1 = graythresh(im);
% idx=binC<t1;
% Ht1=area(binC(idx),p1(idx));
% set(Ht1(1),'FaceColor',[0.7 0.7 0.5]);
for i = 1:length(thres)-1
    idx= find(binC>thres(i) & binC<thres(i+1));
    Ht1=area(binC(idx),p1(idx));
    
    if ~isempty(Ht1)
        set(Ht1(1),'FaceColor',[1/i 0.2 0.4]);
    end
    hold on;
end
xlabel('Normalized counts')
ylabel('Probability density')
title('Multilevel Otsu thresholding')
axis('square')
xlim([0 1]);
% legend('nucleus', 'background');
hold on;
end

function noisePlotter(imAvTemp, imStdTemp)
meanArr = reshape(imAvTemp, [numel(imAvTemp), 1]);
varArr = reshape(imStdTemp.^2, [numel(imStdTemp), 1]);
scatter(meanArr, varArr, '.r');
[sortMean,indSort] = sort(meanArr);
sortVar = varArr(indSort);
[lB, uB]= bounds(meanArr(2:end));
binMeanInd = discretize(sortMean,linspace(lB, uB, 10));
end