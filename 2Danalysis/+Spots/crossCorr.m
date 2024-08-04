function [corrStruct] = crossCorr(im, imBin, metadataDS, plotOn)
if strcmp(plotOn, 'yes')
    plotOn = 1;
else
    plotOn = 0;
end

xPixSize = metadataDS.imagingInfo.XpixelSize;
yPixSize = metadataDS.imagingInfo.YpixelSize;

imBin = imbinarize(imBin);
imBin = bwareaopen(imBin, 900);
corrStruct = [];
if plotOn
    figure('Color', 'w')
    subplot(2,2,1);
    imshow(rescale(im));
end
cent = regionprops(imBin, 'Centroid');
cent = vertcat(cent(:).Centroid); % Cent: First element = x, second element = y
bounds = bwboundaries(imBin); % bound: 2nd element, x; 1st element, y.
bb = regionprops(imBin, "BoundingBox");
bb = vertcat(bb(:).BoundingBox);
% diagonal = sqrt(bb(:,3).^2 + bb(:,4).^2);
minD = min(bb(:,3), bb(:,4))./2;
line1Ends = horzcat(cent(:,1) + minD, cent(:,2) - minD, cent(:,1) - minD, cent(:,2) + minD); % [Top x, top y, bot x, bot y]
line2Ends = horzcat(cent(:,1) - minD, cent(:,2) - minD, cent(:,1) + minD, cent(:,2) + minD); % [Top x, top y, bot x, bot y]
% hold on; plot([line1Ends(1,1), line1Ends(1,3)], [line1Ends(1,2), line1Ends(1,3)], 'k--');
% hold on; plot([line2Ends(1,1), line2Ends(1,3)], [line2Ends(1,2), line2Ends(1,4)], 'k--');
line1 = cell(1, size(cent, 1));
line2 = cell(1, size(cent, 1));
lineSegIn1 = cell(1, size(cent, 1));
lineSegIn2 = cell(1, size(cent, 1));
pgonTemp = cell(1, size(cent, 1));
boxX = cell(1, size(cent, 1));
boxY = cell(1, size(cent, 1));
pgon = cell(1, size(cent, 1));
nucCrop = cell(1, size(cent, 1));
maskNucBox = cell(1, size(cent, 1));
maskAll = zeros(size(im, [1,2]));
corrX = cell(1, size(cent, 1));
corrY = cell(1, size(cent, 1));
xLags = cell(1, size(cent, 1));
yLags = cell(1, size(cent, 1));
meanCorrX = cell(1, size(cent, 1));
stdCorrX = cell(1, size(cent, 1));
meanCorrY = cell(1, size(cent, 1));
stdCorrY = cell(1, size(cent, 1));
% im = 2^8.*ones(size(im));
% im = imnoise(im, 'salt & pepper');
% im = uint16(im);
im = rescale(im);

for i = 1:size(cent, 1)
    line1{i} = horzcat(linspace(line1Ends(i,1),line1Ends(i,3), 100)', ...
        linspace(line1Ends(i,2),line1Ends(i,4), 100)');
    line2{i} = horzcat(linspace(line2Ends(i,1),line2Ends(i,3), 100)', ...
        linspace(line2Ends(i,2),line2Ends(i,4), 100)');
    pgonTemp{i} = polyshape(bounds{i}(:,2), bounds{i}(:,1));
    [lineSegIn1{i}, ~] = intersect(pgonTemp{i}, line1{i});
    [lineSegIn2{i}, ~] = intersect(pgonTemp{i}, line2{i});
    boxX{i} = [max(min(lineSegIn2{i}(:,1)), min(lineSegIn1{i}(:,1))), ...
        min(max(lineSegIn2{i}(:,1)), max(lineSegIn1{i}(:,1)))];
    boxY{i} = [max(min(lineSegIn2{i}(:,2)), min(lineSegIn1{i}(:,2))), ...
        min(max(lineSegIn2{i}(:,2)), max(lineSegIn1{i}(:,2)))];    
    pgon{i} = polyshape([horzcat(boxX{i}, fliplr(boxX{i}))], repelem(boxY{i},2));
    maskNucBox{i} = poly2mask([horzcat(boxX{i}, fliplr(boxX{i}))], repelem(boxY{i},2), size(im,1), size(im,2));
    maskAll = or(maskAll, maskNucBox{i});     

    nucCrop{i} = imcrop(im, [min(boxX{i}), min(boxY{i}), ...
        max(boxX{i}) - min(boxX{i})+1, max(boxY{i}) - min(boxY{i})+1]);
    [corrX{i}, xLags{i}, corrY{i}, yLags{i}] = corrFunc(nucCrop{i});
    meanCorrX{i} = mean(corrX{i}, 1)';    
    stdCorrX{i} = std(corrX{i}, 0, 1)';
    meanCorrX{i} = meanCorrX{i}(((length(meanCorrX{i})+1)/2):end);
    meanCorrY{i} = mean(corrY{i}, 1);
    stdCorrY{i} = std(corrX{i}, 0, 1);
    meanCorrY{i} = meanCorrY{i}(((length(meanCorrY{i})+1)/2):end)';
    if plotOn
        subplot(2,2,2); hold on;
        scatter(bounds{i}(:,2), bounds{i}(:,1), 'b.');
        set(gca, 'YDir','reverse')
        hold on;
        plot(cent(i,1), cent(i,2), '*r');
        hold on; plot(line2{i}(:,1), line2{i}(:,2), 'm--')
        hold on; plot(line1{i}(:,1), line1{i}(:,2), 'm--')
        hold on; plot(pgon{i})
        hold on;
        axis equal;
        axis square;
        subplot(2,2,3); hold on;
        imshow(double(maskAll).*double(im), []);
    end
end

if plotOn
    hold off;
end
corrXLim = min(cellfun(@length, meanCorrX));
corrYLim = min(cellfun(@length, meanCorrY));
allCorrX = cellfun(@(x)x(1:corrXLim), meanCorrX, 'un', 0);
meanAllCorrX = mean(horzcat(allCorrX{:}),2);
allCorrY = cellfun(@(x)x(1:corrYLim), meanCorrY, 'un', 0);
meanAllCorrY = mean(horzcat(allCorrY{:}),2);

meanAllCorrX = meanAllCorrX(1: min(length(meanAllCorrX), length(meanAllCorrY)));
meanAllCorrY = meanAllCorrY(1: min(length(meanAllCorrX), length(meanAllCorrY)));

allCorr = mean(horzcat(meanAllCorrX, meanAllCorrY), 2);
semCorr = std(horzcat(meanAllCorrX, meanAllCorrY), 0, 2)./sqrt(size(horzcat(meanAllCorrX, meanAllCorrY),2));

if plotOn
    subplot(2,2,4);
    % plot(0.043*(1:length(meanAllCorr)), meanAllCorrX); hold on; plot(0.043*(1:length(meanAllCorrY)), meanAllCorrY);
    % plot(0.043*(1:length(allCorr)), allCorr); 
    hold on; 
    eb = errorbar(xPixSize*(1:length(allCorr)), allCorr, semCorr, 'Color', 'k');
    eb.CapSize = 0;
    hold on; plot(xPixSize*(1:length(allCorr)), zeros(1,length(allCorr)), 'Color',[0.6 0.6 0.6], 'LineStyle',':');
    ylabel('Correlation of pixel values'); xlabel('lag ({\mu}m)');    
end
corrStruct.corrX = allCorrX;
corrStruct.corrY = allCorrY;
end

% corrFunc()
% 
% imBoxNuc = im;
% imBoxNuc(~maskAll) = 0;
% 
% 
% imshow(imBin);
% imSeg = im;
% imSeg(imBin==0) = NaN;
% subplot(2,2,4);
% imshow(rescale(imSeg));

function [cx, lagx, cy, lagy] = corrFunc(im)

numLags = min(length(im)-1, 30);
% cx = zeros(size(im, 1), 2*size(im, 2)-1);
% lagx = zeros(size(im, 1), 2*size(im, 2)-1);
% cy = zeros(size(im, 2), 2*size(im, 1)-1);
% lagy = zeros(size(im, 2), 2*size(im, 1)-1);

% x correlation
for i = 1:size(im,1)
    x = double(im(i,:));
%     [cx(i,:),lagx(i,:)] = crosscorr(x, x, NumLags=60);
    [cx(i,:),lagx(i,:)] = crosscorr(x, x, NumLags=numLags);
end
for i = 1:size(im,2)
    y = double(im(:,i));
% [cy(i,:),lagy(i,:)] = crosscorr(y, y, NumLags=60);
    [cy(i,:),lagy(i,:)] = crosscorr(y, y, NumLags=numLags);
end
% errorbar(0.043.*lagy(1,:), mean(cy,1), std(cy,0,1)./sqrt(size(cy,1)));
% hold on;
% errorbar(0.043.*lagx(1,:), mean(cx,1), std(cx,0,1)./sqrt(size(cx,1)));
% xlim([0, inf]);
end

function [corrStruct] = crossCorrOld(imUse, registerNucPropStruct, globalHotSpotStruct, registeredFrames, breakingTimeFrame, metadataDS, plotOn)
im = imUse;
hsDS = globalHotSpotStruct;
nucDS = registerNucPropStruct;
frames = registeredFrames;
tFinal = breakingTimeFrame;
delT = metadataDS.imagingInfo.timeResolution; %in seconds
%~~~~~~~~~~~~~~~~~~~
inParam.wholeFrame = 'no'; %   'Options: 'yes', 'no'. Yes is for whole frame
inParam.wholeNuc = 'yes'; % 'Options: 'yes', 'no'. Yes is for whole nuc vs background. use only imclearborder.
inParam.window = 12; % Hardcoded (about one micron): Actual window size is 2*window + 1
inParam.hotSpotBB = 'off'; % 'on' for hotspot bb, 'off' for fixes window around centroid.
inParam.filterType = 'detrend'; %   Oprions 'none', 'medfilt', 'detrend'

%________________________________________________________________________________
% whole frame autocorrelation
corrCoeffFrame = zeros(numel(frames),1);
meanFrame = zeros(numel(frames),1);
    t = 1;
    while t <=numel(frames)
        ff = frames(t);        
        meanFrame(t) = mean(im(:,:,t),'all');
        corrCoeffFrame(t) = corr2(im(:,:,1), im(:,:,ff));
        t = t+1;
    end
    timeFrame = frames'.*delT;
    fitFrame = fit(timeFrame, meanFrame, 'exp1');  
    
    corrStruct.frameMean = meanFrame;
    corrStruct.frameCorr = corrCoeffFrame;
    corrStruct.timeFrame = timeFrame;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     figure('Color','w');
%     hold on;
%     [value1, edge1] = histcounts(im(:,:,1), 'Normalization', 'probability');
%     [valNuc1, edgeNuc1] = histcounts(imNuc1, 'Normalization', 'probability');
%     [valBg1, edgeBg1] = histcounts(imBg1, 'Normalization', 'probability');
%     [value2, edge2] = histcounts(stack(:,:,end), 'Normalization', 'probability');
%     [valNuc2, edgeNuc2] = histcounts(imNuc2, 'Normalization', 'probability');
%     [valBg2, edgeBg2] = histcounts(imBg2, 'Normalization', 'probability');
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     labelMatGlobal = imclearborder(globalHotSpotStruct.globalNucLabel);
%     labelMatGlobal = bwlabel(imbinarize(labelMatGlobal));        
%     bgMask = zeros(size(im(:,:,1)));
%     [~, bgMask] = Nucleus.removeNuc(im(:,:,1), 1,bgMask);

    nNuc = numel(unique(nucDS{1}.labelMat)) - 1; %  Use with imclearborder only
    meanNucWhole = cell(1, (nNuc));
    corrCoeffNucWhole = cell(1, (nNuc));        
    mockIm1 = zeros(size(im, 2), size(im, 1));
    for i=1:nNuc
        t = 1;
        imWin = cell(1, length(frames));
        figure('Color', 'w');
        while t <=length(frames)
            ff = frames(t);    
            labelMat = (nucDS{ff}.labelMat == i);
            imWin{ff} = im(:,:,ff);
            imWin{ff}(~labelMat) = 0;
            p1 = plot(histcounts(nonzeros(imWin{ff}), 'Normalization', 'probability'));
            if ff ==1
                p1.Color = [1, 0, 0];
            elseif ff == length(frames)
                p1.Color = [0, 0, 1];
            else
            p1.Color = [ff, ff, ff]/(2*length(frames));
            end
            hold on;
            t = t+1;
        end
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (strcmp(plotOn,'yes'))
        figure('Color', 'w')
        plot(timeFrame, corrCoeffFrame, 'k-', 'Linewidth', 1.5);
        xlabel('Time (s)');
        ylabel('{\sigma}');
        title('whole frame correlation');
        figure('Color', 'w');
        plot(fitFrame, timeFrame, meanFrame);
        xlabel('Time (s)');
        ylabel('Fluorescence intensity (a. u)');
        title('whole frame bleaching');
    end
end