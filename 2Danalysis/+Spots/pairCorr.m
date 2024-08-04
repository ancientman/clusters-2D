function [pairCorrStruct] = pairCorr(im, globalNucLabel, globalRegSpotStruct, registerNucPropStruct, registeredFrames, breakingTimeFrame, metaDataDS, pixSize, plotOn)
prunePix = 3;
rMax = 150;
cropZeroPad = 100; 
nNuc = length(registerNucPropStruct{1}.pixList);
tTime = length(globalRegSpotStruct);
imCentIn = cell(1, nNuc);
imCentOut = cell(1, nNuc);
imSpotIn = cell(1, nNuc);
imSpotOut = cell(1, nNuc);
imSpotBinIn = cell(1, nNuc);
imSpotBinOut = cell(1, nNuc);
corrIn = cell(1, nNuc);
corrOut = cell(1, nNuc);
corrSpotRadAvg = cell(1, nNuc);
corrSpotRadSem = cell(1, nNuc);
corrCentRadAvg = cell(1, nNuc);
corrCentRadSem = cell(1, nNuc);
lambdaHalf = cell(1, nNuc);
accumSpotCent = cell(1, nNuc);
centCorr = cell(1, nNuc);

%________________________________________________________
for i = 1:nNuc
    accumSpotCent{i} = zeros(size(im(1, 2)));
    for t = 1:tTime
        imTemp = im(:,:,t);
        tempCellIn = globalRegSpotStruct{t}.idxInHotspot{i};
        if ~isempty(globalRegSpotStruct{t}.idxInHotspot{i})
            tempCellIn = tempCellIn(~(cellfun(@isempty, tempCellIn)));
        else
            tempCellIn = {};
        end        
        
        tempCellOut = globalRegSpotStruct{t}.idxOutHotspot{i};
        if ~isempty(globalRegSpotStruct{t}.idxOutHotspot{i})    
            tempCellOut = tempCellOut(~(cellfun(@isempty, tempCellOut)));
        else
            tempCellOut = {};
        end
        
        tempLengthCellIn = cellfun(@length , tempCellIn, 'un', 0);
        tempPrunedCellIn = ...
            tempCellIn(cellfun(@(x) any(x>prunePix, 'all'), tempLengthCellIn));
        tempLengthCellOut = cellfun(@length , tempCellOut, 'un', 0);
        tempPrunedCellOut = ...
            tempCellOut(cellfun(@(x) any(x>prunePix, 'all'), tempLengthCellOut));
        tempPixArrIn = vertcat(tempPrunedCellIn{:});
        tempPixArrOut = vertcat(tempPrunedCellOut{:});
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        imBuffSpot = zeros(size(imTemp));
        imBuffSpot(tempPixArrIn) = 1;
        imSpotBinIn{i}(:,:,t) = double(imBuffSpot);
        imSpotIn{i}(:,:,t) = imSpotBinIn{i}(:,:,t).*im(:,:,t);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        imBuffSpot = zeros(size(imTemp));
        imBuffSpot(tempPixArrOut) = 1;
        imSpotBinOut{i}(:,:,t) = double(imBuffSpot);
        imSpotOut{i}(:,:,t) = imSpotBinOut{i}(:,:,t).*im(:,:,t);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
                
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        imBuffCent = zeros(size(imTemp));
        if ~isempty(tempCellIn)
            tempCentInInd = sub2ind(size(imTemp), ...
                globalRegSpotStruct{t}.centroidInHotspot{i}(:, 2), ...
                globalRegSpotStruct{t}.centroidInHotspot{i}(:, 1));
            imBuffCent(tempCentInInd) = 1;
        end        
        imCentIn{i}(:,:,t) = double(imBuffCent);        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        imBuffCent = zeros(size(imTemp));
        if ~isempty(tempCellOut)
            tempCentOutInd = sub2ind(size(imTemp), ...
                globalRegSpotStruct{t}.centroidOutHotspot{i}(:, 2), ...
                globalRegSpotStruct{t}.centroidOutHotspot{i}(:, 1));
            imBuffCent(tempCentOutInd) = 1;
        end        
        imCentOut{i}(:,:,t) = double(imBuffCent);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sumAllPix = size(im, 1)*size(im, 2);
        sumSpotIn = length(tempPixArrIn);
        sumSpotOut = length(tempPixArrOut);
        
        corrIn{i}(:,:,t) = (sumAllPix^2/sumSpotIn^2)^(-1).*(xcorr2(imSpotIn{i}(:,:,t), imSpotIn{i}(:,:,1)));
        corrOut{i}(:,:,t) = (sumAllPix^2/sumSpotOut^2)^(-1).*(xcorr2(imSpotOut{i}(:,:,t), imSpotOut{i}(:,:,1)));
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        corrCut = imcrop(corrOut{i}(:,:,t), [floor(size(im, 2))-floor(rMax), floor(size(im, 2))-floor(rMax), 2*rMax, 2*rMax]); 
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        [valAvg, valSem] = cart2PolCorr(rMax, corrCut);
        corrSpotRadAvg{i}(:,t) = valAvg;
        corrSpotRadSem{i}(:,t) = valSem;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        accumSpotCent{i} = accumSpotCent{i} | imCentIn{i}(:,:,t);
        accumSpotCent{i} = accumSpotCent{i} | imCentOut{i}(:,:,t);
        
    end
    [row, col] = find(accumSpotCent{i});
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentNuc = globalNucLabel;
    currentNuc(currentNuc~=i) = 0;
    bounds = bwboundaries(currentNuc,'noholes');
    cent = regionprops(currentNuc, 'Centroid');
    cent = vertcat(cent(:).Centroid); % Cent: First element = x, second element = y
    bb = regionprops(currentNuc, "BoundingBox");
    bb = vertcat(bb(:).BoundingBox);
    minD = min(bb(:,3), bb(:,4))./2;
    line1Ends = horzcat(cent(:,1) + minD, cent(:,2) - minD, cent(:,1) - minD, cent(:,2) + minD); % [Top x, top y, bot x, bot y]
    line2Ends = horzcat(cent(:,1) - minD, cent(:,2) - minD, cent(:,1) + minD, cent(:,2) + minD); % [Top x, top y, bot x, bot y]

    line1 = horzcat(linspace(line1Ends(i,1),line1Ends(i,3), 100)', ...
        linspace(line1Ends(i,2),line1Ends(i,4), 100)');
    line2 = horzcat(linspace(line2Ends(i,1),line2Ends(i,3), 100)', ...
        linspace(line2Ends(i,2),line2Ends(i,4), 100)');
    pgonTemp = polyshape(bounds{1}(:,2), bounds{1}(:,1));
    [lineSegIn1, ~] = intersect(pgonTemp, line1);
    [lineSegIn2, ~] = intersect(pgonTemp, line2);
    boxX = [max(min(lineSegIn2(:,1)), min(lineSegIn1(:,1))), ...
        min(max(lineSegIn2(:,1)), max(lineSegIn1(:,1)))];
    boxY = [max(min(lineSegIn2(:,2)), min(lineSegIn1(:,2))), ...
        min(max(lineSegIn2(:,2)), max(lineSegIn1(:,2)))];    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure; 
    imshow(accumSpotCent{i}); 
    hold on; 
    plot(bounds{1}(:,2), bounds{1}(:,1), 'g', 'LineWidth', 1)
    hold on;
    pgon = polyshape([horzcat(boxX, fliplr(boxX))], repelem(boxY,2));
    hold on; plot(pgon); hold on;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
%     plot([min(col), max(col), max(col), min(col), min(col)], ...
%         [max(row), max(row), min(row), min(row), max(row)]);
%     spotCentCrop = imcrop(accumSpotCent{i}, [min(col), min(row), (max(col) - min(col)), (max(row) - min(row))]);
    spotCentCrop = imcrop(accumSpotCent{i}, [min(boxX), min(boxY), ...
        max(boxX) - min(boxX)+1, max(boxY) - min(boxY)+1]);
    centCropPad = zeros(size(spotCentCrop, 1) + 2*cropZeroPad,  size(spotCentCrop, 2) + 2*cropZeroPad);
    centCropPad(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = spotCentCrop;
    mask = zeros(size(centCropPad));
    mask(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = 1;
    nParticles = sum(spotCentCrop, 'all');  % number of particles within mask
    nPixMask = sum(mask, 'all');      % area of mask
    NP = real(fftshift(ifft2(abs(fft2(mask)).^2))); % Normalization for correct boundary conditions
    G1 = nPixMask^2/nParticles^2*real(fftshift(ifft2(abs(fft2(centCropPad)).^2)))./NP; % 2D G(r) with proper normalization
    G = imcrop(G1, [floor(size(G1, 2)/2+1)-cropZeroPad, floor(size(G1, 1)/2+1)-cropZeroPad, 2*cropZeroPad, 2*cropZeroPad]);  %only return valid part of G
    [valAvg, valSem] = cart2PolCorr(cropZeroPad, G);
    corrCentRadAvg{i} = valAvg;
    corrCentRadSem{i} = valSem;
    lambdaHalf{i} = plotErr(corrCentRadAvg{i}, corrCentRadSem{i}, pixSize);
end
pairCorrStruct.avg = corrCentRadAvg;
pairCorrStruct.sem = corrCentRadSem;
pairCorrStruct.lambdaHalf = lambdaHalf;
end

function [valRadBinAvg, valRadBinSem] = cart2PolCorr(rMax, corrCut)
xVals = ones(1, 2*rMax+1)'*(-rMax:rMax);
yVals = (-rMax:rMax)'*ones(1, 2*rMax+1);
[theta, rho, valTemp] = cart2pol(xVals,yVals,  corrCut);  % convert x, y to polar coordinates
aR = reshape(rho,1, (2*rMax+1)^2);
aValTemp = reshape(valTemp,1, (2*rMax+1)^2);
[rrTemp,ind] = sort(aR);
corrRadVal = aValTemp(ind);
rad= 0:floor(max(rrTemp));
[n, bin] = histcounts(rrTemp, rad);
valRadBinAvg = zeros(length(bin)-1, 1);
valRadBinSem = zeros(length(bin)-1, 1);
for j = 1:length(bin)-1
    valRadBinAvg(j) = sum(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sum(rrTemp>bin(j) & rrTemp<=bin(j+1));
    valRadBinSem(j) = std(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sqrt(sum(rrTemp>bin(j) & rrTemp<=bin(j+1)));                        
end
end

function halfLength = plotErr(val, err, pixSize)
halfLength = [];
xArr = (1:length(val)).*pixSize;
figure('color', 'w');

% ----------- Patch plot -----------------
% patchTop = val+err;
% patchTop = reshape(patchTop,1,[]);
% patchBot = val-err;
% patchBot = reshape(patchBot, 1, []);
% yPatch=[patchBot,fliplr(patchTop)];
% xPatch=[xArr,fliplr(xArr)];
% pt = patch(xPatch, yPatch, 1);
% pt.FaceColor = [0.6 0.6 0.6];
% pt.EdgeColor = 'none';
% pt.FaceAlpha = 0.6;
% hold on;
% pl = plot(xArr, val);
% pl.LineWidth = 1;
% pl.Color = [0.4 0.4 0.4];
% pl.LineStyle = '-';
% hold on;
% ----------------------------------------

% ----------- Err plot -----------------
pt = errorbar(xArr, val, err);
pt.Marker = 'none';
pt.MarkerSize = 1;
pt.LineWidth = 1;
pt.Color = [0 0 0];
pt.LineStyle = '-';
pt.CapSize = 1;
hold on;
% ----------------------------------------

% ----------- Fit -----------------
[~, limMin] = max(val(1:40));
lims = limMin:40;
tbl = table(xArr(lims)', val(lims));
modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));  
beta0 = [0.2, 5, 0.5]; 
try
    mdl = fitnlm(tbl, modelfun, beta0);
    coeffs= mdl.Coefficients{:, 'Estimate'};
    yFitted = coeffs(1) + coeffs(2) * exp(-coeffs(3)*xArr);
    halfLength = xArr(3)+ log(2)/coeffs(3);
    hold on;
    pf = plot(xArr, yFitted, 'k-', 'LineWidth', 1);
    hold on;
catch
    aaa = 1
end
% ---------------------------------

% ----------- Line -----------------

pp = plot([0, length(val)].*pixSize, [1 1]);
pp.LineWidth = 1;
pp.Color = [0.8 0 0];
pp.LineStyle = '--';
% ----------------------------------

ylabel('Pair correlation G(r)');
xlabel('r ({\mu}m)');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
xlim([0 3]);
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end