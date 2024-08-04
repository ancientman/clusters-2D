function startAnalysis(metaDataDS)
RawFolder  = metaDataDS.expInfo.rawImagesFolderName;
ProcFolder = metaDataDS.expInfo.procImagesFolderName;
XpixelSize = metaDataDS.imagingInfo.XpixelSize/1000;% changed to microns
metaDataDS.imagingInfo.XpixelSize = XpixelSize;%%%%%%%% Assignment here
YpixelSize = metaDataDS.imagingInfo.YpixelSize/1000;% changed to microns
metaDataDS.imagingInfo.YpixelSize = YpixelSize;%%%%%%%% Assignment here
bitDepth = metaDataDS.imagingInfo.bitDepth;
XYpsf = metaDataDS.imagingInfo.XYpsf/1000;% changed to microns
Zpsf = metaDataDS.imagingInfo.Zpsf/1000;% changed to microns
imUseType = metaDataDS.imagingInfo.imUseType;
overrideTimePoints = metaDataDS.imagingInfo.overrideTimePoints;
elementSize = metaDataDS.imagingInfo.elementSize;
clearBorder = metaDataDS.imagingInfo.clearBorder;
spotDetect = metaDataDS.imagingInfo.spotDetect;
spotFilter = metaDataDS.imagingInfo.spotFilter;
nucleusDetectMode = metaDataDS.imagingInfo.nucleusDetectMode;
nucThres = metaDataDS.imagingInfo.nucThres;
nucleusFeatureSize = metaDataDS.imagingInfo.nucleusFeatureSize; % already in microns
minSpotSize = metaDataDS.imagingInfo.minSpotSize/1000; % changed to microns from nm
centroidShiftConstraint = metaDataDS.imagingInfo.centroidShiftConstraint; % permissible shift between centroids in consecutive time points
if overrideTimePoints
    totalFrames = overrideTimePoints;
else
    totalFrames = metaDataDS.imagingInfo.T;
end
rawImages = dir(append(RawFolder,filesep,'*.tif'));
% Movie size
imDim = [metaDataDS.imagingInfo.Y,metaDataDS.imagingInfo.X,...
   totalFrames];
DeltaT  = metaDataDS.imagingInfo.DeltaT;
smoothFilterType = metaDataDS.imagingInfo.nucSmoothFilter;
if smoothFilterType == 7
    warning ('using abandoned filter type: DoG');
end
smoothingParam = metaDataDS.imagingInfo.smoothingParam;
deconvolutionType = metaDataDS.imagingInfo.deconvolutionType;
T = imDim(3);
imRawSmooth = zeros(imDim(1),imDim(2), T,'single'); %custom smoothing of raw frames
imRawSharp = zeros(imDim(1),imDim(2), T,'single'); %custom sharpening of raw frames
imRaw = zeros(imDim(1),imDim(2), T,'single'); %frames read from the raw data
imUse = zeros(imDim(1),imDim(2), T,'single'); %frames for use
imThres = zeros(imDim(1),imDim(2), T,'single'); %frames for use
nucMask = zeros(imDim(1),imDim(2), T, 'single'); %binary image with nucleus pixels
nucLabel = zeros(imDim(1),imDim(2), T, 'single'); %labeled image o nucleus pixels
spotMask = zeros(imDim(1),imDim(2), T, 'single'); %masks on hub spots
nucSpotLabel = zeros(imDim(1),imDim(2), T, 'single'); %labels on hub spots, each label for a nucleus
overlaySpotNuc = zeros(imDim(1),imDim(2), T, 'single'); %Visualization Output: nuc masks plus hub spot masks
overlaySpotPerim = zeros(imDim(1),imDim(2), T, 'single'); %Visualization Output: spot masks plus spots
onlySpot = zeros(imDim(1),imDim(2), T, 'single'); %Visualization Output: spots grayscale only
nucMeanRel = cell(1,T); %mean normalized pixel intensity value of each labeled nucleus per frame
nucMeanAbs = cell(1,T); %mean absolute pixel intensity value of each labeled nucleus per frame
nucMax = cell(1,T); %max pixel intensity value of each labeled nucleus per frame
nucCentroid = cell (1,T); %centroid of each labeled nucleus per frame
nucArea = cell(1,T); %area of each labeled nucleus per frame
framePropStruct = struct([]); %properties of connected components of each binarized frame
primeNucPropStruct = struct([]); %properties of labeled nuclei
alignedNucPropStruct = struct([]); %properties of aligned nuclei
nucPropStruct = struct([]); %properties of all nuclei
nucCentroidNL = {}; %centroid coordinates of all labeled nuclei
tempSpotPropStruct = struct([]); %properties of detected spots, overwritten every frame
spotPropStruct = struct([]); %properties of detected spots
spotIntensityNL = {}; %mean pixel intensity value of spots NL = labeled by nucleus
spotIntensityTL = cell (1,T); %mean pixel intensity value of spots NL = labeled by time
spotAreaNL = {}; %mean area of spots NL = labeled by nucleus
spotAreaTL = cell (1,T); %mean area of spots NL = labeled by time
spotCentroidNL = {}; %centroid coordinates of all spots within a labeled nucleus
relSpotCentroidNL = {}; %spot centroid minus nuc centroid of all spots within a labeled nucleus
saveTiffOptions.append = 1; %options for saveTiff function
saveTiffOptions.message = 0; %options for saveTiff function
stopMark = 0; %for in loop debugging
finalTimePoint = T; %for loop breaking
minNucSize = ceil(nucleusFeatureSize^2/(XpixelSize*YpixelSize)); %all in microns
metaDataDS.analysisInfo.minNucSize = minNucSize; %%%%%%%% Assignment here
removeNucMask = zeros(imDim(1),imDim(2));
if minSpotSize==0 || minSpotSize<XYpsf% Convert to pixels here
    minSpotSize = ceil((pi*(XYpsf/2)^2)/(XpixelSize*YpixelSize));
else
    minSpotSize = ceil((pi*(minSpotSize/2)^2)/(XpixelSize*YpixelSize));
end
metaDataDS.analysisInfo.minSpotSize = minSpotSize; %%%%%%%% Assignment here
if mod(ceil(sqrt(minSpotSize)),2)==1
    spotElementSize = ceil(sqrt(minSpotSize));
else
    spotElementSize = ceil(sqrt(minSpotSize))+1;
end
metaDataDS.analysisInfo.spotElementSize = spotElementSize; %%%%%%%% Assignment here

plotOn = 'no'; %   Options: 'yes' to turn on plots, 'no' to keep them off.
corrT = [];
autoCorrStruct = [];
pairCorrStruct = [];

corrT = cell(1, totalFrames);
%% Main image processing loop
for t=1:totalFrames
    fprintf('time iteration = %d of %d\n',t,totalFrames);
    imRaw(:,:,t) = (imread(append(RawFolder, filesep, rawImages(t).name)));
    if imUseType==1 % use "raw"
        imUse(:,:,t) = imRaw(:,:,t);
    elseif imUseType==2 % use "smooth"
        imUse(:,:,t) = Preprocess.smoothRawImage(imRaw(:,:,t), smoothFilterType, smoothingParam);
    elseif imUseType==3 % use "sharp"
        sharpParam = ceil(XYpsf/(2.0*XpixelSize));
        imUse(:,:,t) = Preprocess.sharpRawImage(imRaw(:,:,t), deconvolutionType, sharpParam);
    else % use "raw"
        imUse(:,:,t) = imRaw(:,:,t);
    end    
    %% Remove undesired nuclei by drawing polygons around them
%     if clearBorder == 0
        [imNucRemoved, removeNucMask] = Nucleus.removeNuc(imUse(:,:,t), t, removeNucMask);
        imUse(:,:,t) = imNucRemoved;
%     end    
    %% Segmentation and labeling of nuclei
    nucMask(:,:,t) = Nucleus.nucDetect(imUse(:,:,t), nucThres, nucleusDetectMode, ...
        elementSize, minNucSize, clearBorder, t, stopMark);
    if(isempty(nucMask(:,:,t)))
        fprintf('\nNo nuc found. breaking at t = %d\n', t);
        break;
    end
    [framePropStruct{t}] = Nucleus.nucLabeler(nucMask(:,:,t), imUse(:,:,t));
    %% Assigning properties of the labeled nuclei
    if t==1
        [primeNucPropStruct{t}] = Nucleus.nucProp(framePropStruct{t}, imUse(:,:,t), t, stopMark);
        nucLabel(:,:,t) = primeNucPropStruct{t}.labelMat;
        nucMeanRel{t} = primeNucPropStruct{t}.meanIntensity;
        nucMeanAbs{t} = primeNucPropStruct{t}.meanAbsIntensity;
        nucMax{t} = primeNucPropStruct{t}.maxIntensity;
        nucCentroid{t} = primeNucPropStruct{t}.centroid;
        nucArea{t} = primeNucPropStruct{t}.area;
        %% Break if the total nuclei in the current frames doesn't match the previous
    else
        if framePropStruct{t}.NumObjects ~= framePropStruct{t-1}.NumObjects
            fprintf('\nbreaking at t = %d\n', t);
            figure; imshowpair(label2rgb(labelmatrix(framePropStruct{t})),...
                label2rgb(labelmatrix(framePropStruct{t-1})), 'montage');
            caption = sprintf('image #%d and #%d', t, t-1);
            title(caption, 'FontSize', 14);
            break;
        else
            [alignedNucPropStruct{t}] = Nucleus.alignedNucProp(framePropStruct{t}, ...
                nucPropStruct{t-1}, imUse(:,:,t), t, stopMark);
            nucLabel(:,:,t) = alignedNucPropStruct{t}.labelMat;
            nucMeanRel{t} = alignedNucPropStruct{t}.meanRelIntensity;
            nucMeanAbs{t} = alignedNucPropStruct{t}.meanAbsIntensity;
            nucMax{t} = alignedNucPropStruct{t}.maxIntensity;
            nucCentroid{t} = alignedNucPropStruct{t}.centroid;
            nucArea{t} = alignedNucPropStruct{t}.area;
        end
    end
    if t==1
        plotOn = "no";
    end
    [corrT{t}] = Spots.crossCorr(imUse(:,:,t), nucMask(:,:,t), metaDataDS, plotOn);
    plotOn = "no";
    %   All properties of the nucleus labeled images
    %______________________________________________________________________________________
    nucPropStruct{t}.labelMat = nucLabel(:,:,t);
    nucPropStruct{t}.meanRelIntensity = nucMeanRel{t};
    nucPropStruct{t}.meanAbsIntensity = nucMeanAbs{t};
    nucPropStruct{t}.maxIntensity = nucMax{t};
    nucPropStruct{t}.centroid = nucCentroid{t};
    nucPropStruct{t}.area = nucArea{t};
    nNuc = max(nucPropStruct{t}.labelMat, [], 'all');
    
    %%   Multiple thresholding    
    %______________________________________________________________________________________
    [imThres(:,:,t), totalSpots(t,:)] = Nucleus.adapThres(imUse(:,:,t), nucLabel(:,:,t));
%     [imThres(:,:,t), totalSpots(t,:)] = Nucleus.adapThresVis(imUse(:,:,t), nucLabel(:,:,t));
%     imUse(:,:,t) = imThres(:,:,t);
    
    %% Segment and labeled spots. (all spots within a nucleus are labeled the same)
%     [spotMask(:,:,t), nucSpotLabel(:,:,t)] = Spots.spotDetect(imUse(:,:,t), ...
%         nucLabel(:,:,t), spotElementSize, spotDetect, spotFilter, minSpotSize);
%     [tempSpotPropStruct{t}.sProp] = Spots.spotProp(imUse(:,:,t), ...
%         nucSpotLabel(:,:,t), nNuc, t, stopMark);
%______________________________________________________________________________________
    [spotMask(:,:,t), nucSpotLabel(:,:,t)] = Spots.spotDetect(imThres(:,:,t), ...
        nucLabel(:,:,t), spotElementSize, 5, spotFilter, minSpotSize);
    [tempSpotPropStruct{t}.sProp] = Spots.spotProp(imUse(:,:,t), ...
        nucSpotLabel(:,:,t), nNuc, t, stopMark);
%______________________________________________________________________________________

    %% Assign the properties of the labeled spots
    for i=1:nNuc
        spotPropStruct{t}.spotRelMean{i} = (tempSpotPropStruct{t}.sProp{i}.meanSpotRelIntensity - nucPropStruct{t}.meanRelIntensity(i))';
        spotPropStruct{t}.spotAbsMean{i} = (tempSpotPropStruct{t}.sProp{i}.meanSpotAbsIntensity - nucPropStruct{t}.meanAbsIntensity(i))';
        spotPropStruct{t}.spotAbsTotVal{i} = tempSpotPropStruct{t}.sProp{i}.sumSpotAbsIntensity';
        spotPropStruct{t}.spotMols{i} = tempSpotPropStruct{t}.sProp{i}.spotMols;
        spotPropStruct{t}.spotArea{i} = tempSpotPropStruct{t}.sProp{i}.spotArea;
        spotPropStruct{t}.spotCentroid{i} = tempSpotPropStruct{t}.sProp{i}.spotCentroid;
        spotPropStruct{t}.spotDev{i} = (tempSpotPropStruct{t}.sProp{i}.stDevSpotIntensity)';
        spotPropStruct{t}.spotIdx{i} = (tempSpotPropStruct{t}.sProp{i}.spotIdx);
    end
   overlaySpotPerim(:,:,t) = rgb2gray(imoverlay(rescale(imUse(:,:,t)), bwperim(spotMask(:,:,t))));
   onlySpotTemp = imUse(:,:,t);
   onlySpotTemp(~spotMask(:,:,t)) = 0;
   onlySpot(:,:,t) = onlySpotTemp;
end
%%   Main loop ends
%_________________________________________________________________________________________________
%
if totalFrames>t
    breakingTimeFrame = t-1;
else
    breakingTimeFrame = totalFrames;
end
metaDataDS.analysisInfo.breakingTimeFrame = breakingTimeFrame; 
%_________________________________________________________________________________________________

corrXAll = [];
corrYAll = [];
for t = 1:breakingTimeFrame
    corrXAll = horzcat(corrXAll, horzcat(corrT{t}.corrX{:}));
    corrYAll = horzcat(corrYAll, horzcat(corrT{t}.corrX{:}));    
end
crossCorrStruct.corrXAll = corrXAll;
crossCorrStruct.corrYAll = corrYAll;
corrAll = horzcat(corrXAll, corrYAll);
meanCorrAll = mean(corrAll, 2, 'omitnan');
stdCorrAll = std(corrAll, 0, 2, 'omitnan');
semCorrAll = stdCorrAll./sqrt(size(corrAll,2));
figure('Color', 'w')
eb = errorbar(XpixelSize.*(1:size(semCorrAll, 1)), meanCorrAll, stdCorrAll, 'Color', 'k');
eb.CapSize = 0;
hold on; plot(XpixelSize.*(1:size(semCorrAll, 1)), zeros(1,size(semCorrAll, 1)), 'Color',[0.6 0.6 0.6], 'LineStyle',':');
ylabel('Correlation of pixel values'); xlabel('lag ({\mu}m)');
ax = axis;
x0 = 200;
y0 = 200;
plotWidth=250;
plotHeight=250;
box("off");
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);

% for i=1:breakingTimeFrame
%     
%     ssi(i) = ssim(imThres(:,:,1), imThres(:,:,i));
%     aCorrNuc(i) = corr2(imThres(:,:,1), imThres(:,:,i));
% end
% figure; plot(((1:breakingTimeFrame).*DeltaT), ssi); title('ssim vs time'); xlabel('Time(s)'); ylabel('ssim index');
% figure('Color', 'w'); plot(DeltaT*(1:totalFrames), totalSpots); title('Total spots vs time'); xlabel('Time(s)'); ylabel('Total spots');

%% Postprocessing
%_________________________________________________________________________________________________
[registerNucPropStruct, registeredFrames] = Nucleus.registeredNucProp(...
    nucPropStruct, imUse, centroidShiftConstraint, breakingTimeFrame);

[imNormAdd, imSpotAdd, imSpotMaskAdd] = Nucleus.globalNormalizedIntensity(...
    registerNucPropStruct, imUse, onlySpot, breakingTimeFrame);

[globalHotSpotStruct] = Spots.globalHotSpotFinder(...
    imNormAdd, nucPropStruct, nucThres, elementSize, spotElementSize, ...
    nucleusDetectMode, minNucSize, minSpotSize, clearBorder);

imAddDS.imNormAdd = imNormAdd;
imAddDS.imNormAddMask = globalHotSpotStruct.globalNucLabel;
imAddDS.imSpotAdd = imSpotAdd;
imAddDS.imSpotMaskAdd = imSpotMaskAdd;

%_________________________________________________________________________________________________
% % time intengration spatial cross correlation

imTimeAddTemp = zeros(size(imUse, 2, 1));
for i = 1:breakingTimeFrame
    imTimeAddTemp = imadd(imTimeAddTemp, double(rescale(imUse(:,:,i))));
    maskTimeAddTemp = Nucleus.nucDetect(imTimeAddTemp, nucThres, nucleusDetectMode, ...
        elementSize, minNucSize, clearBorder, 1, stopMark);
    [intTimeCrossCorrStruct{i}] = Spots.crossCorr(imTimeAddTemp, single(maskTimeAddTemp), metaDataDS, plotOn);
end



%_________________________________________________________________________________________________

plotOn = "no";
[allTimeCrossCorrStruct] = Spots.crossCorr(imAddDS.imNormAdd, imAddDS.imNormAddMask, metaDataDS, plotOn);
plotOn = "no";
crossCorrStruct.allTimeCrossCorr = allTimeCrossCorrStruct;
crossCorrStruct.intTimeCrossCorr = intTimeCrossCorrStruct;
%_________________________________________________________________________________________________

save([ProcFolder,'/xCorrDS.mat'],'crossCorrStruct'); 


try 
    [globalRegSpotStruct, globalSpotHotSpotStruct] = Spots.regSpotProp(nucPropStruct, ...
        registerNucPropStruct, spotPropStruct, globalHotSpotStruct.hotspotBW, ...
        registeredFrames, 1, breakingTimeFrame);
catch
    globalRegSpotStruct = [];
    globalSpotHotSpotStruct = [];
end

try
    [hotSpotPropNL] = Spots.timeBinHotSpot(...
        metaDataDS, registerNucPropStruct, imUse, ...
        nucPropStruct, spotPropStruct, registeredFrames);
catch
    hotSpotPropNL = [];
end

% plotOn = 'yes'; %   Options: 'yes' to turn on plots, 'no' to keep them off.
% crossCorrStruct = [];
% [crossCorrStruct] = Spots.crossCorr(...
%     imUse, registerNucPropStruct, globalHotSpotStruct, ...
%     registeredFrames, breakingTimeFrame, metaDataDS, plotOn);

autoCorrStruct = [];
% [autoCorrStruct] = Spots.autoCorr(...
%     imUse, registerNucPropStruct, globalHotSpotStruct, ...
%     registeredFrames, breakingTimeFrame, metaDataDS, plotOn);

pairCorrStruct = [];
[pairCorrStruct] = Spots.pairCorr(...
    imUse, globalHotSpotStruct.globalNucLabel, globalRegSpotStruct, registerNucPropStruct, ...
    registeredFrames, breakingTimeFrame, metaDataDS, XpixelSize, plotOn);

%_________________________________________________________________________________________________
% [nucBleachStruct, spotHotspotBleachStruct] = Spots.intensityCompare(imUse, metaDataDS, breakingTimeFrame, hotSpotPropNL);






%% Output datastructure saving
%_____________________________________________________________________________________
% PlottingFunctions.saveSomeStuff(...
%     ProcFolder, ...
%     registerNucPropStruct, ...
%     globalRegSpotStruct, ...
%     globalHotSpotStruct, ...
%     globalSpotHotSpotStruct, ...
%     imAddDS, ...
%     hotSpotPropNL, ...
%     crossCorrStruct, ...
%     pairCorrStruct, ...
%     metaDataDS);
%______________________________________________________________________________________

%% Visualization
nNuc = length(globalRegSpotStruct{1}.spotCentroid);
figure('Color', 'k');
% imDisp1 = imfuse(rescale(imNormAdd), bwperim(globalHotSpotStruct.hotspotBW));
imDisp1 = imoverlay(rescale(imNormAdd), bwperim(globalHotSpotStruct.hotspotBW));
% imDisp1 = imoverlay(rescale(imUse(:,:,1)), bwperim(globalHotSpotStruct.hotspotBW));
[row, col] = find(imAddDS.imNormAddMask);

imDisp1Crop = imcrop(imDisp1, [min(col), min(row), ...
                (max(col) - min(col)), (max(row) - min(row))]);
imshow(rescale(imDisp1Crop));

imNormAddCrop = imcrop(imNormAdd, [min(col), min(row), ...
                (max(col) - min(col)), (max(row) - min(row))]);
hold on;

% B = bwboundaries(rescale(imbinarize(imNormAdd)),'noholes');
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
% end
axis square;
axis equal;
axis image;
set(gca,'LooseInset',get(gca,'TightInset'));
x0=10;
y0=10;
width=250;
height=250;
axis square;
axis equal;
axis image;
set(gcf,'position',[x0,y0,width,height]);
% set(gca, 'color', 'k');
set(gca,'visible','off');
set(gca, 'YDir','reverse')
figure('Color', 'w');
%%
hold on;
for i = 1:nNuc
%     cellfun(@(x) plot(x.spotCentroid{i}(:, 1), x.spotCentroid{i}(:, 2), 'w.'), globalRegSpotStruct);    
    for t = 1:length(globalRegSpotStruct)
        if ~isempty(globalRegSpotStruct{t}.centroidInHotspot{i})
            scat = scatter(globalRegSpotStruct{t}.centroidInHotspot{i}(:, 1), globalRegSpotStruct{t}.centroidInHotspot{i}(:, 2));
            scat.Marker = 'o';
            scat.CData = [1 0 1 ];
%           scat.CData = [1 1 1];
            scat.MarkerFaceColor = scat.CData;
            scat.SizeData = 2;
            scat.MarkerFaceAlpha = 1;
            scat.MarkerEdgeAlpha = 1;
            hold on;
        end
        if ~isempty(globalRegSpotStruct{t}.centroidOutHotspot{i})
            scat = scatter(globalRegSpotStruct{t}.centroidOutHotspot{i}(:, 1), globalRegSpotStruct{t}.centroidOutHotspot{i}(:, 2));
            scat.Marker = 'o';
            scat.CData = [0 1 1];
%             scat.CData = [1 1 1];
            scat.MarkerFaceColor = scat.CData;
            scat.SizeData = 2;
            scat.MarkerFaceAlpha = 1;
            scat.MarkerEdgeAlpha = 1;
            hold on;
        end
    end
    hold on;
end
B = bwboundaries(globalHotSpotStruct.hotspotBW,'noholes');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'color', [0.2 0.2 0.2], 'LineWidth', 2)
end
hold on;
B = bwboundaries(rescale(imbinarize(imNormAdd)),'noholes');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'color', [0.2 0.2 0.2], 'LineWidth', 2)
end
%%
x0=10;
y0=10;
width=250;
height=250;
axis square;
axis equal;
axis image;
set(gcf,'position',[x0,y0,width,height]);
% set(gca, 'color', 'k');
set(gca,'visible','off');
set(gca, 'YDir','reverse')
hold off;
%______________________________________________________________________________________

%% Write a movie
%  v = VideoWriter('G:\Tyrone_analysis\temp\newfile.avi','Uncompressed AVI');
%  v.FrameRate = 6;
%  open(v);
%  for k=1:size(onlySpot, 3)     
%      writeVideo(v,imoverlay(rescale(onlySpot(:,:,k)), bwperim(nucLabel(:,:,k))));
% %      writeVideo(v,imoverlay(rescale(onlySpot(:,:,k)), bwperim(globalHotSpotStruct.hotspotBW)));
%  end
%  close(v);
%______________________________________________________________________________________

% HelperFunctions.saveTiff(overlaySpotPerim, append(ProcFolder, filesep, 'spotStack', '.tif'), saveTiffOptions);
% HelperFunctions.saveTiff((bwperim(nucMask) | spotMask), append(ProcFolder, filesep, 'spotNucStack', '.tif'), saveTiffOptions);
% hotSpotStack = zeros([size(hotSpotPropNL.hotSpotBW{1}) size(hotSpotPropNL.hotSpotBW, 2)]);
% for i=1:size(hotSpotPropNL.hotSpotBW, 2)
%     hotSpotStack(:,:,i) = 256.*uint8(hotSpotPropNL.hotSpotBW{i});
% end
% % hotSpotStack = (cellfun(@(x) cat(3,x), hotSpotPropNL.hotSpotBW, 'UniformOutput', false));
% HelperFunctions.saveTiff(uint8(hotSpotStack), append(ProcFolder, filesep, 'hotSpotStack', '.tif'), saveTiffOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting functions
%______________________________________________________________________________________
if(breakingTimeFrame>10)% Only call plotting functions for data of more than 10 timepoints
%     PlottingFunctions.spotPropPlot1(ProcFolder); % nucleus intensity properties
    PlottingFunctions.spotPropPlot3(ProcFolder); % hotspot occupancy 
    PlottingFunctions.spotPropPlot4(ProcFolder); % spots occurance rates in and outside hotspots
    PlottingFunctions.spotPropPlot7(ProcFolder); % area and intensity of spots inside and outside hotspots
end
%______________________________________________________________________________________
end