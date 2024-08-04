function zStacker()
filePath = 'E:\Analysis\input\0628_bcdXhb5_nc14_em3_4.tif';
stackInfo = imfinfo(filePath);
imFrame1 = imread(filePath, 1);
% stackInfo = imfinfo('C:\xData\Rawdata_RM\bcdxhb\imStack14.tif');
% imFrame1 = imread('C:\xData\Rawdata_RM\bcdxhb\imStack14.tif', 1);
totalSlice = size(stackInfo, 1);
totalTimePoint = 8;
channel = 2;
slice = 9;
firstChannel = 'green';
secondChannel = 'red';
greenStack = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
reScaleGreenStack = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
maxGreenProject = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
greenProjNucRaw = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
greenProjNucFilt = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwGreenProjNucRaw = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwGreenProjNucFilt = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);

redStack = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
reScaleRedStack = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
minRedProject = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
redNucProjRaw = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
redNucProjFilt = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwRedNucProjRaw = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwRedNucProjFilt = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwBothProjNucRaw = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
bwBothProjNucFilt = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);

greenSpotMask = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
greenSpotNL = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
greenSpotProject = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);

redSpotMask = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
redSpotNL = zeros([size(imFrame1, 2), size(imFrame1, 1), slice, totalTimePoint]);
redSpotProject = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);

bothSpotMask = zeros([size(imFrame1, 2), size(imFrame1, 1), 3, slice, totalTimePoint]);

removeNucMask = zeros([size(imFrame1, 2), size(imFrame1, 1)]);

projNucMask = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
convHullProjNuc = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
greenK = 1;
redK = 2;

projNucLabelMat = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
hullLabelMat = zeros([size(imFrame1, 2), size(imFrame1, 1), totalTimePoint]);
projNucCent = cell(1, totalTimePoint);
redCent = cell(1, totalTimePoint);
spotElementSize = 3;
spotDetect = 1;
spotFilter = 1;
minSpotSize = 20;
% if minSpotSize==0 || minSpotSize<XYpsf% Convert to pixels here
%     minSpotSize = ceil((pi*(XYpsf/2)^2)/(XpixelSize*YpixelSize));
% else
%     minSpotSize = ceil((pi*(minSpotSize/2)^2)/(XpixelSize*YpixelSize));
% end

for j=1:totalTimePoint
    fprintf('time point = %d\n', j);
    
    for i=1:slice
        greenStack(:,:,i, j) = imread(filePath, greenK);
        reScaleGreenStack(:,:,i,j) = rescale(greenStack(:,:,i, j));
        maxGreenProject(:,:,j) = max(reScaleGreenStack(:,:,:,j), [], 3);
        
        
        redStack(:,:,i, j) = imread(filePath, redK);
        reScaleRedStack(:,:,i,j) = rescale(redStack(:,:,i, j));
        minRedProject(:,:,j) = min(reScaleRedStack(:,:,:,j), [], 3);
                
        greenK = greenK + channel;
        redK = redK + channel;
    end
    
    greenProjNucRaw(:,:,j) = (imadjust(maxGreenProject(:,:,j)));
    bwGreenProjNucRaw(:,:,j) = imbinarize(greenProjNucRaw(:,:,j));
    
    greenProjNucFilt(:,:,j) = imgaussfilt(greenProjNucRaw(:,:,j));
    bwGreenProjNucFilt1 = imbinarize(greenProjNucFilt(:,:,j));
    bwGreenProjNucFilt2 = imfill(bwGreenProjNucFilt1, 'holes');
    bwGreenProjNucFilt(:,:,j) = (bwareaopen(bwGreenProjNucFilt2, 200)); 

    redNucProjRaw(:,:,j) = imcomplement(imadjust(minRedProject(:,:,j)));
    bwRedNucProjRaw1 = ~imbinarize(minRedProject(:,:,j));
    bwRedNucProjRaw2 = imfill(bwRedNucProjRaw1, 'holes');
    bwRedNucProjRaw2 = (bwareaopen(bwRedNucProjRaw2, 200));
    bwRedNucProjRaw(:,:,j) = bwRedNucProjRaw2;
    
    redNucProjFilt(:,:,j) = HelperFunctions.wienerFilter2D(redNucProjRaw(:,:,j), 9, 9);
    bwRedNucProjFilt1 = imbinarize(redNucProjFilt(:,:,j));
    bwRedNucProjFilt2 = imfill(bwRedNucProjFilt1,'holes');
    bwRedNucProjFilt(:,:,j) = (bwareaopen(bwRedNucProjFilt2, 200)); 
    
    bwBothProjNucRaw(:,:,j) = bwGreenProjNucRaw(:,:,j).*bwRedNucProjRaw(:,:,j);
    
    bwBothProjNucFilt(:,:,j) = bwGreenProjNucFilt(:,:,j).*bwRedNucProjFilt(:,:,j);
    
    %% Remove undesired nuclei by drawing polygons around them
    [imNucRemoved, removeNucMask] = Nucleus.removeNuc(bwBothProjNucRaw(:,:,j), j, removeNucMask);
    bw = imNucRemoved;    
%     bw = bwBothNucFilt(:,:,j);
    windowSize = 21; 
    expFactor = 1.2;
    kernel = ones(windowSize) / (windowSize^expFactor);
    bw = conv2(bw, kernel, 'same');
    bw2 = bw > 0.95;
    bw3 = imfill(bw2, 'holes');
    projNucMask(:,:,j) = (bwareaopen(bw3, 2000));
    projNucLabelMat(:,:,j) =  bwlabel( projNucMask(:,:,j));
    
    %% Find approx centroid of z project nuclei
    convHullProjNuc(:,:,j) = bwconvhull(projNucMask(:,:,j),'objects');
    hullLabelMat(:,:,j) =  bwlabel(convHullProjNuc(:,:,j));   
    hullProp = regionprops(bwconncomp(convHullProjNuc),'Centroid');
    projNucCent{j} = cat(1,hullProp.Centroid);   
    
    %% Find spots for all frames
    for i=1:slice  
        spotDetect = 3;
        spotFilter = 2;
        [redSpotMask(:,:,i,j), redSpotNL(:,:,i,j)] = Spots.spotDetect(reScaleRedStack(:,:,i,j), ...
            projNucLabelMat(:,:,j), spotElementSize, spotDetect, spotFilter, minSpotSize);
        redSpotMask2 = bwareaopen(redSpotMask(:,:,i,j),60);
        redSpotMask(:,:,i,j) = redSpotMask2;
        redSpotProject(:,:,j) = redSpotProject(:,:,j) + redSpotMask(:,:,i,j); 
        spotDetect = 1;
        spotFilter = 1;
        [greenSpotMask(:,:,i,j), greenSpotNL(:,:,i,j)] = Spots.spotDetect(reScaleGreenStack(:,:,i,j), ...
            projNucLabelMat(:,:,j), spotElementSize, spotDetect, spotFilter, minSpotSize);
        greenSpotProject(:,:,j) = greenSpotProject(:,:,j) + greenSpotMask(:,:,i,j); 
        bothSpotMask(:,:,:,i,j) = imfuse(redSpotMask(:,:,i,j),greenSpotMask(:,:,i,j),  'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    end
   CCred = bwconncomp(redSpotMask(:,:,:,j));
   cent = regionprops3(CCred,'Centroid');
   redCent{j} = cat(1, cent.Centroid);
   
end

HelperFunctions.imshow3D(projNucMask);
end