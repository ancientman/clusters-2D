function [hotSpotPropNL] = timeBinHotSpot(metaDataDS, registerNucPropStruct, imUse, nucPropStruct, spotPropStruct, registeredFrames)
analysisTime = metaDataDS.analysisInfo.analysisTime; % in seconds
timeWindow = metaDataDS.analysisInfo.timeWindow; % in seconds
timeSlide = metaDataDS.analysisInfo.timeSlide; % in seconds
startEndTimeCell = {};
DeltaT  = metaDataDS.imagingInfo.DeltaT;
nucThres = metaDataDS.imagingInfo.nucThres;
elementSize = metaDataDS.imagingInfo.elementSize;
nucleusDetectMode = metaDataDS.imagingInfo.nucleusDetectMode;
clearBorder = metaDataDS.imagingInfo.clearBorder;
XpixelSize = metaDataDS.imagingInfo.XpixelSize; % in microns
spotElementSize = metaDataDS.analysisInfo.spotElementSize;
minNucSize = metaDataDS.analysisInfo.minNucSize;
minSpotSize = metaDataDS.analysisInfo.minSpotSize;
breakingTimeFrame = metaDataDS.analysisInfo.breakingTimeFrame;
i = 1;
while (timeSlide*(i-1) + timeWindow)<=analysisTime
    startEndTimeCell{i} = [timeSlide*(i-1)+DeltaT; timeSlide*(i-1) + (timeWindow+DeltaT)];
    i = i+1;
end
startEndIterCell = cell(1,numel(startEndTimeCell));
imNormAv = cell(1,numel(startEndTimeCell));
hotSpotPropTL = cell(1,numel(startEndTimeCell));
allSpotProp = cell(1,numel(startEndTimeCell));
hotSpotStruct = cell(1,numel(startEndTimeCell));
spotHotSpotStruct = cell(1,numel(startEndTimeCell));
hotSpotBW = cell(1, numel(startEndTimeCell));

for i=1:numel(startEndTimeCell)
    startEndIterCell{i}(1,:) = ceil(startEndTimeCell{i}(1)/DeltaT);
    startEndIterCell{i}(2,:) = ceil(startEndTimeCell{i}(2)/DeltaT);
    if startEndIterCell{i}(2,:)>breakingTimeFrame
        startEndIterCell{i}(2,:) = breakingTimeFrame;
    end
    imNormAv{i} = Nucleus.imNormAv(registerNucPropStruct, imUse, ...
        startEndIterCell{i}(1,:), startEndIterCell{i}(2,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     imNormAv{i} = imgaussfilt(imNormAv{i}, [5,5]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hotSpotStruct{i} = Spots.hotSpotFinder(imNormAv{i}, nucPropStruct, ...
        nucThres, elementSize, spotElementSize, ...
        nucleusDetectMode, minNucSize, minSpotSize, clearBorder);
    
    hotSpotPropTL{i} = Spots.hotSpotProp(hotSpotStruct{i});
    
    [~, spotHotSpotStruct{i}] = Spots.regSpotProp(nucPropStruct, registerNucPropStruct, ...
        spotPropStruct, hotSpotStruct{i}.hotspotBW, ...
        registeredFrames, startEndIterCell{i}(1,:), startEndIterCell{i}(2,:));
    
    allSpotProp{i} = Spots.allNucSpotProp(hotSpotStruct{i}, spotHotSpotStruct{i}, XpixelSize, timeWindow);
    
    hotSpotBW{i} = hotSpotStruct{i}.hotspotBW;  
    if startEndIterCell{i}(2,:)>=breakingTimeFrame
        break;
    end
end

hotSpotAreaNL = cell(1, numel(hotSpotStruct{1}.hotspotAreaUniq));
meanHotSpotAreaNL = cell(1, numel(hotSpotPropTL{1}.meanHotSpotArea));
hotSpotsPerNucNL = cell(1, numel(hotSpotPropTL{1}.hotSpotsPerNuc));
hotSpotCentroidUniqNL = cell(1, numel(hotSpotPropTL{1}.hotSpotsPerNuc));
hotSpotIdxUniqNL = cell(1, numel(hotSpotPropTL{1}.hotSpotsPerNuc));

for i=1:numel(hotSpotStruct{1}.hotspotAreaUniq) % total nuclei
    hotSpotCentUniqTemp = {};
    hotSpotIdxUniqTemp = {};
    for j=1:numel(hotSpotStruct) % total time points
        if ~isempty(hotSpotPropTL{j}) && ~isempty(hotSpotStruct{j})
            hotSpotAreaNL{i} = vertcat(hotSpotAreaNL{i}, hotSpotStruct{j}.hotspotAreaUniq{i});
            meanHotSpotAreaNL{i} = vertcat(meanHotSpotAreaNL{i}, hotSpotPropTL{j}.meanHotSpotArea{i});
            hotSpotsPerNucNL{i} = vertcat(hotSpotsPerNucNL{i}, hotSpotPropTL{j}.hotSpotsPerNuc{i});
            hotSpotCentUniqTemp{end+1} = hotSpotStruct{j}.hotspotCentroidUniq{i}; 
            hotSpotIdxUniqTemp{end+1} = hotSpotStruct{j}.hotspotIdxListUniq{i};
        end
    end
    hotSpotCentroidUniqNL{i} = hotSpotCentUniqTemp';
    hotSpotIdxUniqNL{i} = hotSpotIdxUniqTemp';
end

hotSpotRadNL = cellfun(@(C) sqrt(C*(XpixelSize)^2/pi), hotSpotAreaNL, 'UniformOutput', false); % Area in um square
hotSpotPropNL.hotSpotBW = hotSpotBW; % in pixels
hotSpotPropNL.hotSpotAreaNL = hotSpotAreaNL; % in pixels
hotSpotPropNL.hotSpotRadNL = hotSpotRadNL; % in pixels
hotSpotPropNL.meanHotSpotAreaNL = meanHotSpotAreaNL; % in pixels
hotSpotPropNL.hotSpotsPerNucNL = hotSpotsPerNucNL; % in pixels
hotSpotPropNL.hotSpotCentroidNL = hotSpotCentroidUniqNL; % in pixels
hotSpotPropNL.hotSpotIdxNL = hotSpotIdxUniqNL; % in pixels
hotSpotPropNL.timeWindow = startEndTimeCell;
%%%%%%%%%%%%%Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = num2cell(1: numel(hotSpotRadNL)); % create a label cell
% x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], hotSpotRadNL, y, 'UniformOutput', 0);
% X = vertcat(x{:});
% boxplot(X(:,1), X(:, 2));
end