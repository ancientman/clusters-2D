function [nucBleachStruct, spotHotspotBleachStruct] = intensityCompare(im, metaDataDS, breakFrame, hotSpotS)
nNuc = size(hotSpotS.hotSpotsPerNucNL, 2);
deltaT  = metaDataDS.imagingInfo.DeltaT; %in seconds
% analysisTime = metaDataDS.analysisInfo.analysisTime; % in seconds
% analysisFrame = floor(analysisTime/deltaT); %in frames
% finalFrame = min(breakFrame, analysisFrame); % in frame
% finalTime = deltaT*finalFrame; % in seconds
% timeWindow = metaDataDS.analysisInfo.timeWindow; % in seconds
% frameWindow = round(timeWindow/deltaT); % in seconds
% timeSlide = metaDataDS.analysisInfo.timeSlide; % in seconds
% frameSlide = round(timeSlide/deltaT);
% startEndTimeCell = {};
% startEndFrameCell = {};
% 
% i=1;
% while (timeSlide*(i-1) + timeWindow)<=finalTime
%     startEndTimeCell{i} = [timeSlide*(i-1)+deltaT; timeSlide*(i-1) + (timeWindow+deltaT)];
%     i = i+1;
% end
% 
% i=1;
% while (frameSlide*(i-1) +frameWindow)<=finalFrame
%     startEndFrameCell{i} = [frameSlide*(i-1)+1; frameSlide*(i-1) + (frameWindow+1)];
%     i = i+1;
% end


startEndIterCell = cell(1,numel(hotSpotS.timeWindow));

for j = 1:nNuc
    for i=1:numel(hotSpotS.timeWindow)
        startEndIterCell{i}(1,:) = ceil(hotSpotS.timeWindow{i}(1)/deltaT);
        startEndIterCell{i}(2,:) = ceil(hotSpotS.timeWindow{i}(2)/deltaT);
        if startEndIterCell{i}(2,:)>breakFrame
            startEndIterCell{i}(2,:) = breakFrame;
        end
        
       meanPixelValues(im, hotSpotS.hotSpotIdxNL{j}{i}, ...
        startEndIterCell{i}(1,:), startEndIterCell{i}(2,:));
        
        
        if startEndIterCell{i}(2,:)>=breakFrame
            break;
        end
    end
    
end


%% Nuclear intensities w spots
% imNucAll = zeros(size(im, 2), size(im, 1), tFinal);
% imNucNoSpot = zeros(size(im, 2), size(im, 1), tFinal);
% imNucNoHotSpot = zeros(size(im, 2), size(im, 1), tFinal);
% imSpotInHotSpot = zeros(size(im, 2), size(im, 1), tFinal);
% imSpotOutHotSpot = zeros(size(im, 2), size(im, 1), tFinal);
% meanNucAll = cell(1, nNuc);
% meanNucNoSpot = cell(1, nNuc);
% meanNucNoHotSpot = cell(1, nNuc);
% meanSpotInHotSpot = cell(1, nNuc);
% meanSpotOutHotSpot = cell(1, nNuc);
% for t = 1:tFinal
%     for i = 1:nNuc
%         imNucTemp = im(:,:,t);
%         imNucTemp(~(regNucS{t}.labelMat==i)) = 0;
%         imHotSpotTemp = imNucTemp;
%         imHotSpotTemp((hotSpotS.hotspotNL==i)) = 0;
%         imSpotTemp = imNucTemp;
%         imSpotTemp((regSpotS{t}.labelMat==i)) = 0;        
%         meanNucAll{i}(t) = mean(imNucTemp(~(imNucTemp==0)),'all');   
%         meanNucNoHotSpot{i}(t) = mean(imHotSpotTemp(~(imHotSpotTemp==0)),'all');
%         meanNucNoSpot{i}(t) = mean(imSpotTemp(~(imSpotTemp==0)),'all');
%         imNucAll(:,:,t) = imNucAll(:,:,t) + imNucTemp;
%         imNucNoSpot(:,:,t) = imNucNoSpot(:,:,t) + imSpotTemp;
%         imNucNoHotSpot(:,:,t) = imNucNoHotSpot(:,:,t) + imHotSpotTemp;
%         imSpotInHotSpotTemp = imNucTemp;
%         if numel(regSpotS{t}.idxInHotspot{i})>0
%             imSpotInHotSpotTemp2(vertcat(regSpotS{t}.idxInHotspot{i}{1,:})) = 1;
%             imSpotInHotSpotTemp(imSpotInHotSpotTemp2~=1) = 0;
%             meanSpotInHotSpot{i}(t) = mean(imSpotInHotSpotTemp(~(imSpotInHotSpotTemp==0)),'all');   
%         end
%         imSpotOutHotSpotTemp = imNucTemp;
%         if numel(regSpotS{t}.idxOutHotspot{i})>0
%             imSpotOutHotSpotTemp3(vertcat(regSpotS{t}.idxOutHotspot{i}{1,:})) = 1;
%             imSpotOutHotSpotTemp(imSpotOutHotSpotTemp3~=1) = 0;
%             meanSpotOutHotSpot{i}(t) = mean(imSpotOutHotSpotTemp(~(imSpotOutHotSpotTemp==0)),'all');
%         end
%     end    
% end

end

function meanPixelValues(im, hotSpotPix, startFrame, endFrame)
imSpot = im(:,:,startFrame:endFrame);
mask = zeros(size(imSpot, 2), size(imSpot, 1));
mask1 = zeros(size(imSpot, 2), size(imSpot, 1));
meanRandRelIntensity = zeros(size(imSpot, 3), 1);
meanRandAbsIntensity = zeros(size(imSpot, 3), 1);
mask(hotSpotPix) = 1;
CC = bwconncomp(mask);
if CC.NumObjects ~=0    
    S = regionprops(CC, 'centroid');
    spotCentroid = round(cat(1, S.Centroid));
    binSize = 11;
    offset = 25;
    for ii = spotCentroid (1, 2) + offset - (binSize-1)/2: spotCentroid (1, 2) + offset + (binSize-1)/2
        for jj = spotCentroid (1, 1) - (binSize-1)/2: spotCentroid (1, 1) + (binSize-1)/2
            mask1(ii, jj) = 1;
        end
    end    
    
    meanSpotRelIntensity = cell(1, CC.NumObjects);
    meanSpotAbsIntensity = cell(1, CC.NumObjects);
    spotCorr = cell(1, CC.NumObjects);
    spotDetCorr = cell(1, CC.NumObjects);
    
    for k = 1:CC.NumObjects
        
        for i=startFrame:endFrame
            imSpotTemp = imSpot(:,:,i);    
            if k==1
                imSpotTemp1 = imSpotTemp;
                imSpotTemp1(~mask1) = 0;
                meanRandRelIntensity(i-startFrame+1) = mean(nonzeros(rescale( imSpotTemp1)), 'all');
                meanRandAbsIntensity(i-startFrame+1) = mean(nonzeros(imSpotTemp1), 'all');                
            end
            meanSpotRelIntensity{k} = vertcat(meanSpotRelIntensity{k}, mean(nonzeros(rescale(imSpotTemp(CC.PixelIdxList{k})))));
            meanSpotAbsIntensity{k} = vertcat(meanSpotAbsIntensity{k}, mean(nonzeros(imSpotTemp(CC.PixelIdxList{k}))));
        end
        spotCorr{k} = xcorr(meanSpotAbsIntensity{k});
        spotDetCorr{k} = xcorr(detrend(meanSpotAbsIntensity{k}));
        if k==1
            randCorr{k} = xcorr(meanRandAbsIntensity);
            randDetCorr{k} = xcorr(detrend(meanRandAbsIntensity));
        end
    end
end

end