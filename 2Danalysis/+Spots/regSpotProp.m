function [regSpotStruct, globalSpotPropStruct] = regSpotProp(nucPropStruct, regNucStruct, spotStruct, globalHotspot, registeredFrames, startFrame, endFrame)
nNuc = max(regNucStruct{1}.labelMat, [], 'all');
regSpotStruct = struct([]);
globalSpotPropStruct = struct([]);
tempLabelMat = zeros(size(regNucStruct{1}.labelMat));
for t = startFrame:endFrame
    regSpotStruct{t}.labelMat = zeros(size(regNucStruct{1}.labelMat)); %Nuclear labels per time frames of global hotspot regions
    regSpotStruct{t}.spotCentroid = cell(1,nNuc);
    regSpotStruct{t}.centroidInHotspot = cell(1,nNuc);
    regSpotStruct{t}.centroidOutHotspot = cell(1,nNuc);
    regSpotStruct{t}.areaInHotspot = cell(1,nNuc);
    regSpotStruct{t}.areaOutHotspot = cell(1,nNuc);
    regSpotStruct{t}.idxInHotspot = cell(1,nNuc);
    regSpotStruct{t}.idxOutHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityTotInHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityRelInHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityAbsInHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityTotOutHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityRelOutHotspot = cell(1,nNuc);    
    regSpotStruct{t}.intensityAbsOutHotspot = cell(1,nNuc);   
    regSpotStruct{t}.intensityDevInHotspot = cell(1,nNuc);
    regSpotStruct{t}.intensityDevOutHotspot = cell(1,nNuc);    
    regSpotStruct{t}.timeUsed = cell(1,nNuc);
    regSpotStruct{t}.inMolCount = cell(1,nNuc);
    regSpotStruct{t}.inSpotCount = cell(1,nNuc);
    regSpotStruct{t}.outMolCount = cell(1,nNuc);
    regSpotStruct{t}.outSpotCount = cell(1,nNuc);
end
centroidInHotspotNL = cell(1,nNuc); %accumulate for all time points spot centroids, only within global hotspots, labeled by nuc
areaInHotspotNL = cell(1,nNuc); %accumulate for all time points spot area, only within global hotspots, labeled by nuc
idxInHotspotNL = cell(1,nNuc); %accumulate for all time points spot indices, outside global hotspots, labeled by nuc
intensityTotInHotspotNL = cell(1,nNuc); % accumulate for all time points spot total intensity (sum of all pixel values in a spot), only within global hotspots, labeled by nuc
intensityRelInHotspotNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, only within global hotspots, labeled by nuc
nucIntensityInNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, only within global hotspots, labeled by nuc
intensityAbsInHotspotNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, only within global hotspots, labeled by nuc
intensityDevInHotspotNL = cell(1,nNuc); %accumulate for all time points spot intensity deviation, only within global hotspots, labeled by nuc
centroidOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot centroids, outside global hotspots, labeled by nuc
areaOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot area, outside global hotspots, labeled by nuc
idxOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot indices, outside global hotspots, labeled by nuc
intensityTotOutHotspotNL = cell(1,nNuc); % accumulate for all time points spot total intensity (sum of all pixel values in a spot), only outside global hotspots, labeled by nuc
intensityRelOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, outside global hotspots, labeled by nuc
nucIntensityOutNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, only within global hotspots, labeled by nuc
intensityAbsOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot mean intensity, outside global hotspots, labeled by nuc
intensityDevOutHotspotNL = cell(1,nNuc); %accumulate for all time points spot intensity deviation, outside global hotspots, labeled by nuc
timeUsedNL = cell(1,nNuc); %this is a logical array tagging when a label was considered for analysis
inMolCountNL = cell(1,nNuc); %total mols counted within hotspots it each time point per label
inSpotCountNL = cell(1,nNuc); %total spots counted within hotspots it each time point per label
outMolCountNL = cell(1,nNuc); %total mols counted outside hotspots it each time point per label
outSpotCountNL = cell(1,nNuc); %total spots counted outside hotspots it each time point per label
for t = startFrame:endFrame
    for k = 1:nNuc
        tempRegSpot = [];
        if ismember(t, registeredFrames)
            tempLabelMat(globalHotspot==1) = k;
            tempLabelMat(regNucStruct{t}.labelMat~=k) = 0; 
            regSpotStruct{t}.labelMat = regSpotStruct{t}.labelMat + tempLabelMat;    
%             nonzeroIdx = find(sum(spotStruct{t}.spotCentroid{k},2));
            if ~(isempty(find(sum(spotStruct{t}.spotCentroid{k},2), 1)))
%                 regSpotStruct{t}.spotCentroid{k}(:,1) = round(spotStruct{t}.spotCentroid{k}(:,1) - regNucStruct{t}.shift(k,1));
%                 regSpotStruct{t}.spotCentroid{k}(:,2) = round(spotStruct{t}.spotCentroid{k}(:,2) - regNucStruct{t}.shift(k,2));
                tempRegSpot(:,1) = round(spotStruct{t}.spotCentroid{k}(:,1) - regNucStruct{t}.shift(k,1));
                tempRegSpot(:,2) = round(spotStruct{t}.spotCentroid{k}(:,2) - regNucStruct{t}.shift(k,2));
                j = 1;
                for i = 1:size(tempRegSpot,1) 
                    if tempRegSpot(i,1)>0 && tempRegSpot(i,2)>0 && tempRegSpot(i,1)<size(regSpotStruct{t}.labelMat,2) && tempRegSpot(i,2)<size(regSpotStruct{t}.labelMat,2)
                        regSpotStruct{t}.spotCentroid{k}(j,1) = tempRegSpot(i,1);
                        regSpotStruct{t}.spotCentroid{k}(j,2) = tempRegSpot(i,2);
                        j = j+1;                    
                    end
                end
                if ~isempty(regSpotStruct{t}.spotCentroid{k})
                    regSpotCentIdx = sub2ind(size(regNucStruct{1}.labelMat), regSpotStruct{t}.spotCentroid{k}(:,1), regSpotStruct{t}.spotCentroid{k}(:,2));
                    centroidInHotspotIdx = intersect(regSpotCentIdx, find(tempLabelMat'), 'stable');
                    centroidOutHotspotIdx = setdiff(regSpotCentIdx, find(tempLabelMat'), 'stable');
                else
                    centroidInHotspotIdx = [];
                    centroidOutHotspotIdx = [];
                end
                if ~(isempty(centroidOutHotspotIdx))
                    [regSpotStruct{t}.centroidOutHotspot{k}(:,1), regSpotStruct{t}.centroidOutHotspot{k}(:,2)] =...
                        ind2sub(size(regNucStruct{1}.labelMat), centroidOutHotspotIdx);        
                    [~, outHotspotIdx, ~] = intersect(regSpotStruct{t}.spotCentroid{k}, regSpotStruct{t}.centroidOutHotspot{k}, 'rows');
                    regSpotStruct{t}.areaOutHotspot{k} = spotStruct{t}.spotArea{k}(outHotspotIdx);
                    regSpotStruct{t}.idxOutHotspot{k} = spotStruct{t}.spotIdx{k}(outHotspotIdx);
                    regSpotStruct{t}.intensityTotOutHotspot{k} = spotStruct{t}.spotAbsTotVal{k}(outHotspotIdx);
                    regSpotStruct{t}.intensityRelOutHotspot{k} = spotStruct{t}.spotRelMean{k}(outHotspotIdx);
                    regSpotStruct{t}.intensityAbsOutHotspot{k} = spotStruct{t}.spotAbsMean{k}(outHotspotIdx);
                    regSpotStruct{t}.intensityDevOutHotspot{k} = spotStruct{t}.spotDev{k}(outHotspotIdx);
                    regSpotStruct{t}.outSpotCount{k} = length(regSpotStruct{t}.areaOutHotspot{k});
                    regSpotStruct{t}.outMolCount{k} = sum(spotStruct{t}.spotMols{k}(outHotspotIdx));
                    
                else
                    regSpotStruct{t}.outSpotCount{k} = 0;
                end
                if ~(isempty(centroidInHotspotIdx))
                    [regSpotStruct{t}.centroidInHotspot{k}(:,1), regSpotStruct{t}.centroidInHotspot{k}(:,2)] =...
                        ind2sub(size(regNucStruct{1}.labelMat), centroidInHotspotIdx);        
                    [~, inHotspotIdx, ~] = intersect(regSpotStruct{t}.spotCentroid{k}, regSpotStruct{t}.centroidInHotspot{k}, 'rows');
                    regSpotStruct{t}.areaInHotspot{k} = spotStruct{t}.spotArea{k}(inHotspotIdx);
                    regSpotStruct{t}.idxInHotspot{k} = spotStruct{t}.spotIdx{k}(inHotspotIdx);
                    regSpotStruct{t}.intensityTotInHotspot{k} = spotStruct{t}.spotAbsTotVal{k}(inHotspotIdx);   
                    regSpotStruct{t}.intensityRelInHotspot{k} = spotStruct{t}.spotRelMean{k}(inHotspotIdx);   
                    regSpotStruct{t}.intensityAbsInHotspot{k} = spotStruct{t}.spotAbsMean{k}(inHotspotIdx);   
                    regSpotStruct{t}.intensityDevInHotspot{k} = spotStruct{t}.spotDev{k}(inHotspotIdx);  
                    regSpotStruct{t}.timeUsed{k} = 1;                    
                    regSpotStruct{t}.inSpotCount{k} = length(regSpotStruct{t}.areaInHotspot{k});
                    regSpotStruct{t}.inMolCount{k} = sum(spotStruct{t}.spotMols{k}(inHotspotIdx));
                else
                    regSpotStruct{t}.timeUsed{k} = 0;
                    regSpotStruct{t}.inMolCount{k} = 0;
                    regSpotStruct{t}.inSpotCount{k} = 0;
                end
            else
                regSpotStruct{t}.timeUsed{k} = 0;
                regSpotStruct{t}.inMolCount{k} = 0;
                regSpotStruct{t}.inSpotCount{k} = 0;
            end
            tempLabelMat(:) = 0;
        else
            regSpotStruct{t}.timeUsed{k} = 0;
            regSpotStruct{t}.inMolCount{k} = 0;
            regSpotStruct{t}.inSpotCount{k} = 0;
        end
        
        nucIntensityInNL{1, k} = vertcat(nucIntensityInNL{1, k}, nucPropStruct{t}.meanAbsIntensity(k).*ones(length(regSpotStruct{t}.intensityRelInHotspot{k}), 1));
        nucIntensityOutNL{1, k} = vertcat(nucIntensityOutNL{1, k}, nucPropStruct{t}.meanAbsIntensity(k).*ones(length(regSpotStruct{t}.intensityRelOutHotspot{k}), 1));
        intensityTotInHotspotNL{1,k} = vertcat(intensityTotInHotspotNL{1,k}, (regSpotStruct{t}.intensityTotInHotspot{k}));
        intensityRelInHotspotNL{1,k} = vertcat(intensityRelInHotspotNL{1,k}, (regSpotStruct{t}.intensityRelInHotspot{k}));
        intensityAbsInHotspotNL{1,k} = vertcat(intensityAbsInHotspotNL{1,k}, (regSpotStruct{t}.intensityAbsInHotspot{k}));
        intensityDevInHotspotNL{1,k} = vertcat(intensityDevInHotspotNL{1,k}, (regSpotStruct{t}.intensityDevInHotspot{k}));
        areaInHotspotNL{1,k} = vertcat(areaInHotspotNL{1,k}, (regSpotStruct{t}.areaInHotspot{k}));
        centroidInHotspotNL{1,k} = vertcat(centroidInHotspotNL{1,k}, regSpotStruct{t}.centroidInHotspot{k});   
        intensityTotOutHotspotNL{1,k} = vertcat(intensityTotOutHotspotNL{1,k}, (regSpotStruct{t}.intensityTotOutHotspot{k}));
        intensityRelOutHotspotNL{1,k} = vertcat(intensityRelOutHotspotNL{1,k}, (regSpotStruct{t}.intensityRelOutHotspot{k}));
        intensityAbsOutHotspotNL{1,k} = vertcat(intensityAbsOutHotspotNL{1,k}, (regSpotStruct{t}.intensityAbsOutHotspot{k}));
        intensityDevOutHotspotNL{1,k} = vertcat(intensityDevOutHotspotNL{1,k}, (regSpotStruct{t}.intensityDevOutHotspot{k}));
        areaOutHotspotNL{1,k} = vertcat(areaOutHotspotNL{1,k}, (regSpotStruct{t}.areaOutHotspot{k}));
        centroidOutHotspotNL{1,k} = vertcat(centroidOutHotspotNL{1,k}, regSpotStruct{t}.centroidOutHotspot{k});        
        timeUsedNL{1,k} = vertcat(timeUsedNL{1,k}, regSpotStruct{t}.timeUsed{k});
        inMolCountNL{1,k} = vertcat(inMolCountNL{1,k}, regSpotStruct{t}.inMolCount{k});
        inSpotCountNL{1,k} = vertcat(inSpotCountNL{1,k}, regSpotStruct{t}.inSpotCount{k});
        outMolCountNL{1,k} = vertcat(outMolCountNL{1,k}, regSpotStruct{t}.outMolCount{k});
        outSpotCountNL{1,k} = vertcat(outSpotCountNL{1,k}, regSpotStruct{t}.outSpotCount{k});
    end    
end
for k=1:nNuc
    globalSpotPropStruct{k}.nucIntensityInNL = nucIntensityInNL{k};
    globalSpotPropStruct{k}.nucIntensityOutNL = nucIntensityOutNL{k};
    globalSpotPropStruct{k}.centroidInHotspotNL = centroidInHotspotNL{k};
    globalSpotPropStruct{k}.areaInHotspotNL = areaInHotspotNL{k};
    globalSpotPropStruct{k}.intensityTotInHotspotNL = intensityTotInHotspotNL{k};    
    globalSpotPropStruct{k}.intensityRelInHotspotNL = intensityRelInHotspotNL{k};    
    globalSpotPropStruct{k}.intensityAbsInHotspotNL = intensityAbsInHotspotNL{k};    
    globalSpotPropStruct{k}.intensityDevInHotspotNL = intensityDevInHotspotNL{k}; 
    globalSpotPropStruct{k}.centroidOutHotspotNL = centroidOutHotspotNL{k};
    globalSpotPropStruct{k}.areaOutHotspotNL = areaOutHotspotNL{k};
    globalSpotPropStruct{k}.intensityTotOutHotspotNL = intensityTotOutHotspotNL{k};  
    globalSpotPropStruct{k}.intensityRelOutHotspotNL = intensityRelOutHotspotNL{k};  
    globalSpotPropStruct{k}.intensityAbsOutHotspotNL = intensityAbsOutHotspotNL{k};      
    globalSpotPropStruct{k}.intensityDevOutHotspotNL = intensityDevOutHotspotNL{k};
    globalSpotPropStruct{k}.timeUsedNL = timeUsedNL{k};
    globalSpotPropStruct{k}.inMolCountNL = inMolCountNL{k};
    globalSpotPropStruct{k}.inSpotCountNL = inSpotCountNL{k};
    globalSpotPropStruct{k}.outMolCountNL = outMolCountNL{k};
    globalSpotPropStruct{k}.outSpotCountNL = outSpotCountNL{k};
end
end