function spotCentroidPlot2(DSfolder)
regSpotStruct = load(append(DSfolder, filesep, 'spotDS.mat'));
gSpotDS = load(append(DSfolder, filesep, 'globalSpotDS.mat'));
figure;
imshow(bwperim(regSpotStruct.globalRegSpotStruct{1}.labelMat),[]);
hold on;
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
totalFrames = length(regSpotStruct.globalRegSpotStruct);
for k=1:nNuc
    plot(gSpotDS.globalSpotHotSpotStruct{k}.centroidInHotspotNL(:,1), gSpotDS.globalSpotHotSpotStruct{k}.centroidInHotspotNL(:,2), ...
        '.g');%, 'color', [1/(k) 1/(2*k) 1/(1.5*k)]);
    hold on;
    plot(gSpotDS.globalSpotHotSpotStruct{k}.centroidOutHotspotNL(:,1), gSpotDS.globalSpotHotSpotStruct{k}.centroidOutHotspotNL(:,2), ...
        '.r');%, 'color', [1 0 0]);
    hold on;
end
xlim([0 size(regSpotStruct.globalRegSpotStruct{1}.labelMat, 2)]);
ylim([0 size(regSpotStruct.globalRegSpotStruct{1}.labelMat, 2)]);
axis equal
end