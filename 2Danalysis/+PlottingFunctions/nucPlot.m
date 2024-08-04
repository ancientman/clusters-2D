function nucPlot(nucPropStruct, metaDataDS, breakingTime)
XpixelSize = metaDataDS.imagingInfo.XpixelSize/1000; %changed to microns
YpixelSize = metaDataDS.imagingInfo.YpixelSize/1000; %changed to microns
nucleusFeatureSize = metaDataDS.imagingInfo.nucleusFeatureSize; %in microns
imDim = [metaDataDS.imagingInfo.Y,metaDataDS.imagingInfo.X];
nNuc = max(nucPropStruct{1}.labelMat(), [], 'all');
templateMat = zeros(imDim(1,2),imDim(1,1));
idx = cell(1,nNuc);
boundaries = cell(1,nNuc);
thisBoundary = cell(1,nNuc);
for i = 1:nNuc
    figure;
    for t = 1:breakingTime
        pixelIndexList = label2idx(nucPropStruct{t}.labelMat);
        idx{i} = pixelIndexList{i};
        templateMat(idx{i}) = 1;
        boundaries{i} = bwboundaries(templateMat, 'noholes');
        thisBoundary{i} = boundaries{i};  
        %Turn on to check where the boundaries at
%         if max(thisBoundary{i}{1}(:,1),[],'all')>floor(250)
%             fprintf('\nproblem time is %d\n',t);
%             break;
%         elseif max(thisBoundary{i}{1}(:,2),[],'all')>floor(250)
%             fprintf('\nproblem time is %d\n',t);
%             break;
%         end
        plot(thisBoundary{i}{1}(:,2), thisBoundary{i}{1}(:,1), 'b', 'LineWidth', 2);
        hold on;
        plot(nucPropStruct{t}.centroid(i,1), nucPropStruct{t}.centroid(i,2), 'r*'); 
        hold on;
        templateMat(:) = 0;
    end    
    xlim([1 imDim(1,2)]); ylim([1 imDim(1,1)]);
end
end