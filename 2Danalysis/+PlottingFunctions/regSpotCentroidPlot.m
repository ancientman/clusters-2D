function regSpotProp(regNucPropStruct, spotPropStruct, metaDataDS, bwSpot, registeredFrames, breakingTime)
% nNuc = max(nucPropStruct{1}.labelMat, [], 'all');
% tTime = breakingTime;
% xScale = metaDataDS.imagingInfo.XpixelSize/1000; %changed to microns;
% yScale = metaDataDS.imagingInfo.YpixelSize/1000; %changed to microns;
% timeScale = metaDataDS.imagingInfo.timeResolution; %time resolution in seconds;
% nucleusFeatureSize = metaDataDS.imagingInfo.nucleusFeatureSize;
% centroidShiftConstraint = ceil(nucleusFeatureSize/(5*xScale)); %permissible shift between centroids in consecutive time points
% tIter = linspace(1, tTime, tTime);
% spotCount = zeros(1,breakingTime);
% templateMat = zeros(size(nucPropStruct{1}.labelMat));
% idx = cell(1,nNuc);
% boundaries = cell(1,nNuc);
% thisBoundary = cell(1,nNuc);
% for k = 1:nNuc 
%     figure;
%     for t = 1:tTime        
%         pixelIndexList = label2idx(nucPropStruct{t}.labelMat);
%         idx{k} = pixelIndexList{k};
%         templateMat(idx{k}) = 1;
%         boundaries{k} = bwboundaries(templateMat, 'noholes');
%         thisBoundary{k} = boundaries{k};     
%         if t==1
%             xShift = 0;
%             yShift = 0;
% %             plot(spotPropStruct{t}.spotCentroid{:,k}(:,1), spotPropStruct{t}.spotCentroid{:,k}(:,2), '.b');
%             hold on;
%             plot(thisBoundary{k}{1}(:,2), thisBoundary{k}{1}(:,1), 'r', 'LineWidth', 2);
%             hold on;
%             plot(nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2), 'r*'); 
%             hold on;
% %             templateMat(:) = 0;
%         else
%             if (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
%                     nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], 'euclidean')...
%                     >centroidShiftConstraint)
%                 %this centroid doesnot satisfy constraint with the previous,
%                 if tIter(t-1)==0 
%                     %the previous time point is skipped
%                     tTemp = t-1;
%                     while tIter(tTemp)==0 %find the last non skipped time point
%                         tTemp = tTemp-1;
%                     end
%                     %tTemp = tTemp-1; %last non-skipped time point
%                     if (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
%                         nucPropStruct{tTemp}.centroid(k,1), nucPropStruct{tTemp}.centroid(k,2)], 'euclidean')...
%                         >=centroidShiftConstraint)
%                     %if current still does not satisfy the constraint with the
%                     %last non-skipped point: skip current
%                         fprintf('\n skipping t = %d\n', t);
%                         tIter(t) = 0;
%                     elseif(pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
%                         nucPropStruct{tTemp}.centroid(k,1), nucPropStruct{tTemp}.centroid(k,2)], 'euclidean')...
%                         <centroidShiftConstraint)   
%                         %if satisfies constraint with the last non-skipped time point: plot current
%                         xShift = nucPropStruct{t}.centroid(k,1) - nucPropStruct{tTemp}.centroid(k,1);
%                         yShift = nucPropStruct{t}.centroid(k,2) - nucPropStruct{tTemp}.centroid(k,2);
%                         plot(spotPropStruct{t}.spotCentroid{:,k}(:,1)-xShift, ...
%                             spotPropStruct{t}.spotCentroid{:,k}(:,2)-yShift, '.b');
%                         hold on;
%                         plot(nucPropStruct{t}.centroid(k,1)-xShift, nucPropStruct{t}.centroid(k,2)-yShift, '*black');
%                         hold on;
%                         plot(thisBoundary{k}{1}(:,2)-xShift, thisBoundary{k}{1}(:,1)-yShift, 'r', 'LineWidth', 2);                        
%                         hold on;                        
%                         fprintf('\nhere\n');
%                     end
%                 elseif tIter(t-1)~=0
%                     %this centroid doesnot satisfy constraint with the previous,
%                     %however the previous time point is non-skipped: skip
%                     %current
%                     fprintf('\n skipping t = %d\n', t); 
%                     tIter(t) = 0;
%                 end
%             elseif (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
%                     nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], 'euclidean')...
%                     <centroidShiftConstraint) 
%                 if tIter(t-1)==0 
%                     %if current satisfies the constraint with previous
%                     %however the last time point is skipped: skip current
%                     fprintf('\n skipping t = %d\n', t); %skip current
%                     tIter(t) = 0;
%                 elseif tIter(t-1)~=0 
%                     %if current satisfies the constraint with the previous
%                     %however the last time point is non-skipped: plot current
%                     xShift = nucPropStruct{t}.centroid(k,1) - nucPropStruct{t-1}.centroid(k,1);
%                     yShift = nucPropStruct{t}.centroid(k,2) - nucPropStruct{t-1}.centroid(k,2);
%                     plot(spotPropStruct{t}.spotCentroid{:,k}(:,1)-xShift, ...
%                         spotPropStruct{t}.spotCentroid{:,k}(:,2)-yShift, '.b');
%                     hold on;
%                     plot(nucPropStruct{t}.centroid(k,1)-xShift, nucPropStruct{t}.centroid(k,2)-yShift, '*black');
%                     hold on;     
%                     plot(thisBoundary{k}{1}(:,2)-xShift, thisBoundary{k}{1}(:,1)-yShift, 'r', 'LineWidth', 2);
%                     hold on;
%                     fprintf('dist - %d and time is %d\n', pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2); nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], 'euclidean'), tIter(t));
%                 end
%             end
%         end
%     end
%     xlim([0, size(nucPropStruct{1}.labelMat(:,2), 1)]);
%     ylim([0, size(nucPropStruct{1}.labelMat(:,1), 1)]);
%     xBound = [-50 0];
%     yBound = [0 25];
%     figure;
%     for t=1:tTime
%         spotDist = zeros(size(size(spotPropStruct{t}.spotCentroid{:,k}, 1), 2));
%         for i=1:size(spotPropStruct{t}.spotCentroid{:,k}, 1)
%             spotDist(i,1) = spotPropStruct{t}.spotCentroid{:,k}(i,1)-nucPropStruct{t}.centroid(k,1);
%             spotDist(i,2) = spotPropStruct{t}.spotCentroid{:,k}(i,2)-nucPropStruct{t}.centroid(k,2);
%             if(spotPropStruct{t}.spotCentroid{:,k}(i,1)-nucPropStruct{t}.centroid(k,1)>xBound(1,1) && ...
%                     spotPropStruct{t}.spotCentroid{:,k}(i,1)-nucPropStruct{t}.centroid(k,1)<xBound(1,2) && ...
%                     spotPropStruct{t}.spotCentroid{:,k}(i,2)-nucPropStruct{t}.centroid(k,2)>yBound(1,1) && ...
%                     spotPropStruct{t}.spotCentroid{:,k}(i,2)-nucPropStruct{t}.centroid(k,2)<yBound(1,2))            
%                 spotCount(t) = spotCount(t) + 1;
%                 plot(spotPropStruct{t}.spotCentroid{:,k}(i,1), spotPropStruct{t}.spotCentroid{:,k}(i,2), 'black.');
%                 hold on;
%     %             fprintf('here at%d\n',t);
%             else
%                 plot(spotPropStruct{t}.spotCentroid{:,k}(i,1), spotPropStruct{t}.spotCentroid{:,k}(i,2), 'r.');
%     %             fprintf('not here at%d\n',t);
%             end
%             line(linspace(nucPropStruct{t}.centroid(k,1)+xBound(1,1),(nucPropStruct{t}.centroid(k,1)+xBound(1,2)), 50),...
%                 linspace(nucPropStruct{t}.centroid(k,2)+yBound(1,1),(nucPropStruct{t}.centroid(k,2)+yBound(1,1)), 50));
%             hold on;
%             line(linspace(nucPropStruct{t}.centroid(k,1)+xBound(1,1),(nucPropStruct{t}.centroid(k,1)+xBound(1,1)), 50),...
%                 linspace(nucPropStruct{t}.centroid(k,2)+yBound(1,1),(nucPropStruct{t}.centroid(k,2)+yBound(1,2)), 50));
%             hold on;
%             line(linspace(nucPropStruct{t}.centroid(k,1)+xBound(1,2),(nucPropStruct{t}.centroid(k,1)+xBound(1,2)), 50),...
%                 linspace(nucPropStruct{t}.centroid(k,2)+yBound(1,1),(nucPropStruct{t}.centroid(k,2)+yBound(1,2)), 50));
%             hold on;
%             line(linspace(nucPropStruct{t}.centroid(k,1)+xBound(1,1),(nucPropStruct{t}.centroid(k,1)+xBound(1,2)),50),...
%                 linspace(nucPropStruct{t}.centroid(k,2)+yBound(1,2),(nucPropStruct{t}.centroid(k,2)+yBound(1,2)), 50));
%             hold on;
%         end
%     end
% end
% xlim([0, size(nucPropStruct{1}.labelMat(:,2), 1)]);
% ylim([0, size(nucPropStruct{1}.labelMat(:,1), 1)]);
% figure;
% plot(0.13.*linspace(1,tTime,tTime), spotCount);
end