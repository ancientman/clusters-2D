function [registerNucPropStruct, regFrames] = registeredNucProp(nucPropStruct, im, maxDistShift, totalFrames)
nNuc = max(nucPropStruct{1}.labelMat, [], 'all');
tIter = linspace(1, totalFrames, totalFrames);
templateMat = cell(1,nNuc);
for k=1:nNuc
    templateMat{k} = zeros(size(nucPropStruct{1}.labelMat));
end
bw = zeros(size(nucPropStruct{1}.labelMat));
registerStruct = struct([]);
for t = 1:totalFrames
    registerStruct{t}.labelMat = zeros(size(nucPropStruct{1}.labelMat));
    nucPropStruct{t}.centroid = (nucPropStruct{t}.centroid);
end
for k = 1:nNuc 
    for t = 1:totalFrames
        bw(:) = 0;
        templateMat{k}(:) = 0;
        pixelIndexList = label2idx(nucPropStruct{t}.labelMat);
        templateMat{k}(pixelIndexList{k}) = 1;
        boundaries = bwboundaries(templateMat{k}, 'noholes');
        if t==1
            registerStruct{1}.shift{:,k}(:,1) = 0;
            registerStruct{1}.shift{:,k}(:,2) = 0;
            registerStruct{1}.centroid{:,k}(:,1) = nucPropStruct{t}.centroid(k,1)-...
                registerStruct{1}.shift{:,k}(:,1);
            registerStruct{1}.centroid{:,k}(:,2) = nucPropStruct{t}.centroid(k,2)-...
                registerStruct{1}.shift{:,k}(:,2);
            bw = poly2mask(boundaries{1}(:,2)-registerStruct{1}.shift{:,k}(:,1),...
                boundaries{1}(:,1)-registerStruct{1}.shift{:,k}(:,2), ...
                size(nucPropStruct{1}.labelMat,2), size(nucPropStruct{1}.labelMat,1));
            registerStruct{1}.labelMat = registerStruct{1}.labelMat + k.*bw;
        else
            if (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
                    nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], 'euclidean')...
                    >maxDistShift)
                %this centroid doesnot satisfy constraint with the previous,
                if tIter(t-1)==0 
                    %the previous time point is skipped
                    tTemp = t-1;
                    while tIter(tTemp)==0 %find the last non skipped time point
                        tTemp = tTemp-1;
                    end
                    %tTemp = tTemp-1; %last non-skipped time point
                    if (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
                        nucPropStruct{tTemp}.centroid(k,1), nucPropStruct{tTemp}.centroid(k,2)], 'euclidean')...
                        >=maxDistShift)
                    %if current still does not satisfy the constraint with the
                    %last non-skipped point: skip current
                        fprintf('\n skipping t = %d\n', t);
                        tIter(t) = 0;
                    elseif(pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
                        nucPropStruct{tTemp}.centroid(k,1), nucPropStruct{tTemp}.centroid(k,2)], 'euclidean')...
                        <maxDistShift)   
                        %if satisfies constraint with the last non-skipped time point: plot current
                        registerStruct{t}.shift{:,k}(:,1) = (nucPropStruct{t}.centroid(k,1)-nucPropStruct{tTemp}.centroid(k,1));
                        registerStruct{t}.shift{:,k}(:,2) = (nucPropStruct{t}.centroid(k,2)-nucPropStruct{tTemp}.centroid(k,2));
                        registerStruct{t}.centroid{:,k}(:,1) = nucPropStruct{t}.centroid(k,1)-...
                            registerStruct{t}.shift{:,k}(:,1);
                        registerStruct{t}.centroid{:,k}(:,2) = nucPropStruct{t}.centroid(k,2)-...
                            registerStruct{t}.shift{:,k}(:,2);
                          bw = poly2mask(boundaries{1}(:,2)-registerStruct{t}.shift{:,k}(:,1), ...
                            boundaries{1}(:,1)-registerStruct{t}.shift{:,k}(:,2), ...
                            size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1));
%                             BW = poly2mask(boundaries{k}(:,2), ...
%                                 boundaries{k}(:,1), ...
%                                 size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1));
                        registerStruct{t}.labelMat = registerStruct{t}.labelMat + k.*bw;
%                         fprintf('dist - %d and time is %d\n', ...
%                             round(pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2); ...
%                             nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], ...
%                             'euclidean'),2), tIter(t));
                    end
                elseif tIter(t-1)~=0
                    %this centroid doesnot satisfy constraint with the previous,
                    %however the previous time point is non-skipped: skip
                    %current
                    fprintf('\n skipping t = %d\n', t); 
                    tIter(t) = 0;
                end
            elseif (pdist([nucPropStruct{t}.centroid(k,1), nucPropStruct{t}.centroid(k,2);...
                    nucPropStruct{t-1}.centroid(k,1), nucPropStruct{t-1}.centroid(k,2)], 'euclidean')...
                    <maxDistShift) 
                if tIter(t-1)==0 
                    %if current satisfies the constraint with previous
                    %however the last time point is skipped: skip current
                    fprintf('\n skipping t = %d\n', t); %skip current
                    tIter(t) = 0;
                elseif tIter(t-1)~=0 
                    %if current satisfies the constraint with the previous
                    %however the last time point is non-skipped: plot current
                    registerStruct{t}.shift{:,k}(:,1) = (nucPropStruct{t}.centroid(k,1)-nucPropStruct{t-1}.centroid(k,1));
                    registerStruct{t}.shift{:,k}(:,2) = (nucPropStruct{t}.centroid(k,2)-nucPropStruct{t-1}.centroid(k,2));
                    registerStruct{t}.centroid{:,k}(:,1) = nucPropStruct{t}.centroid(k,1)-registerStruct{t}.shift{:,k}(:,1);
                    registerStruct{t}.centroid{:,k}(:,2) = nucPropStruct{t}.centroid(k,2)-registerStruct{t}.shift{:,k}(:,2);
                    bw = poly2mask(boundaries{1}(:,2)-registerStruct{t}.shift{:,k}(:,1), ...
                        boundaries{1}(:,1)-registerStruct{t}.shift{:,k}(:,2), ...
                        size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1));                    
%                     BW = poly2mask(boundaries{k}(:,2), ...
%                         boundaries{k}(:,1), ...
%                         size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1));                    
                    registerStruct{t}.labelMat = registerStruct{t}.labelMat + k.*bw;
%                     fprintf('dist - %d and time is %d\n', round(pdist([nucPropStruct{t}.centroid(k,1), ...
%                         nucPropStruct{t}.centroid(k,2); nucPropStruct{t-1}.centroid(k,1), ...
%                         nucPropStruct{t-1}.centroid(k,2)], 'euclidean'),2), tIter(t));
                end
            end
        end
    end
end
for t= 1:totalFrames
    tempCentMat = cell2mat(registerStruct{t}.centroid');
    tempShiftMat = cell2mat(registerStruct{t}.shift');
    registerStruct{t}.centroid = round(tempCentMat);
    registerStruct{t}.shift = round(tempShiftMat);
end
% totalRegisteredFrames = nnz(tIter);
labelMatRegistered = zeros(size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1), totalFrames);
labelMatUnregistered = zeros(size(nucPropStruct{1}.labelMat, 2), size(nucPropStruct{1}.labelMat, 1), totalFrames);
for t=1:totalFrames
    labelMatUnregistered(:,:,t) = nucPropStruct{t}.labelMat;
    labelMatRegistered(:,:,t) = registerStruct{t}.labelMat;
end

labelMat = cell(1, totalFrames);
% CC = cell(1, nNuc);
f1 = @(x) getfield(x, 'labelMat');
labelMat = cellfun(f1, registerStruct, 'un', 0);
for t= 1:totalFrames
    for i = 1:nNuc        
        CC = bwconncomp(labelMat{t}==i);
        s = regionprops(CC, im(:,:,t), 'PixelValues', 'PixelList', 'PixelIdxList');
        registerStruct{t}.pixValue{i} = s.PixelValues;
        registerStruct{t}.pixList{i} = s.PixelList; 
        registerStruct{t}.pixIdxList{i} = s.PixelIdxList; 
    end
end
registerNucPropStruct = registerStruct;
regFrames = tIter;
end
