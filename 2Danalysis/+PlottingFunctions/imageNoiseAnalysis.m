function imageNoiseAnalysis()

% filePath = 'E:\Analysis\1c1P\singlePlane\bead\bd4um_3.tif';
% xAxis = 160:200;
% yAxis = 160:220;

filePath = 'E:\Dropbox (Princeton)\Data_from_Rahul\bcd_integration_time\f3_bcd2x_gfp_532ms.tif';
% xAxis = 60:120;
% yAxis = 60:120;

xAxis = 100:120;
yAxis = 100:120;

% filePath = 'E:\Analysis\1c1P\singlePlane\gfp\20210412_nls_gfp_5.tif';
% xAxis = 110:160;
% yAxis = 170:220;

% filePath = 'E:\Analysis\1c1P\singlePlane\2xa\bcd2xa_ant_nc14_em1_530ms_1.tif';


stack = HelperFunctions.loadTiff(filePath); 
stack1 = rescale(stack(:,:,1));
stack1 = medfilt2(stack1, [5, 5]);

%------------------------------------------------------------
%   Use fixed ROI
stack1Mod = stack1(xAxis, yAxis); 
stackMod = stack(xAxis, yAxis, :); 
[X, Y] = meshgrid(xAxis, yAxis);
roi(X,Y) = 1;
stack1Roi = zeros(size(stack1));
stack1Roi(X,Y) = 1; 
%------------------------------------------------------------
% Draw manual ROI
% imagesc(stack(:,:,1));
% axis('square');
% axis('off');
% roi = roipoly;
%------------------------------------------------------------
figure('Color', 'w')
imagesc(stack1+ max(stack1, [], 'all').*(bwperim(stack1Roi)));
axis('square');
% hold off;
%-----------------------------------------------------------
% stackMod = stack.*uint16(stack1Roi);
% stackMod = stack(yAxis, xAxis, :);
%% -----------------------------------------------------------
% stackMod = zeros(size(stack));
el = 5;
se = strel('disk',el);
nucSize = 3;
blurControl = 0.65;
for i=1:size(stack, 3)
%     nucBin = Nucleus.seedNuclearMask3(stack(:,:,i), el, nucSize, blurControl, i, 100);
%     nucBin = Nucleus.seedNuclearMask3(stack(:,:,i), el, nucSize, blurControl, i, 100); % for gfp
%     nucBin = imclearborder(nucBin);
    nucBin = imbinarize(stack(:,:,i));
%     stackMod(:,:,i) = stack(:,:,i).*uint16(nucBin);
end
%%
filtWinSize = 15;
meanArr = 0;
varArr = 0;
for i=1:size(stackMod, 3)
%     imMean = conv2(stackMod(:,:,i), ones(filtWinSize)/filtWinSize^2, 'same');
    imMean = mean2(stackMod(:,:,i));
%     imStd =  stdfilt(stackMod(:,:,i), ones(filtWinSize));    
    imStd = std2(stackMod(:,:,i));
    meanArr = vertcat(meanArr, reshape(imMean, [numel(imMean), 1]));
    varArr = vertcat(varArr, reshape(imStd.^2, [numel(imStd), 1]));
end
scatter(meanArr(2:end), varArr(2:end), '.r');


[sortMean,indSort] = sort(meanArr);
sortVar = varArr(indSort);
[lB, uB]= bounds(meanArr(2:end));
binMeanInd = discretize(sortMean,linspace(lB, uB, 10));

binMean = accumarray(binMeanInd(2:end),sortMean(2:end))./accumarray(binMeanInd(2:end),1);
binVar = accumarray(binMeanInd(2:end),sortVar(2:end))./accumarray(binMeanInd(2:end),1);



end