function outDS = corrPlotCombine(txtFilePath, type)
txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end

nucDS = cell(1, totalSubDirs);
gSpotDS = cell(1, totalSubDirs);
procDS = cell(1, totalSubDirs);
pairCorrDS = cell(1, totalSubDirs);
xCorrDS = cell(1, totalSubDirs);
Xpixel = zeros(1, totalSubDirs);
nNuc = zeros(1, totalSubDirs);

pairCorrTemp = [];
pairCorrLambda = [];

xCorrTempX = [];
xCorrTempY = [];
xCorrLambda = [];

xCorrTimeX = cell(1, 60);
xCorrTimeY = cell(1, 60);

for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        nucDS{i} = load(append(dataFiles.(fileID), filesep, 'nucDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        pairCorrDS{i} = load(append(dataFiles.(fileID), filesep, 'pairCorrDS.mat'));
        xCorrDS{i} = load(append(dataFiles.(fileID), filesep, 'xCorrDS.mat'));
        nNuc(i) = length(pairCorrDS{i}.pairCorrStruct.lambdaHalf);
        Xpixel(i) = procDS{i}.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
        for j=1:nNuc(i)
            pairCorrTemp = horzcat(pairCorrTemp, pairCorrDS{i}.pairCorrStruct.avg{j});
            pairCorrLambda = vertcat(pairCorrLambda, pairCorrDS{i}.pairCorrStruct.lambdaHalf{j});
            xCorrTempX = horzcat(xCorrTempX, xCorrDS{i}.crossCorrStruct.corrXAll(:,1));
            xCorrTempY = horzcat(xCorrTempY, xCorrDS{i}.crossCorrStruct.corrYAll(:,1));            
        end

        for t = 1:length(xCorrDS{i}.crossCorrStruct.intTimeCrossCorr)
            xCorrTimeX{t} = horzcat(xCorrTimeX{t}, horzcat(xCorrDS{i}.crossCorrStruct.intTimeCrossCorr{t}.corrX{:}));
            xCorrTimeY{t} = horzcat(xCorrTimeY{t}, horzcat(xCorrDS{i}.crossCorrStruct.intTimeCrossCorr{t}.corrY{:}));
        end
    end
end

xArrPair = (1:size(xCorrTimeX{t}, 1)).*Xpixel(1);
for t = 1:60
    corrTimeX{t} = mean(xCorrTimeX{t}, 2);
    corrTimeY{t} = mean(xCorrTimeY{t}, 2);
    xCorrTimeHalfLength(t) = patchPlot1(horzcat(corrTimeX{t},  corrTimeY{t}) ,xArrPair);
end


pairCorrMean = mean(pairCorrTemp, 2, 'omitnan');
pairCorrSem = std(pairCorrTemp, 1, 2, 'omitnan')./sqrt(size(pairCorrTemp, 2));
xArrPair = (1:size(pairCorrTemp, 1)).*Xpixel(1);
pairHalfLength = patchPlot1(pairCorrTemp,xArrPair);

xCorrAll = horzcat(xCorrTempY, xCorrTempX);
meanXCorrAll = mean(xCorrAll, 2, 'omitnan');
stdXCorrAll = std(xCorrAll, 0, 2, 'omitnan');
semXCorrAll = stdXCorrAll./sqrt(size(xCorrAll,2));
xArrCross = (1:size(xCorrTempX, 1)).*Xpixel(1);

xHalfLength = patchPlot1(xCorrAll,xArrCross);

outDS.pairCorrMat = pairCorrTemp;
outDS.xCorrXMat = xCorrTempX;
outDS.xCorrYMat = xCorrTempY;
outDS.corrTimeX = xCorrTimeX;
outDS.corrTimeY = xCorrTimeY;
outDS.xCorrTimeHalfLength = xCorrTimeHalfLength;
outDS.name = type;
procFolder = 'G:\Tyrone_analysis\2D\all2DCombine\all';
outDSName = (['allCorr', type]);
save([procFolder, filesep, outDSName, '.mat'],'outDS');
end

function halfLength = patchPlot1(dataMat, xArr)
halfLength = [];
figure('Color','w');
val = mean(dataMat, 2, 'omitnan');
err = std(dataMat, 1, 2, 'omitnan')./sqrt(size(dataMat, 2));

% ----------- Patch plot -----------------
patchTop = val+err;
patchTop = reshape(patchTop,1,[]);
patchBot = val-err;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch=[xArr,fliplr(xArr)];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = [0.6 0.6 0.6];
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.6;
hold on;
pl = plot(xArr, val);
pl.LineWidth = 1;
pl.Color = [0.4 0.4 0.4];
pl.LineStyle = '-';
hold on;
% ----------------------------------------

% ----------- Fit -----------------
hold on;
[~, limMin] = max(val(1:30));
lims = limMin:30;
tbl = table(xArr(lims)', val(lims));
modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));  
beta0 = [0.2, 5, 0.5]; 
    mdl = fitnlm(tbl, modelfun, beta0);
    coeffs= mdl.Coefficients{:, 'Estimate'};
    yFitted = coeffs(1) + coeffs(2) * exp(-coeffs(3)*xArr);
    halfLength = xArr(3)+ log(2)/coeffs(3);
    hold on;
    pf = plot(xArr, yFitted, 'k-', 'LineWidth', 1);
hold on;
% ---------------------------------
% ----------- Line -----------------
hold on;
pp = plot([0, max(xArr)], [1 1]);
pp.LineWidth = 1;
pp.Color = [0.8 0 0];
pp.LineStyle = '--';
% ----------------------------------
axes2 = gca;
ylabel('Paired correlation G(r)');
xlabel('x ({\mu}m)');
title('Paired ring correlation');
box(axes2,'on');
set(axes2,'Color',[1 1 1],'FontSize',10,'LineWidth',1);
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

ylim([0 10]);
xlim([0 3]);

hold off;
close;
end
