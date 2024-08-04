function compareCorr(txtFilePath, type)

txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end

% color = [102,194,165;...
%     252,141,98;...
%     141,160,203;...
%     231,138,195;
%     89, 157, 212];

color = [107,174,214; ...
    66,146,198; ...
    33,113,181; ...
    8,81,156; ...
    8,48,107];

% color = [1,1,230;...
%     30,30,230;...
%     80,80,230;...
%     130,130,230;
%     180,180,230];

procDS = cell(1, totalSubDirs);
corrDS = cell(1, totalSubDirs);
Xpixel = zeros(1, totalSubDirs);
corrMean = cell(1, totalSubDirs);
corrStd = cell(1, totalSubDirs);
corrSem = cell(1, totalSubDirs);
f1 = figure('color', 'w'); hold on;
f2 = figure('color', 'w'); hold on;
hold on;
for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        Xpixel(i) = procDS{i}.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
        corrDS{i} = load(append(dataFiles.(fileID), filesep, 'xCorrDS.mat'));
        corrMean{i} = mean(horzcat(corrDS{i}.crossCorrStruct.corrXAll, corrDS{i}.crossCorrStruct.corrYAll), 2);
        corrStd{i} = std(horzcat(corrDS{i}.crossCorrStruct.corrXAll, corrDS{i}.crossCorrStruct.corrYAll), 0, 2);
        corrSem{i} = corrStd{i}./sqrt(size(horzcat(corrDS{i}.crossCorrStruct.corrXAll, corrDS{i}.crossCorrStruct.corrYAll), 2));
        xArr = Xpixel(i).*[1:length(corrMean{i})];
        set(0, 'CurrentFigure', f1)
        dataLim = 20;
        [halfLengthCross{i}, pl1(i)] = patchPlot1(corrMean{i}, corrStd{i}, xArr, dataLim, color(i,:));
        ylim([0 1]);
        xlim([0 1]);    
        hold on;
    end
end
name = append({'133 ms', '266 ms', '512 ms', '1064 ms'});
legend(pl1, name, "Box", "off");
ylim([-0.05 inf])
set(0, 'CurrentFigure', f1)
axes('Position',[.5 .3 .4 .2])
box on
plotBar(halfLengthCross, color)
ylim([0.1 0.3])
end

function plotBar(dataCell, color)
color = color./255;
mat = horzcat(dataCell{:});
data = mat(1,:);
err = mat(2,:);
b = bar(1:length(dataCell),data);
for i = 1:length(b)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = 'flat';
    b(i).FaceAlpha = 0.4;
    b(i).BarWidth = 0.6;
    b(i).LineStyle = 'none';
    for j = 1:length(dataCell)
        b(i).CData(j,:) = color(j,:); % Color for first data coloumn
    end
end

hold on

er = errorbar(1:length(data), data, err); 
er.Color = [0.1 0.1 0.1];                            
er.LineStyle = 'none';
er.LineWidth = 1;
er.CapSize = 0;
ylabel('{\mu} m')
end

function [halfLength, pt] = patchPlot1(val, err, xArr, dataLim, color)
color = color./255;
% val = mean(dataMat, 2, 'omitnan');
% err = std(dataMat, 1, 2, 'omitnan')./sqrt(size(dataMat, 2));

% err = std(dataMat, 1, 2, 'omitnan');

% % ----------- Patch plot -----------------
patchCropIdx = min(find(isnan(err)));
if isempty(patchCropIdx)
    patchCropIdx = length(xArr);
end

patchTop = val(1:patchCropIdx-1)+err(1:patchCropIdx-1);
patchTop = reshape(patchTop,1,[]);
patchBot = val(1:patchCropIdx-1)-err(1:patchCropIdx-1);
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch=[xArr(1:patchCropIdx-1),fliplr(xArr(1:patchCropIdx-1))];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = color;
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.3;
hold on;
pl = plot(xArr, val);
pl.LineWidth = 1;
pl.Color = color;
pl.LineStyle = '-';
hold on;
%% -----------------------------------------

% hold on;
% pl = errorbar(xArr, val, err);
% pl.LineWidth = 1;
% pl.Color = color;
% pl.LineStyle = '-';
% pl.CapSize = 0;
% hold on;
% ----------------------------------------

% ----------- Fit -----------------
hold on;
[~, limMin] = max(val(1:dataLim));
lims = limMin:dataLim;
tbl = table(xArr(lims)', val(lims));
modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));  
beta0 = [0.2, 5, 0.5]; 
mdl = fitnlm(tbl, modelfun, beta0);
coeffs= mdl.Coefficients{:, 'Estimate'};
yFitted = coeffs(1) + coeffs(2) * exp(-coeffs(3)*xArr);
halfLength = xArr(3)+ log(2)/coeffs(3);
halfLengthStd = halfLength*(mdl.Coefficients.SE(3))/coeffs(3);
halfLength = [halfLength; halfLengthStd];
hold on;
% pf = plot(xArr, yFitted);
% pf.LineWidth = 1;
% pf.Color = color;
% pf.LineStyle = ':';
% hold on;
% ---------------------------------
% ----------- Line -----------------
% hold on;
% pp = plot([0, max(xArr)], [1 1]);
% pp.LineWidth = 1;
% pp.Color = [0 0 0];
% pp.LineStyle = '--';
% ----------------------------------
axes2 = gca;
ylabel('Pixel cross correlation');
xlabel('x ({\mu}m)');
box(axes2,'on');
set(axes2,'Color',[1 1 1],'FontSize',10,'LineWidth',1);
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
xlim([0 1]);
hold off;
end