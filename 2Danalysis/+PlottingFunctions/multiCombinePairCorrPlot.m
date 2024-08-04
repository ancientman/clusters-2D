function multiCombinePairCorrPlot(folderPath)
% Combines pair correlation data from different fly lines
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 
color = [102,194,165;252,141,98;141,160,203;231,138,195];

halfLengthPair = cell(1, length(struct2C));
halfLengthCross = cell(1, length(struct2C));

f1 = figure('color', 'w'); hold on;
f2 = figure('color', 'w'); hold on;
for i = 1:length(struct2C)
    name{i} = struct2C{i}.outDS.name;
%     corrMean{i} = mean(struct2C{i}.outDS.corrMat, 2, 'omitnan');
%     corrStd{i} = std(struct2C{i}.outDS.corrMat, 0, 2, 'omitnan');
%     corrSem{i} = corrStd{i}./sqrt(size(struct2C{1}.outDS.corrMat, 2));
    xArr = 0.043.*(1:size(struct2C{i}.outDS.pairCorrMat, 1));
    set(0, 'CurrentFigure', f1)
    dataLim = 20;
    if ~isempty(struct2C{i}.outDS.pairCorrMat)
        [halfLengthPair{i}, pl1(i)] = patchPlot1(struct2C{i}.outDS.pairCorrMat,xArr, dataLim, color(i,:));
        legName(i) = name(i); 
    end
    ylim([0 10]);
    xlim([0 1]);    
    hold on;
    set(0, 'CurrentFigure', f2)
    dataLim = 20;
    xArr = 0.043.*(1:size(struct2C{i}.outDS.xCorrXMat, 1));
    [halfLengthCross{i}, pl2(i)] = patchPlot1(horzcat(struct2C{i}.outDS.xCorrXMat, struct2C{i}.outDS.xCorrYMat),xArr, dataLim, color(i,:));
    ylim([-0.1 1]);
    xlim([0 0.5]);
    ylabel('Spatial cross correlation')
    title('')
    hold on;
end

[halfLengthPSF, pl2(end+1)] = AnalysisSim.psf();
xlim([0 0.5]);
name{end+1} = 'PSF';
halfLengthCross{end+1} = halfLengthPSF;

for i=1:length(struct2C)
if ~isempty(struct2C{i}.outDS.pairCorrMat)
    legend(pl1(i), legName(i), "Box", "off");
end
end
legend(pl1(:), legName(:), "Box", "off");

legend(pl2, name, "Box", "off");

set(0, 'CurrentFigure', f1)
xlim([0 1])
hold on;
pp = plot([0, max(xArr)], [1 1]);
pp.LineWidth = 1;
pp.Color = [0 0 0];
pp.LineStyle = '--';

set(0, 'CurrentFigure', f2)
axes('Position',[.5 .3 .4 .2])
box on
plotBar(halfLengthCross, color)
ylim([0.15 inf])
ylabel('{\mu}m')
% set(0, 'CurrentFigure', f1)
% hold on;
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
end

function [halfLength, pl] = patchPlot1(dataMat, xArr, dataLim, color)
color = color./255;
val = mean(dataMat, 2, 'omitnan');
err = std(dataMat, 1, 2, 'omitnan')./sqrt(size(dataMat, 2));

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
pt.FaceAlpha = 0.6;
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
ylabel('Paired autocorrelation G(r)');
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
xlim([0 1]);
hold off;
end
