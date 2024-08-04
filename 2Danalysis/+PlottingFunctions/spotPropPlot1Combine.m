function outDS = spotPropPlot1Combine(txtFilePath, type)
txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end

nucDS = cell(1, totalSubDirs);
spotDS = cell(1, totalSubDirs);
hotSpotLMat = cell(1, totalSubDirs);
procDS = cell(1, totalSubDirs);
corrDS = cell(1, totalSubDirs);
nucArea =  cell(1, totalSubDirs);
hotSpotArea =  cell(1, totalSubDirs);
deltaT = zeros(1, totalSubDirs);
Xpixel = zeros(1, totalSubDirs);
nNuc = zeros(1, totalSubDirs);
tFrames = zeros(1, totalSubDirs);
movMeanWindow = zeros(1, totalSubDirs);

spotPerAreaRateInTL = 0;
spotPerAreaRateOutTL = 0;

spotPerAreaRateInNL = 0;
spotPerAreaRateOutNL = 0;

meanSpotIntensity = cell(1, totalSubDirs);
meanNucIntensity = cell(1, totalSubDirs);
meanIntensitySub = cell(1, totalSubDirs);
meanIntensityRatio = cell(1, totalSubDirs);
intensityRatioNuc = cell(1, totalSubDirs);
intensitySubNuc = cell(1, totalSubDirs);
intensityNuc = cell(1, totalSubDirs);
intensitySpot = cell(1, totalSubDirs);

for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        nucDS{i} = load(append(dataFiles.(fileID), filesep, 'nucDS.mat'));
        spotDS{i} = load(append(dataFiles.(fileID), filesep, 'spotDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        deltaT(i) = procDS{i}.conditionDS.imagingInfo.timeResolution; %in seconds
        Xpixel(i) = procDS{i}.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
        nNuc(i) = size(nucDS{i}.registerNucPropStruct{1}.pixValue, 2);
        tFrames(i) = length(spotDS{i}.globalRegSpotStruct);
        
        meanSpotIntensity{i} = cell(1, nNuc(i));
        meanNucIntensity{i} = cell(1, nNuc(i));
        meanIntensitySub{i} = cell(1, nNuc(i));
        meanIntensityRatio{i} = cell(1, nNuc(i));
        intensityRatioNuc{i} = cell(1, nNuc(i));
        intensitySubNuc{i} = cell(1, nNuc(i));
        intensityNuc{i} = cell(1, nNuc(i));
        intensitySpot{i} = cell(1, nNuc(i));
        
        for k = 1:nNuc(i)        
            meanSpotIntensity{i}{k} = cell(1, tFrames(i));
            meanNucIntensity{i}{k} = cell(1, tFrames(i));
            meanIntensitySub{i}{k} = cell(1, tFrames(i));
            meanIntensityRatio{i}{k} = cell(1, tFrames(i));
            for t=1:tFrames(i)
                meanIntensitySpotsTemp = vertcat(vertcat(spotDS{i}.globalRegSpotStruct{t}.intensityAbsInHotspot{k}), ...
                    vertcat(spotDS{i}.globalRegSpotStruct{t}.intensityAbsOutHotspot{k}));
                meanIntensityNucTemp = mean(nucDS{i}.registerNucPropStruct{t}.pixValue{k});
                meanSpotIntensity{i}{k}{t} = mean(meanIntensitySpotsTemp);
                meanNucIntensity{i}{k}{t} = mean(meanIntensityNucTemp);
                meanIntensitySub{i}{k}{t} = (meanSpotIntensity{i}{k}{t} - meanNucIntensity{i}{k}{t});      
                meanIntensityRatio{i}{k}{t} = (meanSpotIntensity{i}{k}{t} - meanNucIntensity{i}{k}{t})./meanNucIntensity{i}{k}{t};            
            end
            ratioTemp = vertcat(meanIntensityRatio{i}{k}{:});
            subTemp = vertcat(meanIntensitySub{i}{k}{:});
            nucTemp = vertcat(meanNucIntensity{i}{k}{:});
            spotTemp = vertcat(meanSpotIntensity{i}{k}{:});
            if k == 1
                intensityRatioNuc{i} = ratioTemp;
                intensitySubNuc{i} = subTemp;
                intensityNuc{i} = nucTemp;
                intensitySpot{i} = spotTemp;
            else
                intensityRatioNuc{i} = horzcat(intensityRatioNuc{i}, ratioTemp);
                intensitySubNuc{i} = horzcat(intensitySubNuc{i}, subTemp);
                intensityNuc{i} = horzcat(intensityNuc{i}, nucTemp);
                intensitySpot{i} = horzcat(intensitySpot{i}, spotTemp);
            end
        end
    end
    intensityRatioNuc{i} = reshape(intensityRatioNuc{i}, [size(intensityRatioNuc{i}, 1)*size(intensityRatioNuc{i}, 2), 1]);
    intensitySubNuc{i} = reshape(intensitySubNuc{i}, [size(intensitySubNuc{i}, 1)*size(intensitySubNuc{i}, 2), 1]);
    intensityNuc{i} = reshape(intensityNuc{i}, [size(intensityNuc{i}, 1)*size(intensityNuc{i}, 2), 1]);
    intensitySpot{i} = reshape(intensitySpot{i}, [size(intensitySpot{i}, 1)*size(intensitySpot{i}, 2), 1]);
end

intensityNuc = vertcat(intensityNuc{:});
[intensityNuc, idx] = sort(intensityNuc);

intensitySpot = vertcat(intensitySpot{:});
intensitySpot = intensitySpot(idx);

intensityRatioNuc = vertcat(intensityRatioNuc{:});
intensityRatioNuc = intensityRatioNuc(idx);

intensitySubNuc = vertcat(intensitySubNuc{:});
intensitySubNuc = intensitySubNuc(idx);

nBins = 11;
[~ , nucInEdges] = discretize(intensityNuc, nBins);

spotNucBin = cell(nBins-1,1);
ratioNucBin = cell(nBins-1,1);
subNucBin = cell(nBins-1,1);

for i = 1:length(intensityNuc)
    for j = 1:nBins-1
        if intensityNuc(i)>nucInEdges(j) && intensityNuc(i)<=nucInEdges(j+1)
            spotNucBin{j} = vertcat(spotNucBin{j}, intensitySpot(i));
            subNucBin{j} = vertcat(subNucBin{j}, intensitySubNuc(i));
            ratioNucBin{j} = vertcat(ratioNucBin{j}, intensityRatioNuc(i));
            continue;
        end
    end
end


countSubNucBin = cellfun(@length, subNucBin);

meanSpotNucBin = cellfun(@(x) mean(x, 'omitnan'), spotNucBin);
semSpotNucBin = cellfun(@(x) std(x, 0, 'omitnan'), spotNucBin);
semSpotNucBin = semSpotNucBin./sqrt(vertcat(countSubNucBin));

meanSubNucBin = cellfun(@(x) mean(x, 'omitnan'), subNucBin);
semSubNucBin = cellfun(@(x) std(x, 0, 'omitnan'), subNucBin);
semSubNucBin = semSubNucBin./sqrt(vertcat(countSubNucBin));

meanRatioNucBin = cellfun(@(x) mean(x, 'omitnan'), ratioNucBin);
semRatioNucBin = cellfun(@(x) std(x, 0, 'omitnan'), ratioNucBin);
semRatioNucBin = semRatioNucBin./sqrt(vertcat(countSubNucBin));


figure('Color', 'w'); cla; hold on;
errorbar(nucInEdges(1:nBins-1), meanRatioNucBin, semRatioNucBin, '-v','MarkerSize',10,...
    'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerFaceColor','w', 'Color', [0.3 0.3 0.3], 'linewidth', 2, 'LineStyle', 'none')
hold on;
f = fit(nucInEdges(1:nBins-1)', meanRatioNucBin, 'exp1', 'Exclude', [1 2]);%vertcat(xCellNan{:})<0.2);
pf = plot(f);
pf(1).Color = [0.3, 0.3, 0.3];
pf.LineStyle = '--';
pf.LineWidth = 2;
hold on;
title('Capicua');
ylabel('Relative enrichment in clusters');
xlabel('Mean nuclear concentration (a. u.)');
ax = gca;
ax.FontSize = 12;
% ax.XLim = [0.1, 0.5];
% ax.YLim = [8, 25];
ax.XGrid = 1;
ax.YGrid = 1;
ax.LineWidth = 1.5;
ax.Box = 1;
grid off;
hold off;

x0 = 100;
y0= 100;
plotWidth=500;
plotHeight=500;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;

outDS.spotRateIn = spotPerAreaRateInTL;
outDS.spotRateOut = spotPerAreaRateOutTL;

% procFolder = 'E:\Dropbox (Princeton)\Data\Analysis_RM\DS for plots';
% outDSName = (['inOutDS', type]);
% save([procFolder, filesep, outDSName, '.mat'],'outDS');
% boxPlotter2(spotPerAreaRateInTL, spotPerAreaRateOutTL , 'light'); % for box plots

% histPlotter1(spotPerAreaRateInNL, spotPerAreaRateOutNL); % For histograms
% boxPlotter1(spotPerAreaRateInTL, spotPerAreaRateOutTL , 'light'); % for box plots
% boxPlotter1(spotPerAreaRateInNL, spotPerAreaRateOutNL , 'light'); % for box plots
end

function boxPlotter1(dataA, dataB, theme)
% Boxplots of distribution fit of Spots/UnitArea/Time with wilcoxon ranksum test values
if strcmp(theme, 'light')
    figure1 = figure('Color', [1 1 1]);
elseif strcmp(theme, 'dark')
        figure1 = figure('Color', [0 0 0]);
else
    figure1 = figure('Color', [1 1 1]);
end
dataCombine = [dataA; dataB];
g1 = repmat({'Inside'},length(dataA),1);
g2 = repmat({'Outside'},length(dataB),1);
g = [g1; g2];
colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];
% hard coded plot parameters
x0=0;
y0=0;
plotWidth=300;
plotHeight=300;
boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'ColorGroup',g, 'LabelOrientation', 'horizontal');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colorPalete(j,:),'FaceAlpha',.5);
end
hold on;
dataMean = [mean(dataA), mean(dataB)];
dataSem = [std(dataA)/sqrt(length(dataA)), std(dataB)/sqrt(length(dataB))];
errorbar(dataMean, dataSem, '.','MarkerSize',10,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
axes1 = gca;
box(axes1,'on');
grid on;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
ylabel({'Rate ({\mu}m^{-2} . s^{-1})',''});
title('Frequency of aggregate detection', 'FontWeight', 'bold');
textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center');
plotWidth=200;
plotHeight=200;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end 

function boxPlotter2(dataA, dataB, theme)
% Boxplots of distribution fit of Spots/UnitArea/Time with wilcoxon ranksum test values
if strcmp(theme, 'light')
    figure('Color', [1 1 1]);
elseif strcmp(theme, 'dark')
        figure('Color', [0 0 0]);
else
    figure('Color', [1 1 1]);
end

dataCombine = [dataA; dataB];
g1 = repmat({'Inside'},length(dataA),1);
g2 = repmat({'Outside'},length(dataB),1);
g = [g1; g2];

colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];

boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'ColorGroup',g, 'LabelOrientation', 'horizontal', 'Symbol', '');
set(findobj(gca,'type','line'),'linew',2);
set(findobj(gcf, 'type', 'line', 'Tag', 'Median'),'Color', [0.7 0.7 0.7]);

set(findobj('-regexp','Tag','(Lower|Upper) (Whisker|Adjacent Value)'),'Color',[0.7, 0.7, 0.7]);

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colorPalete(j,:),'FaceAlpha',.1);
    set(h(j),'LineWidth',2);
    set(h(j),'MarkerSize',10);
    x = get(h(j),'XData');
    y = get(h(j),'YData');
    c = get(h(j),'Color');
    l = get(h(j),'LineWidth');
    ht = y(2)-y(1);
    wd = x(3)-x(1);
    rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',c,'LineWidth',l)
end
delete(h);
hold on;

dataMean = [mean(dataA), mean(dataB)];
dataSem = [std(dataA)/sqrt(length(dataA)), std(dataB)/sqrt(length(dataB))];
errorbar(dataMean, dataSem, '.','MarkerSize',10,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
hold on;

% x1=ones(length(dataA), 1).*(1+(rand(length(dataA), 1)-0.5)/5);
% x2=ones(length(dataB), 1).*(1+(rand(length(dataB), 1)-0.5)/10);
% f1=scatter(x1, dataA, 'filled');
% f1.MarkerFaceColor = colorPalete(2,:);
% f1.MarkerFaceAlpha = 0.2;
% hold on ;
% f2=scatter(x2.*2,dataB,'filled');
% f2.MarkerFaceColor = colorPalete(1,:);
% f2.MarkerFaceAlpha = 0.2;
% hold on;



ylabel({'Frequency of aggregate detection'; ' ({\mu}m^{-2} . s^{-1})'});
% title('Rate of aggregate detection', 'FontWeight', 'bold');
% textString = ['p = ', num2str(ranksum(dataA, dataB))];%, '%2.3fe%05d')];
% text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center', 'FontSize',12);

ax = gca;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;

ylim([0,8]);
end 


function histPlotter1(dataA, dataB)    
%% Histogram and distribution fit of Spots/UnitArea/Time
figure2 = figure('Color', [1 1 1]);% = figure('Color', [0 0 0]);
% [minCombined, maxCombined] = bounds([spotPerAreaRateIn; spotPerAreaRateOut]);
% BE = linspace(minCombined, maxCombined, 20);
hIn = histogram(dataA, 'Normalization', 'probability');
hold on; 
PdIn = fitdist(dataA,'normal');
PdfIn = pdf(PdIn,dataA);
PdfIn = PdfIn*sum(hIn.Values * hIn.BinWidth); %pdf times area under the histogram curve
[dataA, idxIn] = sort(dataA); %sort raw data, det indices
PdfIn = PdfIn(idxIn); %sort y per those indices
pIn = plot(dataA,PdfIn,'-', 'linewidth', 2);
hold on;
hOut = histogram(dataB, 'Normalization', 'probability');
hold on; 
PdOut = fitdist(dataB,'normal');
PdfOut = pdf(PdOut,dataB);
PdfOut = PdfOut*sum(hOut.Values * hOut.BinWidth); %pdf times area under the histogram curve
[dataB, idxOut] = sort(dataB); %sort raw data, det indices
PdfOut = PdfOut(idxOut); %sort y per those indices
pOut = plot(dataB,PdfOut,'-', 'linewidth', 2);
hold off;
%% Parameters from the pdf
% %%mean and sigma for normal
meanIn = PdIn.mu;
varIn = PdIn.sigma;
meanOut = PdOut.mu;
varOut = PdOut.sigma;
%% Plot properties
axes1 = gca;
ylabel('{\rho}');
xlabel({'Frequency of aggregate detection'; ' ({\mu}m^{-2} . s^{-1})'});
title('Frequency of aggregate detection', 'FontWeight', 'bold');
box(axes1,'on');
set(axes1,'FontSize',12,'LineWidth',1.5);
% set(axes1,'Color',[0 0 0],'FontSize',12,'LineWidth',1.5,'XColor',...
% [0.93 0.70 0.13],'YColor', [0.93 0.70 0.13]);
set(hIn,'Parent',axes1,'DisplayName','Inside enrichment zones','LineWidth',1.5,...
    'EdgeColor',[0 0 0], 'FaceColor',[0.5 0 0.5]);
set(pIn, 'Parent',axes1, 'LineWidth',2,'LineStyle','-.', 'Color',[1 0 1]);
set(hOut, 'Parent',axes1,'DisplayName','Outside enrichment zones',...
    'LineWidth',1.5,'EdgeColor',[0 0 0], 'FaceColor',[0 0.5 0.5]);
set(pOut, 'Parent',axes1,'LineWidth',2,'LineStyle','-.', 'Color',[0 1 1]);
set(gca, 'XScale', 'log')
grid on;
hold off;
legend([pIn pOut],strcat('Inside enrichment zones {\mu} = ', num2str(meanIn)), ...
    strcat('Outside enrichment zones {\mu} = ',num2str(meanOut)), ...
    'Textcolor',[0 0 0],'FontSize',12, 'Location', 'northeast');
ylim([0, 0.15]);
hold off;
end