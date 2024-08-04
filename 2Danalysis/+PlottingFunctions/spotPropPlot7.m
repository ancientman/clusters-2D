function spotPropPlot7(DSfolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots areas and intensities of spots inside and outside the global hotspots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gSpotDS = load(append(DSfolder, filesep, 'globalHotSpotDS.mat'));
propStruct = load(append(DSfolder, filesep, 'procConditionsDS.mat'));
corrStruct = load(append(DSfolder, filesep, 'corrDS.mat'));
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
Xpixel = propStruct.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
for i=1:nNuc
    areaInHotspotTL = (Xpixel)^2 .* vertcat(gSpotDS.globalSpotHotSpotStruct{i}.areaInHotspotNL);
    diaInHotspotTL = sqrt(4*areaInHotspotTL /pi);
    iMeanInHotspotTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.intensityAbsInHotspotNL);
    areaOutHotspotTL = (Xpixel)^2 .*vertcat(gSpotDS.globalSpotHotSpotStruct{i}.areaOutHotspotNL);
    diaOutHotspotTL = sqrt(4*areaOutHotspotTL /pi);
    iMeanOutHotspotTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.intensityAbsOutHotspotNL);
    iNucOutTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.nucIntensityOutNL);
    iNucInTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.nucIntensityInNL);
end
boxPlotter1(diaInHotspotTL, diaOutHotspotTL, 'light', 'aggSize');
boxPlotter1(iMeanInHotspotTL, iMeanOutHotspotTL, 'light', 'intensity');
end

function boxPlotter1(dataA, dataB, theme, type)
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
% colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];
colorPalete = [0.5, 0.8, 0.9; 0.9, 0.5, 0.7];
% hard coded plot parameters
x0=0;
y0=0;
plotWidth=300;
plotHeight=300;
% dataLim hardcoded
boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'ColorGroup',g, 'LabelOrientation', 'horizontal');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colorPalete(j,:),'FaceAlpha',.9);
end
hold on;
dataMean = [mean(dataA), mean(dataB)];
dataSem = [std(dataA)/sqrt(length(dataA)), std(dataB)/sqrt(length(dataB))];
errorbar(dataMean, dataSem, '.','MarkerSize',10,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
axes1 = gca;
set(axes1,'Color',[1 1 1],'LineWidth',1.5);
box(axes1,'on');
grid on;
if strcmp(type, 'aggSize')
    ylabel('Bcd agg. size ({\mu}m)');
    title('Mean aggregate size');
    
elseif strcmp(type, 'intensity')
    ylabel('Bcd agg. intensity (a. u)');
    title('Mean aggregate intensity');
end

textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center');
plotWidth=500;
plotHeight=500;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end 
