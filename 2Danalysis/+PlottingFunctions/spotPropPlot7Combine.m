function outDS = spotPropPlot7Combine(txtFilePath, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots areas and intensities of spots inside and outside the global hotspots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end

gSpotDS = cell(1, totalSubDirs);
procDS = cell(1, totalSubDirs);
corrDS = cell(1, totalSubDirs);
Xpixel = zeros(1, totalSubDirs);
nNuc = zeros(1, totalSubDirs);

for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        gSpotDS{i} = load(append(dataFiles.(fileID), filesep, 'globalHotSpotDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
%         corrDS{i} = load(append(dataFiles.(fileID), filesep, 'corrDS.mat'));
        nNuc(i) = length(gSpotDS{i}.globalSpotHotSpotStruct);
        Xpixel(i) = procDS{i}.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
        for j=1:nNuc(i)
            areaInHotspotTLTemp = (Xpixel(i))^2 .* (gSpotDS{i}.globalSpotHotSpotStruct{j}.areaInHotspotNL);        
            iMeanInHotspotTLTemp = (gSpotDS{i}.globalSpotHotSpotStruct{j}.intensityAbsInHotspotNL);
            areaOutHotspotTLTemp = (Xpixel(i))^2 .*(gSpotDS{i}.globalSpotHotSpotStruct{j}.areaOutHotspotNL);        
            iMeanOutHotspotTLTemp = (gSpotDS{i}.globalSpotHotSpotStruct{j}.intensityAbsOutHotspotNL);
            iNucOutTLTemp = vertcat(gSpotDS{i}.globalSpotHotSpotStruct{j}.nucIntensityOutNL);
            iNucInTLTemp = vertcat(gSpotDS{i}.globalSpotHotSpotStruct{j}.nucIntensityInNL);
            
            if i==1 && j==1
                areaInHotspotTL = areaInHotspotTLTemp;        
                iMeanInHotspotTL = iMeanInHotspotTLTemp;
                areaOutHotspotTL = areaOutHotspotTLTemp;        
                iMeanOutHotspotTL = iMeanOutHotspotTLTemp;
                
                iNucOutTL =  iNucOutTLTemp;
                iNucInTL = iNucInTLTemp;
                
                areaInHotspotNL = mean(areaInHotspotTLTemp);        
                iMeanInHotspotNL = mean(iMeanInHotspotTLTemp);
                areaOutHotspotNL = mean(areaOutHotspotTLTemp);        
                iMeanOutHotspotNL = mean(iMeanOutHotspotTLTemp);
            else            
                areaInHotspotTL = vertcat(areaInHotspotTL, areaInHotspotTLTemp);     
                iMeanInHotspotTL = vertcat(iMeanInHotspotTL, iMeanInHotspotTLTemp);
                areaOutHotspotTL = vertcat(areaOutHotspotTL, areaOutHotspotTLTemp);            
                iMeanOutHotspotTL =vertcat(iMeanOutHotspotTL, iMeanOutHotspotTLTemp);   
                
                iNucOutTL =  vertcat(iNucOutTL, iNucOutTLTemp);
                iNucInTL =  vertcat(iNucInTL, iNucInTLTemp);
                areaInHotspotNL = vertcat(areaInHotspotNL, mean(areaInHotspotTLTemp));        
                iMeanInHotspotNL = vertcat( iMeanInHotspotNL, mean(iMeanInHotspotTLTemp));
                areaOutHotspotNL = vertcat(areaOutHotspotNL, mean(areaOutHotspotTLTemp));        
                iMeanOutHotspotNL = vertcat(iMeanOutHotspotNL, mean(iMeanOutHotspotTLTemp));
            end
        end
    end        
end
diaInHotspotTL = sqrt(4*areaInHotspotTL /pi);
diaOutHotspotTL = sqrt(4*areaOutHotspotTL /pi);

diaInHotspotNL = sqrt(4*areaInHotspotNL /pi);
diaOutHotspotNL = sqrt(4*areaOutHotspotNL /pi);

outDS.diaInHotspotTL = diaInHotspotTL;
outDS.diaOutHotspotTL = diaOutHotspotTL;
outDS.intensityInHotspotTL = iMeanInHotspotTL;
outDS.intensityOutHotspotTL = iMeanOutHotspotTL;
outDS.intensityNucOutTL = iNucOutTL;
outDS.intensityNucInTL = iNucInTL;
% outDSName = (['inOutDS', type]);
% procFolder = 'E:\Dropbox (Princeton)\Data\Analysis_RM\DS for plots';
% save([procFolder, filesep, outDSName, '.mat'],'outDS');


boxPlotter2(diaInHotspotTL, diaOutHotspotTL, 'light', 'aggSize');
boxPlotter2(iMeanInHotspotTL, iMeanInHotspotTL, 'light', 'intensity');
boxPlotter2(diaInHotspotNL, diaOutHotspotNL, 'light', 'aggSize');
boxPlotter2(iMeanInHotspotNL, iMeanOutHotspotNL, 'light', 'intensity');

% boxPlotter1(diaInHotspotNL, diaOutHotspotNL, 'light', 'aggSize');
% boxPlotter1(iMeanInHotspotNL, iMeanOutHotspotNL, 'light', 'intensity');
end

function boxPlotter2(dataA, dataB, theme, type)
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
%------------------------------------------------------------------
boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'ColorGroup',g, 'LabelOrientation', 'horizontal', 'Symbol', '');
set(findobj(gca,'type','line'),'linew',1);
set(findobj(gcf, 'type', 'line', 'Tag', 'Median'),'Color', [0.7 0.7 0.7]);

set(findobj('-regexp','Tag','(Lower|Upper) (Whisker|Adjacent Value)'),'Color',[0.7, 0.7, 0.7]);

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colorPalete(j,:),'FaceAlpha',.1);
    set(h(j),'LineWidth',1);
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

% hold on;
% 
% x1=ones(length(dataA), 1).*(1+(rand(length(dataA), 1)-0.5)/5);
% x2=ones(length(dataB), 1).*(1+(rand(length(dataB), 1)-0.5)/10);
% f1=scatter(x1, dataA, 'filled');
% f1.MarkerFaceColor = colorPalete(2,:);
% f1.MarkerFaceAlpha = 0.2;
% hold on ;
% f2=scatter(x2.*2,dataB,'filled');
% f2.MarkerFaceColor = colorPalete(1,:);
% f2.MarkerFaceAlpha = 0.2;
hold on;
%---------------------------------------------------------------------
if strcmp(type, 'aggSize')
    ylabel('Mean cluster size ({\mu}m)');
%     title('Mean aggregate size');
    
elseif strcmp(type, 'intensity')
    ylabel('Mean cluster intensity (a. u)');
%     title('Mean aggregate intensity');
end

textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center', 'FontSize',14);
% hard coded plot parameters
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;

set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
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
colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];
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
if strcmp(type, 'aggSize')
    ylabel('Bcd agg. size ({\mu}m)');
    title('Mean aggregate size');
    
elseif strcmp(type, 'intensity')
    ylabel('Bcd agg. intensity (a. u)');
    title('Mean aggregate intensity');
end

textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center');
plotWidth=200;
plotHeight=200;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end 
