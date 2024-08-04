function spotPropPlot4AllCombine(txtFileFolder)
files=dir(txtFileFolder);
filepath = strings(length(files)-2, 1);
name = strings(length(files)-2, 1);
ext = strings(length(files)-2, 1);
txtFilePath = strings(length(files)-2, 1);
for i = 1:length(files)-2
    [filepath(i),name(i),ext(i)] = fileparts(strcat(files(i+2).folder, filesep, files(i+2).name));
    txtFilePath(i) = fullfile(strcat(files(i+2).folder, filesep, files(i+2).name));
end

txtFilePathChar = convertStringsToChars(txtFilePath);
info = cellfun(@(x) HelperFunctions.readtext(x,'\t'), txtFilePathChar,'un', 0);
totalSubDirs = cellfun(@(x) size(x, 1), info, 'un', 0);
dataFiles = cellfun(@(x) cell2struct(x(:,2), x(:,1),1), info, 'un', 0);


%%%%%%%%% incomplete beyond this%%%%%%%%%%%%%%
nucDS = cell(1, totalSubDirs);
gSpotDS = cell(1, totalSubDirs);
hotSpotLMat = cell(1, totalSubDirs);
procDS = cell(1, totalSubDirs);
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

for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        nucDS{i} = load(append(dataFiles.(fileID), filesep, 'nucDS.mat'));
        gSpotDS{i} = load(append(dataFiles.(fileID), filesep, 'globalHotSpotDS.mat'));
        hotSpotLMat{i} = load(append(dataFiles.(fileID), filesep, 'hotSpotDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        deltaT(i) = procDS{i}.conditionDS.imagingInfo.timeResolution; %in seconds
        Xpixel(i) = procDS{i}.conditionDS.imagingInfo.XpixelSize; %Scaling in microns
        nNuc(i) = length(gSpotDS{i}.globalSpotHotSpotStruct);
        tFrames(i) = length(gSpotDS{i}.globalSpotHotSpotStruct{1}.timeUsedNL);
        movMeanWindow(i) = 10; 
        
        nucArea{i} = regionprops(nucDS{i}.registerNucPropStruct{1}.labelMat, 'Area'); %here 1 is the time
        nucArea{i} = cell2mat(struct2cell(nucArea{i})).*((Xpixel(i))^2);
        hotSpotArea{i} = hotSpotLMat{i}.globalHotSpotStruct.hotspotAreaNL.*((Xpixel(i))^2);
        
        for j=1:nNuc(i)
            if hotSpotArea{i}(j) ~= 0
                spotPerAreaRateInTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.inSpotCountNL)/...
                    (hotSpotArea{i}(j)),movMeanWindow(i));
                spotPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outSpotCountNL)/...
                    (nucArea{i}(j)-hotSpotArea{i}(j)),movMeanWindow(i));
            else
                spotPerAreaRateInTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.inSpotCountNL)...
                    ,movMeanWindow(i));
                spotPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outSpotCountNL)/...
                    (nucArea{i}(j)),movMeanWindow(i));
            end
            
            if spotPerAreaRateInTL == 0
                spotPerAreaRateInTL = spotPerAreaRateInTemp;
            else
                spotPerAreaRateInTL = vertcat(spotPerAreaRateInTL, spotPerAreaRateInTemp);
            end

            if spotPerAreaRateOutTL == 0
                spotPerAreaRateOutTL = spotPerAreaRateOutTemp;
            else
                spotPerAreaRateOutTL = vertcat(spotPerAreaRateOutTL, spotPerAreaRateOutTemp);
            end

            spotPerAreaRateInNL = vertcat(spotPerAreaRateInNL, mean(spotPerAreaRateInTemp));
            spotPerAreaRateOutNL = vertcat(spotPerAreaRateOutNL, mean(spotPerAreaRateOutTemp));
        end
    end
end

spotPerAreaRateInTL = spotPerAreaRateInTL/(deltaT(2));
spotPerAreaRateOutTL = spotPerAreaRateOutTL/(deltaT(2));

spotPerAreaRateInNL = spotPerAreaRateInNL/(deltaT(2));
spotPerAreaRateOutNL = spotPerAreaRateOutNL/(deltaT(2));
boxPlotter1(spotPerAreaRateInNL, spotPerAreaRateOutNL , 'light'); % for box plots
end

function boxPlotter1(dataA, dataB, theme)
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

x1=ones(length(dataA), 1).*(1+(rand(length(dataA), 1)-0.5)/5);
x2=ones(length(dataB), 1).*(1+(rand(length(dataB), 1)-0.5)/10);
f1=scatter(x1, dataA, 'filled');
f1.MarkerFaceColor = colorPalete(2,:);
f1.MarkerFaceAlpha = 0.2;
hold on ;
f2=scatter(x2.*2,dataB,'filled');
f2.MarkerFaceColor = colorPalete(1,:);
f2.MarkerFaceAlpha = 0.2;
hold on;



ylabel({'Rate ({\mu}m^{-2} . s^{-1})',''});
title('Rate of aggregate detection', 'FontWeight', 'bold');
textString = ['p = ', num2str(ranksum(dataA, dataB))];%, '%2.3fe%05d')];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center', 'FontSize',12);

ax = gca;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=500;
plotHeight=500;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;

ylim([0,8]);
end 

