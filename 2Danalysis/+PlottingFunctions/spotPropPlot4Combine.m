function outDS = spotPropPlot4Combine(txtFilePath, type)
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

spotPerAreaRateInTL = [];
molPerAreaRateInTL = [];
spotPerAreaRateOutTL = [];
molPerAreaRateOutTL = [];

spotPerAreaRateInNL = [];
molPerAreaRateInNL = [];
spotPerAreaRateOutNL = [];
molPerAreaRateOutNL = [];

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
%         corrDS{i} = load(append(dataFiles.(fileID), filesep, 'corrDS.mat'));
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

                molPerAreaRateInTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.inMolCountNL)/...
                    (hotSpotArea{i}(j)),movMeanWindow(i));
                
                iNucOutTLTemp = unique(gSpotDS{i}.globalSpotHotSpotStruct{j}.nucIntensityOutNL);
                iNucInTLTemp = unique(gSpotDS{i}.globalSpotHotSpotStruct{j}.nucIntensityInNL);
            
            
                spotPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outSpotCountNL)/...
                    (nucArea{i}(j)-hotSpotArea{i}(j)),movMeanWindow(i));

                molPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outMolCountNL)/...
                    (nucArea{i}(j)-hotSpotArea{i}(j)),movMeanWindow(i));
            else  
                spotPerAreaRateInTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.inSpotCountNL)...
                    ,movMeanWindow(i));
                
                molPerAreaRateInTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.inMolCountNL)...
                    ,movMeanWindow(i));
                
                spotPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outSpotCountNL)/...
                    (nucArea{i}(j)),movMeanWindow(i));

                molPerAreaRateOutTemp = movmean((gSpotDS{i}.globalSpotHotSpotStruct{j}.outMolCountNL)/...
                    (nucArea{i}(j)),movMeanWindow(i));
            end
                        
%             if i ==1 && j==1
%                 spotPerAreaRateInTL = spotPerAreaRateInTemp;
%                 spotPerAreaRateOutTL = spotPerAreaRateOutTemp;
%                 
%                 spotPerAreaRateInNL = mean(spotPerAreaRateInTemp);
%                 spotPerAreaRateOutNL = mean(spotPerAreaRateOutTemp);
%             else
%                 spotPerAreaRateInTL = vertcat(spotPerAreaRateInTL, spotPerAreaRateInTemp);
%                 spotPerAreaRateOutTL = vertcat(spotPerAreaRateOutTL, spotPerAreaRateOutTemp);
%                 
%                 spotPerAreaRateInNL = vertcat(spotPerAreaRateInNL, mean(spotPerAreaRateInTemp));
%                 spotPerAreaRateOutNL = vertcat(spotPerAreaRateOutNL, mean(spotPerAreaRateOutTemp));
%             end
            
            if spotPerAreaRateInTL == 0
                spotPerAreaRateInTL = spotPerAreaRateInTemp;       
                molPerAreaRateInTL = molPerAreaRateInTemp;     
%                 iNucInTL = iNucInTLTemp;
                
            else
                spotPerAreaRateInTL = vertcat(spotPerAreaRateInTL, spotPerAreaRateInTemp);
                molPerAreaRateInTL = vertcat(molPerAreaRateInTL, molPerAreaRateInTemp);
%                 iNucInTL =  vertcat(iNucInTL, iNucInTLTemp);
            end

            if spotPerAreaRateOutTL == 0
%                 iNucOutTL =  iNucOutTLTemp;
                spotPerAreaRateOutTL = spotPerAreaRateOutTemp;
                molPerAreaRateOutTL = molPerAreaRateOutTemp;
            else
                spotPerAreaRateOutTL = vertcat(spotPerAreaRateOutTL, spotPerAreaRateOutTemp);
                molPerAreaRateOutTL = vertcat(molPerAreaRateOutTL, molPerAreaRateOutTemp);
%                 iNucOutTL =  vertcat(iNucOutTL, iNucOutTLTemp);
            end
            spotPerAreaRateInNL = vertcat(spotPerAreaRateInNL, mean(spotPerAreaRateInTemp));
            molPerAreaRateInNL = vertcat(molPerAreaRateInNL, mean(molPerAreaRateInTemp));
            spotPerAreaRateOutNL = vertcat(spotPerAreaRateOutNL, mean(spotPerAreaRateOutTemp));       
            molPerAreaRateOutNL = vertcat(molPerAreaRateOutNL, mean(molPerAreaRateOutTemp));       
        end
    end
end

spotPerAreaRateInTL = spotPerAreaRateInTL./(deltaT(1));
molPerAreaRateInTL = molPerAreaRateInTL./(deltaT(1));
spotPerAreaRateOutTL = spotPerAreaRateOutTL./(deltaT(1));
molPerAreaRateOutTL = molPerAreaRateOutTL./(deltaT(1));


spotPerAreaRateInNL = nonzeros(spotPerAreaRateInNL)./(deltaT(1));
molPerAreaRateInNL = nonzeros(molPerAreaRateInNL)./(deltaT(1));
spotPerAreaRateOutNL = nonzeros(spotPerAreaRateOutNL)./(deltaT(1));
molPerAreaRateOutNL = nonzeros(molPerAreaRateOutNL)./(deltaT(1));

outDS.spotRateIn = spotPerAreaRateInTL;
outDS.molRateIn = molPerAreaRateInTL;
outDS.spotRateOut = spotPerAreaRateOutTL;
outDS.molRateOut = molPerAreaRateOutTL;

% procFolder = 'E:\Dropbox (Princeton)\Data\Analysis_RM\DS for plots';
% outDSName = (['inOutDS', type]);
% save([procFolder, filesep, outDSName, '.mat'],'outDS');
boxPlotter2(spotPerAreaRateInTL, spotPerAreaRateOutTL , 'light'); % for box plots
% boxPlotter2(molPerAreaRateInTL, molPerAreaRateOutTL , 'light'); % for box plots
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
ylabel({'Detection rate ({\mu}m^{-2} . s^{-1})',''});
title('Local maxima detection rate', 'FontWeight', 'bold');
textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center');
plotWidth=250;
plotHeight=250;
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
dataStd = [std(dataA), std(dataB)];
errorbar(dataMean, dataStd, '.','MarkerSize',10,...
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



ylabel({'Detection rate ({\mu}m^{-2} . s^{-1})'});
% title('Rate of aggregate detection', 'FontWeight', 'bold');
% textString = ['p = ', num2str(ranksum(dataA, dataB))];%, '%2.3fe%05d')];
% text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center', 'FontSize',12);

ax = gca;
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
xlabel({'Local maxima detection rate ({\mu}m^{-2} . s^{-1})'});
title('Local maxima detection rate', 'FontWeight', 'bold');
box(axes1,'on');
set(axes1,'FontSize',10,'LineWidth',1);
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