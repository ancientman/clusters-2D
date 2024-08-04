function spotPropPlot4(DSfolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot number of spots per area per second in and outside hotspots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucDS = load(append(DSfolder, filesep, 'nucDS.mat'));
gSpotDS = load(append(DSfolder, filesep, 'globalHotSpotDS.mat'));
hotSpotLMat = load(append(DSfolder, filesep, 'hotSpotDS.mat'));
procDS = load(append(DSfolder, filesep, 'procConditionsDS.mat'));
Xpixel = procDS.conditionDS.imagingInfo.XpixelSize; %Scaling in nm
deltaT = procDS.conditionDS.imagingInfo.timeResolution; %in seconds
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
tFrames = length(gSpotDS.globalSpotHotSpotStruct{1}.timeUsedNL);
movMeanWindow = (1/deltaT);
nucArea = regionprops(nucDS.registerNucPropStruct{1}.labelMat, 'Area'); %here 1 is the time
nucArea = cell2mat(struct2cell(nucArea)).*((Xpixel)^2);
hotSpotArea = hotSpotLMat.globalHotSpotStruct.hotspotAreaNL.*((Xpixel)^2);

for i=1:nNuc
    if i==1
        if hotSpotArea(i) ~= 0
            spotPerAreaRateIn = movmean((gSpotDS.globalSpotHotSpotStruct{i}.inSpotCountNL)/...
                (hotSpotArea(i)),movMeanWindow);
            spotPerAreaRateOut = movmean((gSpotDS.globalSpotHotSpotStruct{i}.outSpotCountNL)/...
                (nucArea(i)-hotSpotArea(i)),movMeanWindow);


            molPerAreaRateIn = movmean((gSpotDS.globalSpotHotSpotStruct{i}.inMolCountNL)/...
                (hotSpotArea(i)),movMeanWindow);
            molPerAreaRateOut = movmean((gSpotDS.globalSpotHotSpotStruct{i}.outMolCountNL)/...
                (nucArea(i)-hotSpotArea(i)),movMeanWindow);
        else
            spotPerAreaRateIn = movmean((gSpotDS.globalSpotHotSpotStruct{i}.inSpotCountNL)...
                ,movMeanWindow);
            spotPerAreaRateOut = movmean((gSpotDS.globalSpotHotSpotStruct{i}.outSpotCountNL)/...
                (nucArea(i)),movMeanWindow);

            molPerAreaRateIn = movmean((gSpotDS.globalSpotHotSpotStruct{i}.inMolCountNL)...
                ,movMeanWindow);
            molPerAreaRateOut = movmean((gSpotDS.globalSpotHotSpotStruct{i}.outMolCountNL)/...
                (nucArea(i)),movMeanWindow);
        end
    else
        if hotSpotArea(i) ~= 0
            spotPerAreaRateIn = vertcat(spotPerAreaRateIn, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.inSpotCountNL)/...
                (hotSpotArea(i)),movMeanWindow));
            spotPerAreaRateOut = vertcat(spotPerAreaRateOut, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.outSpotCountNL)/...
                (nucArea(i)-hotSpotArea(i)),movMeanWindow));

            molPerAreaRateIn = vertcat(molPerAreaRateIn, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.inMolCountNL)/...
                (hotSpotArea(i)),movMeanWindow));
            molPerAreaRateOut = vertcat(molPerAreaRateOut, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.outMolCountNL)/...
                (nucArea(i)-hotSpotArea(i)),movMeanWindow));
        else
            spotPerAreaRateIn = vertcat(spotPerAreaRateIn, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.inSpotCountNL)...
                ,movMeanWindow));
            spotPerAreaRateOut = vertcat(spotPerAreaRateOut, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.outSpotCountNL)/...
                (nucArea(i)),movMeanWindow));

            molPerAreaRateIn = vertcat(molPerAreaRateIn, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.inMolCountNL)...
                ,movMeanWindow));
            molPerAreaRateOut = vertcat(molPerAreaRateOut, ...
                movmean((gSpotDS.globalSpotHotSpotStruct{i}.outMolCountNL)/...
                (nucArea(i)),movMeanWindow));
        end
    end
end

spotPerAreaRateIn = spotPerAreaRateIn/deltaT;
spotPerAreaRateOut = spotPerAreaRateOut/deltaT;

molPerAreaRateIn = molPerAreaRateIn/deltaT;
molPerAreaRateOut = molPerAreaRateOut/deltaT;

% histPlotter1(spotPerAreaRateIn, spotPerAreaRateOut); % For histograms
boxPlotter1(spotPerAreaRateIn, spotPerAreaRateOut , 'light', 'Rate of aggregate detection'); % for box plots
boxPlotter1(molPerAreaRateIn, molPerAreaRateOut , 'light', 'Rate of mol detection'); % for box plots
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
xlabel({'Rate ({\mu}m^{-2} . s^{-1})',''});
title('Rate of aggregate detection', 'FontWeight', 'bold');
box(axes1,'on');
set(axes1,'FontSize',12,'LineWidth',1.5);
% set(axes1,'Color',[0 0 0],'FontSize',12,'LineWidth',1.5,'XColor',...
% [0.93 0.70 0.13],'YColor', [0.93 0.70 0.13]);
set(hIn,'Parent',axes1,'DisplayName','Inside enrichment zones','LineWidth',1.5,...
    'EdgeColor',[0 0 0], 'FaceColor',[0.5, 0.8, 0.9]);
set(pIn, 'Parent',axes1, 'LineWidth',2,'LineStyle','-.', 'Color',[1 0 1]);
set(hOut, 'Parent',axes1,'DisplayName','Outside enrichment zones',...
    'LineWidth',1.5,'EdgeColor',[0 0 0], 'FaceColor',[0.5, 0.8, 0.9]);
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


function boxPlotter1(dataA, dataB, theme, label)
%% Boxplots of distribution fit of Spots/UnitArea/Time with wilcoxon ransum test values
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
set(axes1,'Color',[1 1 1],'LineWidth',1);
axes1.FontSize = 12;
% set(gcf,'position',[x0,y0,plotWidth,plotHeight])
ylabel({'Rate ({\mu}m^{-2} . s^{-1})',''});
title(label, 'FontWeight', 'bold');
textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center');
% hard coded figure prop
plotWidth=250;
plotHeight=250;
grid('off');
box('off');
x0=0;
y0=0;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end 
