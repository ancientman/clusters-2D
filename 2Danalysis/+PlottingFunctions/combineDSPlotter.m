function combineDSPlotter(plotDSFolder)
fileNames = what(plotDSFolder);
fileNames = fileNames.mat;
DS = cellfun(@(x) load(append(plotDSFolder, filesep, x)), fileNames, 'un', 0);
totalFiles = length(DS);
names = strings(totalFiles, 1);
diaIn = cell(totalFiles, 1);
diaOut = cell(totalFiles, 1);
iIn = cell(totalFiles, 1);
iOut = cell(totalFiles, 1);
nucIn = cell(totalFiles, 1);
nucOut = cell(totalFiles, 1);
rateIn = cell(totalFiles, 1);
rateOut = cell(totalFiles, 1);

sizeLimit = 0.1;

for i = 1:totalFiles
    names(i) = DS{i}.DS.type;
    
    temp = DS{i}.DS.diaInHotspotTL;
    temp(DS{i}.DS.diaInHotspotTL<sizeLimit) = [];
    diaIn{i} = temp;
    
    temp = DS{i}.DS.intensityInHotspotTL;
    temp(DS{i}.DS.diaInHotspotTL<sizeLimit) = [];
    iIn{i} = temp;
    
    temp = DS{i}.DS.intensityNucInTL;
    temp(DS{i}.DS.diaInHotspotTL<sizeLimit) = [];
    nucIn{i} = temp;
    
    temp = DS{i}.DS.diaOutHotspotTL;
    temp(DS{i}.DS.diaOutHotspotTL<sizeLimit) = [];
    diaOut{i} = temp;
    
    temp = DS{i}.DS.intensityOutHotspotTL;
    temp(DS{i}.DS.diaOutHotspotTL<sizeLimit) = [];
    iOut{i} = temp;
    
    temp = DS{i}.DS.intensityNucOutTL;
    temp(DS{i}.DS.diaOutHotspotTL<sizeLimit) = [];
    nucOut{i} = temp;
    
    rateIn{i} = DS{i}.DS.spotRateIn;
    rateOut{i} = DS{i}.DS.spotRateOut;
end

%_____________________________________________
%   Box plots
color = [0.5000    0.0000    0.5000];
figure('Color', 'w');
subplot(2, 3, 1);
boxPlotter2(rateIn, names, color, 'rate');
subplot(2, 3, 2)
boxPlotter2(diaIn, names,  color, 'size');
subplot(2, 3, 3)
boxPlotter2(iIn, names,  color, 'intensity');
%-----------------------------------------------------
color = [0.0000    0.5000    0.5000];
subplot(2, 3, 4);
boxPlotter2(rateOut, names, color, 'rate');
subplot(2, 3, 5);
boxPlotter2(diaOut, names,  color, 'size');
subplot(2, 3, 6);
boxPlotter2(iOut, names,  color, 'intensity');
hold off;
%_____________________________________________

%_____________________________________________
%   Histograms
% color = [0.5000    0.0000    0.5000];
% figure('Color', 'w');
% subplot(2, 3, 1);
% histPlotter2(rateIn, names, color, 'rate');
% legend(names, 'Location', 'nw');
% subplot(2, 3, 2);
% histPlotter2(diaIn, names, color, 'size');
% subplot(2, 3, 3);
% histPlotter2(iIn, names, color, 'intensity');
% %-----------------------------------------------------
% color = [0.0000    0.5000    0.5000];
% subplot(2, 3, 4);
% histPlotter2(rateOut, names, color, 'rate');
% legend(names,  'Location', 'nw');
% subplot(2, 3, 5);
% histPlotter2(diaOut, names, color, 'size');
% subplot(2, 3, 6);
% histPlotter2(iOut, names, color, 'intensity');
% hold off;
%_____________________________________________

%_____________________________________________
%   Plot vs nuclear intensity
color = [0.5000    0.0000    0.5000];
plotter2(iIn, nucIn, names, color, 'intensity');
plotter2(diaIn, nucIn, names, color, 'size');
color = [0.0000    0.5000    0.5000];
plotter2(iOut, nucOut, names, color, 'intensity');
plotter2(diaOut, nucOut, names, color, 'size');
%_____________________________________________
end

function plotter2(spotCell, nucCell, names, color, type)
% figure('color', 'w');
totalFiles = length(spotCell);
p = gobjects(1, totalFiles);
ff = gobjects(1, totalFiles);
for i = 1:totalFiles
    figure('color', 'w');
    f = polyfit(nucCell{i}, spotCell{i}, 1);
    x2 = linspace(min(nucCell{i}), max(nucCell{i}), 100);
    y2 = polyval(f,x2);
    p(i) = plot(nucCell{i}, spotCell{i});
    p(i).Marker = '.';
    p(i).Color = color;
    p(i).LineStyle = 'none';
    hold on;
    ff(i) = plot(x2, y2);
    ff(i).Color = 'k';
    ff(i).LineStyle = '--';
    ff(i).LineWidth = 1.5;
    %---------------------------------------------------------------------
%     legend(ff(i), num2str(f(2)), 'Location', 'nw');
    legend(strcat(type, ' (', names(i), ')'), 'Location', 'nw');
    legend('boxoff');
    %---------------------------------------------------------------------
    if strcmp(type, 'size')
        ylabel('Bcd agg. size ({\mu}m)');
    %     title('Mean aggregate size');    
        ylim([0, 0.4]);
    elseif strcmp(type, 'intensity')
        ylabel('Bcd agg. intensity (a. u)');
    %     title('Mean aggregate intensity');
        ylim([0, 250])
    elseif strcmp(type, 'rate')   
        ylabel({'Detection Frequency'; ' ({\mu}m^{-2} . s^{-1})'});
    %     title('Frequency of aggregate detection');
        ylim([0, 10]);
    end
    %---------------------------------------------------------------------
    hold on;
end

xlabel('Nuclear intensity (a. u.)');
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
end

function boxPlotter2(dataCell, names, color, type)
% figure('color', 'w');
dataA = dataCell{1};
dataB = dataCell{2};
nameA = names(1);
nameB = names(2);
dataCombine = [dataA; dataB];
g1 = repmat({nameA},length(dataA),1);
g2 = repmat({nameB},length(dataB),1);
g = [g1; g2];
% colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];
%------------------------------------------------------------------
% boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
%     'ColorGroup',g, 'LabelOrientation', 'horizontal', 'Symbol', '');
boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'LabelOrientation', 'horizontal', 'Symbol', '');
set(findobj(gca,'type','line'),'linew',2);
set(findobj(gcf, 'type', 'line', 'Tag', 'Median'),'Color', [0.7 0.7 0.7]);

set(findobj('-regexp','Tag','(Lower|Upper) (Whisker|Adjacent Value)'),'Color',[0.7, 0.7, 0.7]);

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color,'FaceAlpha',.1);
    set(h(j),'LineWidth',2);
    set(h(j),'MarkerSize',10);
    x = get(h(j),'XData');
    y = get(h(j),'YData');
    c = get(h(j),'Color');
    l = get(h(j),'LineWidth');
    ht = y(2)-y(1);
    wd = x(3)-x(1);
    rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',color,'LineWidth',l)
end
delete(h);
hold on;
dataMean = [mean(dataA), mean(dataB)];
dataSem = [std(dataA)/sqrt(length(dataA)), std(dataB)/sqrt(length(dataB))];
errorbar(dataMean, dataSem, '.','MarkerSize',10,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
hold on;
%---------------------------------------------------------------------
%   Scatter plot
% x1=ones(length(dataA), 1).*(1+(rand(length(dataA), 1)-0.5)/5);
% x2=ones(length(dataB), 1).*(1+(rand(length(dataB), 1)-0.5)/10);
% f1=scatter(x1, dataA, 'filled');
% f1.MarkerFaceColor = color;%colorPalete(2,:);
% f1.MarkerFaceAlpha = 0.2;
% hold on ;
% f2=scatter(x2.*2,dataB,'filled');
% f2.MarkerFaceColor = color;%colorPalete(1,:);
% f2.MarkerFaceAlpha = 0.2;
% hold on;
%---------------------------------------------------------------------
if strcmp(type, 'size')
    ylabel('Bcd agg. size ({\mu}m)');
%     title('Mean aggregate size');    
    ylim([0, 0.4]);
elseif strcmp(type, 'intensity')
    ylabel('Bcd agg. intensity (a. u)');
%     title('Mean aggregate intensity');
    ylim([0, 250])
elseif strcmp(type, 'rate')   
    ylabel({'Detection Frequency'; ' ({\mu}m^{-2} . s^{-1})'});
%     title('Frequency of aggregate detection');
    ylim([0, 10]);
end
%---------------------------------------------------------------------
textString = ['p = ', num2str(ranksum(dataA, dataB))];
text(1.5,double(max(dataCombine)),textString, ...
    'HorizontalAlignment','center', 'FontSize',14);
%---------------------------------------------------------------------
% hard coded plot parameters
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
end 

%______________________________________________________________________
function histPlotter2(dataCell, names, color, type)    % 
dataA = dataCell{1};
dataB = dataCell{2};
nameA = names(1);
nameB = names(2);
% colorPalete = [0.0000    0.5000    0.5000;  0.5000    0.0000    0.5000];
[NA,edgesA] = histcounts(dataA, 20, 'Normalization', 'probability');
[NB,edgesB] = histcounts(dataB, 20, 'Normalization', 'probability');
p1 = plot(edgesA(2:end), NA);
p1.Color = color; %colorPalete(1,:);
p1.LineWidth = 1.5;
p1.LineStyle = '-.';
hold on; 
p2 = plot(edgesB(2:end), NB);
p2.Color = color; %colorPalete(2,:);
p2.LineWidth = 1.5;
p2.LineStyle = '--';
hold on;
axes1 = gca;
%---------------------------------------------------------------------
if strcmp(type, 'size')
    ylabel('Normalized probability');
    xlabel('Bcd agg. size ({\mu}m)');
    xlim([0, 0.4]);
%     title('Mean aggregate size');
    
elseif strcmp(type, 'intensity')
    ylabel('Normalized probability');
    xlabel('Bcd agg. intensity (a. u)');
    xlim([0, 250])
%     title('Mean aggregate intensity');
elseif strcmp(type, 'rate')
    ylabel('Normalized probability');
    xlabel({'Detection frequency'; ' ({\mu}m^{-2} . s^{-1})'});
    xlim([0, 10]);
%     title('Frequency of aggregate detection', 'FontWeight', 'bold');
end
%---------------------------------------------------------------------
box(axes1,'on');
set(axes1,'FontSize',12,'LineWidth',1.5);
set(gca, 'XScale', 'log')
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold on;
% legend([pIn pOut],strcat('Inside enrichment zones {\mu} = ', num2str(meanIn)), ...
%     strcat('Outside enrichment zones {\mu} = ',num2str(meanOut)), ...
%     'Textcolor',[0 0 0],'FontSize',12, 'Location', 'northeast');
% ylim([0, 0.15]);
% hold off;
end