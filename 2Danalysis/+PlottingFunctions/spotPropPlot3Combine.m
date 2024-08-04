function spotPropPlot3Combine(txtFilePath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots hotsopt size distribbtion
% Plot hotspot occupancy: cumulative occupancy probability change with
% time, and mean interval between between consequtive occupancies.
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
hotDS = cell(1, totalSubDirs);
procDS = cell(1, totalSubDirs);
deltaT = zeros(1, totalSubDirs);
nNuc = zeros(1, totalSubDirs);
tFrames = zeros(1, totalSubDirs);
cumTime= cell(1, totalSubDirs);
realTime = cell(1, totalSubDirs);
onTimes = cell(1, totalSubDirs);
offTimes = cell(1, totalSubDirs);
movMeanWindow = zeros(1, totalSubDirs);
kPoint = 7; % points for mean values
data = 0;
hsArea = 0;
hsDia = 0;

for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        hotDS{i} = load(append(dataFiles.(fileID), filesep, 'hotSpotDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        for q = 1:length(hotDS{i}.globalHotSpotStruct.hotspotAreaUniq)
            if ~isempty(hotDS{i}.globalHotSpotStruct.hotspotAreaUniq{q})
                dia = procDS{i}.conditionDS.imagingInfo.XpixelSize.*...
                    sqrt((4/pi).*hotDS{i}.globalHotSpotStruct.hotspotAreaUniq{q});
                hsDia = vertcat(hsDia, dia);
            end
        end  
        
    end
end
figure('Color', 'w')
[NA,edgesA] = histcounts(nonzeros(hsDia), 20, 'Normalization', 'probability');
patchX = [edgesA(2:end), fliplr(edgesA(2:end))];
patchY = [NA, 0.1.*zeros(1, length(NA))];
p1 = patch(patchX, patchY, [0.7 0.7 0.7]);
p1.FaceAlpha = 0.5;
p1.EdgeColor = [0.1 0.1 0.1]; 
p1.LineWidth = 1;
p1.LineStyle = '-';
hold on; 
%---------------------------------------------------------------------
gca;
xlabel({'Confinement dia. ({\mu}m)'});
ylabel('Normalized probability');
box(gca, 'on');
set(gca,'FontSize',10,'LineWidth',1);
% set(gca, 'XScale', 'log')
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
%---------------------------------------------------------------------


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i=1:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        gSpotDS{i} = load(append(dataFiles.(fileID), filesep, 'globalHotSpotDS.mat'));
        hotDS{i} = load(append(dataFiles.(fileID), filesep, 'hotSpotDS.mat'));
        procDS{i} = load(append(dataFiles.(fileID), filesep, 'procConditionsDS.mat'));
        deltaT(i) = procDS{i}.conditionDS.imagingInfo.timeResolution; %in seconds
        nNuc(i) = length(gSpotDS{i}.globalSpotHotSpotStruct);
        for q = 1:length(hotDS{i}.globalHotSpotStruct.hotspotAreaUniq)
            if ~isempty(hotDS{i}.globalHotSpotStruct.hotspotAreaUniq{q})
                hsArea = vertcat(hsArea, hotDS{i}.globalHotSpotStruct.hotspotAreaUniq{q});
            end
        end  
        tFrames(i) = length(gSpotDS{i}.globalSpotHotSpotStruct{1}.timeUsedNL);
        movMeanWindow(i) = ceil(1/deltaT(i));
        k = cell(1, nNuc(i));
        kk = cell(1, nNuc(i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot the moving average of occupancy events
        cumTime{i} = zeros(tFrames(i), nNuc(i));
        realTime{i} = zeros(tFrames(i), nNuc(i));
        onTimes{i} = cell(1, nNuc(i));
        offTimes{i} = cell(1, nNuc(i));
        for j=1:nNuc(i)
            temp1 = 0;
            for t=1:tFrames(i)
                realTime{i}(t,j) = gSpotDS{i}.globalSpotHotSpotStruct{j}.timeUsedNL(t);
                temp1 = temp1 + (movmean(gSpotDS{i}.globalSpotHotSpotStruct{j}.timeUsedNL(t),kPoint));    
                cumTime{i}(t, j) = temp1;                    
            end
            propsOn = regionprops(logical(realTime{i}(:,j)), 'Area');
            propsOff = regionprops(~logical(realTime{i}(:,j)), 'Area');
            % counts = cat([props.Area]);
            onTimes{i}{j} = ([propsOn.Area]);
            offTimes{i}{j} = ([propsOff.Area]);
            temp2 = movmean(gSpotDS{i}.globalSpotHotSpotStruct{j}.timeUsedNL,kPoint);
            k{j} = find(gSpotDS{i}.globalSpotHotSpotStruct{j}.timeUsedNL);
            p = 2;
            while p<length(k{j})
                kk{j}(p-1, 1) = k{j}(p)-k{j}(p-1);
                p = p+1;
            end
        end
        onTimeAll{i} = horzcat(onTimes{i}{:});
        offTimeAll{i} = horzcat(offTimes{i}{:});
    end
    if i==1
        data = cumTime{1}(1:size(cumTime{1}, 1), :);
         kkk = kk{1};
         for jj=2:nNuc(i)-1 
            kkk = vertcat(kkk, kk{jj+1});
        end
    else
        data = horzcat(data, cumTime{i}(1:size(cumTime{1}, 1), :));
        for jj=1:nNuc(i)-1 
            kkk = vertcat(kkk, kk{jj+1});
        end
    end 
end
onTimeAll = horzcat(onTimeAll{:})';
offTimeAll = horzcat(offTimeAll{:})';
edges = 0:1:15;
[onCount, onEdges] = histcounts(onTimeAll.*deltaT(1), edges, 'Normalization','probability');
[offCount, offEdges] = histcounts(offTimeAll.*deltaT(1), edges, 'Normalization','probability');


figure('color', 'w');
[axOn, tauOn] = plotFit(onEdges(2:end), log(onCount), [1 0 1]);
ylabel('log (counts)');
xlabel('T_{on} (s)');
hold on;
[axOff, tauOff] = plotFit(offEdges(2:end), log(offCount), [0.3 0.3 0.3]);
ylabel('log (counts)');
xlabel('T_{off} (s)');
fracOn = (onTimeAll./(onTimeAll+offTimeAll));
[fracCount, fracEdges] = histcounts(fracOn, 'Normalization','probability');
figure('color', 'w')
x = fracOn;
histogram(x,'norm','pdf')
hold on;
xscale = (x - mean(x));
xscale = xscale/max([abs(min(xscale)),max(xscale)])/2 + 1/2;
[fracP] = betafit(xscale);
hold on
fplot(@(X) betapdf(X,fracP(1),fracP(2)),[0 1]);





%% Cumulative probability of enrichment zone occupancy
cumulativeEnrichmentFactor(deltaT(1), tFrames(1), data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interlval of consecutive aggregate  detection in enrichment zones
% on times
enrichmentHistogram(onTimeAll, deltaT(1));
% off times
enrichmentHistogram(kkk*deltaT(1), deltaT(1));
end

function cumulativeEnrichmentFactor(deltaT, tFrames, data)
figure1 = figure('Color',[1 1 1]);
xAxis = deltaT(1).*linspace(0,tFrames(1), tFrames(1));
dataMean = mean(data,2)/(tFrames(1));
dataError =  std(data,0,2)/(tFrames(1));
 x_vector = [xAxis fliplr(xAxis)];
patch = fill(x_vector, [(dataMean+dataError)' fliplr((dataMean-dataError)')], [0.7 0.7 0.7]);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
hold on;
p2 = plot(xAxis, dataMean, 'color', [0.3 0.3 0.3], 'LineWidth', 1);
hold on;
p3 = plot(deltaT(1)*linspace(0,tFrames(1)), linspace(0,1), 'k-.');
% hold on;
% errorbar(deltaT(1).*linspace(0,tFrames(1), tFrames(1)), dataMean, dataError'.');
% xlim([0 deltaT(1)*tFrames(1)]);
figure(figure1);
axes2 = gca;
ylabel('Cumulative probability');
xlabel({'Time (s)',''});
title('Probability of an enrichment zone to be occupied');
box(axes2,'on');
set(axes2,'Color',[1 1 1],'FontSize',10,'LineWidth',1);
set(p2, 'Parent',axes2, 'LineWidth',1, 'Color',[0.3 0.3 0.3]);
set(p3, 'Parent',axes2, 'LineWidth',1,'LineStyle','-.', 'Color',[0.4 0.4 0.4]);
legend(p3,'y = x', 'FontSize',10, 'Location', 'northwest', "Box", "off");
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function enrichmentHistogram(kkk, deltaT)
figure2 = figure('Color',[1 1 1]);
hIn = histogram(kkk, 'Normalization', 'probability', 'BinWidth', deltaT(1));
hold on; 
PdIn = fitdist(kkk,'exponential');
PdfIn = pdf(PdIn,kkk);
PdfIn = PdfIn*sum(hIn.Values * hIn.BinWidth); %pdf times area under the histogram curve
[kkk, idx] = sort(kkk); %sort raw data, det indices
PdfIn = PdfIn(idx); %sort y per those indices
p2 = plot(kkk,PdfIn,'-', 'linewidth', 1);
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
figure(figure2);
axes1 = gca;
ylabel('Normalized probability');
xlabel({'Time intervals (s)',''});
title('Time interval until an enrichment zone is occupied again');
box(axes1,'on');
set(axes1,'Color',[1 1 1],'FontSize',10,'LineWidth',1);
set(hIn,'Parent',axes1,'DisplayName','intervals','LineWidth',1,...
    'EdgeColor',[0.4 0.4 0.4], 'FaceColor', [ 0.7 0.7 0.7]);
set(p2, 'Parent',axes1, 'LineWidth',1,'LineStyle','-.', 'Color',[0.3 0.3 0.3]);
ylim([0 1]);
grid off;
box off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function [ax, tauFit] = plotFit(xVal, yVal, color)
xVal(isnan( yVal ) | isinf( yVal )) = [];
yVal(isnan( yVal ) | isinf( yVal )) = [];
s1 = scatter(xVal, yVal, 'filled');
s1.MarkerEdgeAlpha = 0;
s1.MarkerFaceColor = color;
s1.MarkerFaceAlpha = 0.6;
s1.SizeData = 8;
ph = s1;
hold on;
dataTable = table(xVal', yVal');
mdl = fitlm(dataTable);
rSq = mdl.Rsquared.Adjusted;
hold on;
fitVar = mdl.Coefficients.Variables;
lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
slope = [fitVar(2, 1), fitVar(2,2)];

xArr = linspace(0, max(xVal));
yArr = fitVar(2, 1).*xArr + fitVar(1, 1);
hold on;
pl = plot(xArr, yArr);
pl.LineStyle = '--';
pl.LineWidth = 1;
pl.Color = color;
hold on;
ph = pl;

lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
slope = [fitVar(2, 1), fitVar(2,2)];
tauFit = lambda;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w'); 
end
