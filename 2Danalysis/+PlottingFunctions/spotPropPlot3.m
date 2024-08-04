function spotPropPlot3(DSfolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot hotspot occupancy: cumulative occupancy probability change with
% time, and mean interval between between consequtive occupancies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gSpotDS = load(append(DSfolder, filesep, 'globalHotSpotDS.mat'));
propStruct = load(append(DSfolder, filesep, 'procConditionsDS.mat'));
Xpixel = propStruct.conditionDS.imagingInfo.XpixelSize; %Scaling in nm
deltaT = propStruct.conditionDS.imagingInfo.timeResolution; %in seconds
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
tFrames = length(gSpotDS.globalSpotHotSpotStruct{1}.timeUsedNL);
movMeanWindow = ceil(1/deltaT);
k = cell(1, nNuc);
kk = cell(1, nNuc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the moving average of occupancy events
cumTime = zeros(tFrames, nNuc);
kPoint = 7; % points for mean values
for i=1:nNuc
    temp = 0;
    for t=1:tFrames
        temp = temp + (movmean(gSpotDS.globalSpotHotSpotStruct{i}.timeUsedNL(t),kPoint));    
        cumTime(t, i) = temp;
    end    
end
%% Cumulative probability of enrichment zone occupancy
figure1 = figure('Color',[1 1 1]);
xAxis = deltaT.*linspace(0,tFrames, tFrames);
data = mean(cumTime,2)/(tFrames);
errors =  std(cumTime,0,2)/(tFrames);
 x_vector = [xAxis fliplr(xAxis)];
patch = fill(x_vector, [(data+errors)' fliplr((data-errors)')], [0.9290, 0.6940, 0.1250]);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
hold on;
p2 = plot(xAxis, data, 'color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);
hold on;
p3 = plot(deltaT*linspace(0,tFrames), linspace(0,1), 'k-.');
% hold on;
% errorbar(deltaT.*linspace(0,tFrames, tFrames), mean(cumTime,2)/(tFrames), std(cumTime,0,2)/(tFrames)'.');
% xlim([0 deltaT*tFrames]);
figure(figure1);
axes2 = gca;
ylabel('Cumulative probability ({\rho})');
xlabel({'Time (s)',''});
title('Probability of an enrichment zone to be occupied');
box(axes2,'on');
set(axes2,'Color',[1 1 1],'LineWidth',1.5);
set(p2, 'Parent',axes2, 'LineWidth',2, 'Color',[0.2, 0.2, 0.2]);%][0, 0.4470, 0.7410]);
set(p3, 'Parent',axes2, 'LineWidth',2,'LineStyle','-.', 'Color',[0.5 0.5 0.5]);
legend(p3,'y = x', 'Location', 'northwest');
grid on;
x0=0;
y0=0;
plotWidth = 300;
plotHeight = 300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interlval of consecutive aggregate  detection in enrichment zones
figure2 = figure('Color',[1 1 1]);
for i=1:nNuc    
    p1 = plot(movmean(gSpotDS.globalSpotHotSpotStruct{i}.timeUsedNL,movMeanWindow));
%     plot(gSpotDS.spotHotSpotStruct{i}.timeUsedNL, '-o');    
%     stem(deltaT*linspace(0,tFrames), gSpotDS.globalSpotHotSpotStruct{i}.timeUsedNL, 'Marker', '.', 'Color', 'k')
    hold on;
    k{i} = find(gSpotDS.globalSpotHotSpotStruct{i}.timeUsedNL);
    p = 2;
    while p<length(k{i})
        kk{i}(p-1, 1) = k{i}(p)-k{i}(p-1);
        p = p+1;
    end
end
hold off;
kkk = kk{1};
for i=2:nNuc-1 
    kkk = vertcat(kkk, kk{i+1});
end
kkk = kkk.*deltaT;
hIn = histogram(kkk, 'Normalization', 'probability');
hold on; 
PdIn = fitdist(kkk,'exponential');
PdfIn = pdf(PdIn,kkk);
PdfIn = PdfIn*sum(hIn.Values * hIn.BinWidth); %pdf times area under the histogram curve
[kkk, idx] = sort(kkk); %sort raw data, det indices
PdfIn = PdfIn(idx); %sort y per those indices
p2 = plot(kkk,PdfIn,'-', 'linewidth', 2);
hold off;
figure(figure2);
axes1 = gca;
ylabel('Probability ({\rho})');
xlabel({'time until the next  event (s)',''});
title('Time interval until an enrichment zone is occupied again');
box(axes1,'on');
set(axes1,'Color',[1 1 1],'FontSize',12,'LineWidth',1.5);
set(hIn,'Parent',axes1,'DisplayName','intervals','LineWidth',1.5,...
    'EdgeColor',[0, 0.4470, 0.7410], 'FaceColor', [0.9290, 0.6940, 0.1250]);
set(p2, 'Parent',axes1, 'LineWidth',2,'LineStyle','-.', 'Color',[0, 0.4470, 0.7410]);
ylim([0 1]);
hold off;
end

