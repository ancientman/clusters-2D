function spotPropPlot8(DSfolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot histogram of intensities inside and outside the global hotspots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gSpotDS = load(append(DSfolder, filesep, 'globalHotSpotDS.mat'));
propStruct = load(append(DSfolder, filesep, 'procConditionsDS.mat'));
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
Xpixel = propStruct.conditionDS.imagingInfo.XpixelSize; %Scaling in nm
deltaT = propStruct.conditionDS.imagingInfo.timeResolution; %in seconds
for i=1:nNuc
    areaInHotspotTL = (Xpixel/1000)^2.*vertcat(gSpotDS.globalSpotHotSpotStruct{i}.areaInHotspotNL);
    iMeanInHotspotTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.intensityAbsInHotspotNL);
    areaOutHotspotTL = (Xpixel/1000)^2.*vertcat(gSpotDS.globalSpotHotSpotStruct{i}.areaOutHotspotNL);
    iMeanOutHotspotTL = vertcat(gSpotDS.globalSpotHotSpotStruct{i}.intensityAbsOutHotspotNL);
end
figure1 = figure('Color',[1 1 1]);
[minCombined, maxCombined] = bounds([iMeanInHotspotTL; iMeanOutHotspotTL]);
BE = linspace(minCombined, maxCombined, 10);
hIn = histogram(iMeanInHotspotTL, BE, 'Normalization', 'probability');
hold on;
pdIn = fitdist(iMeanInHotspotTL,'normal');
pdfIn = pdf(pdIn,iMeanInHotspotTL);
pdfIn = pdfIn*sum(hIn.Values*hIn.BinWidth); %pdf times area under the histogram curve
[iMeanInHotspotTL,idxMeanIn] = sort(iMeanInHotspotTL); %sort raw data, det indices
pdfIn = pdfIn(idxMeanIn); %sort y per those indices
pIn = plot(iMeanInHotspotTL,pdfIn);
hold on;
hOut = histogram(iMeanOutHotspotTL, BE, 'Normalization', 'probability');
hold on;
pdOut = fitdist(iMeanOutHotspotTL,'normal');
pdfOut = pdf(pdOut,iMeanOutHotspotTL);
pdfOut = pdfOut*sum(hOut.Values*hOut.BinWidth); %pdf times area under the histogram curve
[iMeanOutHotspotTL, idxMeanOut] = sort(iMeanOutHotspotTL); %sort raw data, det indices
pdfOut = pdfOut(idxMeanOut); %sort y per those indices
pOut = plot(iMeanOutHotspotTL,pdfOut);
%% Parameters from the pdf
% %%mean and sigma for normal
meanIn = pdIn.mu;
varIn = pdIn.sigma;
meanOut = pdOut.mu;
varOut = pdOut.sigma;
%% Plot properties
axes1 = gca;
ylabel('{\rho}');
xlabel({'Fluorescence intensity (a. u.)',''});
title('Mean Aggregate Intensity Distribution', 'Color',[0 0 0]);%'Color',[0.93 0.70 0.13]);
box(axes1,'on');
set(axes1,'Color',[1 1 1],'FontSize',12,'LineWidth',1.5,'XColor',...
[0 0 0],'YColor', [0 0 0]);
set(hIn,'Parent',axes1,'DisplayName','Inside enrichment zones','LineWidth',1.5,...
    'EdgeColor',[0 0 0], 'FaceColor',[0.5 0 0.5]);
set(pIn, 'Parent',axes1, 'LineWidth',2,'LineStyle','-.', 'Color',[1 0 1]);
set(hOut, 'Parent',axes1,'DisplayName','Outside enrichment zones',...
    'LineWidth',1.5,'EdgeColor',[0 0 0], 'FaceColor',[0 0.5 0.5]);
set(pOut, 'Parent',axes1,'LineWidth',2,'LineStyle','--', 'Color',[0 1 1]);
hold off;
%% Add legend with parameters
legend([pIn pOut],strcat('Inside enrichment zones {\mu} = ',num2str(round(meanIn,3))),...
    strcat('Outside enrichment zones {\mu}=',num2str(round(meanOut,3))),...
    'Textcolor',[0 0 0],'FontSize',12);
hold off;
ylim([0, 0.35]);
end
