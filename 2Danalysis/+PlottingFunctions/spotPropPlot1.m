function spotPropPlot1(DSfolder, plotOn)
%______________________________________________________________________________________
%
%	Plots the change in nuclear intensity with time
%
%______________________________________________________________________________________
gSpotDS = load(append(DSfolder, filesep, 'globalHotSpotDS.mat'));
hSpotDS = load(append(DSfolder, filesep, 'hotSpotDS.mat'));
nucPropDS = load(append(DSfolder, filesep, 'nucDS.mat'));
propStruct = load(append(DSfolder, filesep, 'procConditionsDS.mat'));
deltaT = propStruct.conditionDS.imagingInfo.timeResolution; %in seconds
nNuc = length(gSpotDS.globalSpotHotSpotStruct);
tFrames = length(gSpotDS.globalSpotHotSpotStruct{1}.timeUsedNL);

f1 = @(x) getfield(x, 'pixValue');
nucValues = cellfun(f1, nucPropDS.registerNucPropStruct, 'un', 0);

f1 = @(x) getfield(x, 'pixIdxList');
nucIdx = cellfun(f1, nucPropDS.registerNucPropStruct, 'un', 0);
nucIdx = vertcat(nucIdx{:});
nucVal = vertcat(nucValues{:});

hsIdx = hSpotDS.globalHotSpotStruct.hotspotIdxListUniq;
hsIdx = repmat(hsIdx, tFrames, 1);

[~,hsElements, ~] = cellfun(@intersect, nucIdx, hsIdx, 'un', 0);

hsVal = cell(tFrames, nNuc);
for i = 1:numel(hsElements)
    hsVal{i} = nucVal{i}(hsElements{i});
end
hsValMean = cellfun(@mean, hsVal, 'un', 0);
hsValMeanT = cellfun(@(x) vertcat(x{:}), num2cell(hsValMean, 1), 'un', 0);
hsValMeanNormT = cellfun(@(x) x/(x(1)), hsValMeanT, 'un', 0);
hsValMeanNormTMean = mean(horzcat(hsValMeanNormT{:}), 2);%          cellfun(@(x) x/(x(1)), hsValMeanT, 'un', 0);
hsValTNormSem = std(horzcat(hsValMeanNormT{:}), 0, 2)./sqrt(nNuc);
nucMeanT = cell(tFrames,1);
timeCell = cell(tFrames,nNuc);
for t = 1:tFrames
    nucMeanT{t} = cellfun(@mean, nucValues{t}, 'un', 0);
    timeCell(t,:) = cellfun(@(x) (t*deltaT), timeCell(t,:),'un', 0);
end
nucValuesTime = cellfun(@vertcat, (nucMeanT{:,1}), 'un', 0);
nucValuesTNorm = cellfun(@(x) (x/x(1)), nucValuesTime, 'un', 0);
nucValTNormMean =  mean(horzcat(nucValuesTNorm{:}),2);
nucValTNormSem = std(horzcat(nucValuesTNorm{:}), 0, 2)./sqrt(nNuc);

if (strcmp(plotOn, 'yes'))
    plotAllBleach(hsValMeanNormTMean, hsValTNormSem, timeCell, 'rich');
    plotAllBleach(nucValTNormMean, nucValTNormSem, timeCell, 'nuc')

    %   Plot individual nuclei
    %_______________________________________________________________________________
    % figure('Color' ,'w');
    % for i = 1:nNuc
    %     hold on;
    %     fitCell = cellfun(@(x, y) fit(x , y, 'exp1'), {vertcat(timeCell{:,i})}, nucValuesTNorm(:,i), 'un', 0);  
    %     p1 = cellfun(@(x, y) plot(x, y), {vertcat(timeCell{:,i})}, nucValuesTNorm(:,i));
    %     p1.Marker = 'o';
    %     p1.MarkerSize = 6;
    %     p1.MarkerFaceColor = [0.8 0.8 0.8];
    %     p1.MarkerEdgeColor = [0.3 0.3 0.3];
    %     p1.LineWidth = 0.5;
    %     p1.LineStyle = 'none';
    %     
    %     hold on;
    %     f1 = cellfun(@plot, fitCell);
    %     f1.LineStyle = '--';
    %     f1.LineWidth = 1.5;
    %     f1.Color = [0.2 0.2 0.2];
    % end
    % legend('off');
    % hold off
    %_________________________________________________________________________________
end
end

%   Plot mean of all nuclei
%_______________________________________________________________________________
% function plotAllBleach(hsValMeanNormTMean, hsValTNormSem, timeCell)
function plotAllBleach(dataMean, dataSem, xCell, flag)

if strcmp(flag, 'rich')
    yLab = 'Enrighment zone intensity (a. u)';
elseif strcmp(flag, 'nuc')
    yLab = 'Nuclear intensity (a. u)';
end
yy = smooth(dataMean);
timeAx = vertcat(xCell{:,1});
figure('Color' ,'w');
curve1 = yy + dataSem;
curve2 = yy - dataSem;
xx = [timeAx; flipud(timeAx)];
inBetween = [curve1; flipud(curve2)];
p1 = fill(xx, inBetween,'g');
p1.FaceColor = [0.8 0.8 0.8];
p1.EdgeColor = 'none';
hold on;
l1 = plot(timeAx, yy);
l1.Color = [0.1 0.1 0.1];
l1.Marker = '.';
l1.LineWidth = 1;
l1.LineStyle = 'none';
hold on;
ff = fit(timeAx , yy, 'exp1');
f1 = plot(ff);
f1.LineWidth = 2;
f1.LineStyle = '--';
f1.Color = [0.8 0.2 0.8];
legend('off');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
grid ('on');
xlabel('Time (s)');
ylabel(yLab)

%_______________________________________________________________________________
end

