function spotPropPlot2(DSfolder)
%______________________________________________________________________________________
%
%	Plots the change of intensity within each hotspot with time
%	time, and mean interval between between consequtive occupancies.
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
end