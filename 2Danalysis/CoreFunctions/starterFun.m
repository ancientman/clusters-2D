function starterFun(analysisFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the starting function to run
% Create experiment information data structure
% run analysis file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaDataDS = makeMetaDataDS(analysisFolder);% constructs datastruct from expinfo.txt with user input
save(append(analysisFolder,filesep,'metaDataDS.mat'),'metaDataDS')
createOutputFolder(metaDataDS);% creates a copy of the input image with split frames
startAnalysis(metaDataDS);% main analysis function
end