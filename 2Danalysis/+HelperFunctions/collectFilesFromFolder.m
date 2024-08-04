function collectFilesFromFolder(dataFolderPath, textFile)
% collects filenames from within a folder and writes them sequentially in a
% text file. 
% dataFolderPath: folder containing all datastructure folders
% textFile: text with all the DS folder names. 
dirInfo = dir(dataFolderPath);
dirInfo = dirInfo(~cellfun('isempty', {dirInfo.date})); 
dirInfo = dirInfo(~ismember({dirInfo.name},{'.','..'}));
fileID = fopen(textFile,'w');
fullFilePaths = fullfile(dataFolderPath, {dirInfo.name});
formatSpec = 'file_%03d\t%s\n'; % file number tag (file_001 etc)
for i=1:length(fullFilePaths)
    fprintf(fileID, formatSpec,i,fullFilePaths{i});
end
end