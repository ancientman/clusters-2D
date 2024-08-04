function createOutputFolder(metaDataDS)
rawfname = metaDataDS.expInfo.rawImagesFolderName;
procfname = metaDataDS.expInfo.procImagesFolderName;
overrideTimePoints = (metaDataDS.imagingInfo.overrideTimePoints);
startFrame = (metaDataDS.imagingInfo.startingFrame);
nucSize = metaDataDS.imagingInfo.nucleusFeatureSize;
pixelSize = metaDataDS.imagingInfo.XpixelSize;
padSize = metaDataDS.imagingInfo.padSize;
if metaDataDS.imagingInfo.bitDepth == 8
    bitDepth = 8;
elseif metaDataDS.imagingInfo.bitDepth == 16
    bitDepth = 16;
else
    warning('bit depth neither 8, nor 16, using 16');
    bitDepth = 16;
end  
if ~exist(procfname,'file')
    mkdir(procfname);
else
    dirInfo = dir(procfname);
    dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
    dirInfo([dirInfo.isdir]) = [];% skip subdirectories
    filenames = fullfile(procfname, {dirInfo.name});
    if size(dirInfo, 1)>0
        delete (filenames{:})
    end
end
if ~exist(rawfname,'file')
    mkdir(rawfname);
else
    dirInfo = dir(rawfname);
    dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
    dirInfo([dirInfo.isdir]) = [];% skip subdirectories
    filenames = fullfile(rawfname, {dirInfo.name});
    tic;
    if size(dirInfo, 1)>0
        delete (filenames{:})
    end
    toc;
end
stackInfo = imfinfo(append(metaDataDS.expInfo.rawDataFilePath));
options.append = 0;
options.message = 0;
if startFrame == 1
    if overrideTimePoints
        if overrideTimePoints<length(stackInfo)
            totalFrames = overrideTimePoints;
        else
            totalFrames = length(stackInfo);
        end
    else
        totalFrames = length(stackInfo);
    end
elseif startFrame ~= 1
    fprintf('\nstarting frame is %d\n\n', startFrame);
    if overrideTimePoints == 0
        totalFrames = (length(stackInfo)-startFrame+1);
    elseif overrideTimePoints<=((length(stackInfo)-startFrame+1))...
            &&overrideTimePoints>=1
        totalFrames = overrideTimePoints;
    elseif overrideTimePoints>((length(stackInfo)-startFrame+1))
        totalFrames = (length(stackInfo)-startFrame+1);
        fprintf("\nOverrideframes larger than remaining frames, using %d\n",totalFrames);
    end
end
if padSize>ceil(nucSize/pixelSize)
    warning('padding bigger than on nucleus limit');
end
for i = startFrame:(startFrame + totalFrames - 1)
    frameSeeker = startFrame - 1 + i;
    imFrame = imread(append(metaDataDS.expInfo.rawDataFilePath), frameSeeker);
    if padSize>0% If padding is used blacken out boundary pixels.
        imFrame(1:(padSize+1), :) = 0;% Top boundary
        imFrame(size(imFrame, 1)-(padSize) : size(imFrame, 1), :) = 0;% Bottom boundary
        imFrame(:, 1:(padSize+1)) = 0;% Left boundary
        imFrame(:, size(imFrame, 1)-(padSize) : size(imFrame, 1)) = 0;% Right boundary
    end
    fname = append(metaDataDS.expInfo.rawImagesFolderName,filesep,'file_',num2str(i,'%04.f'), '.tif');
    if bitDepth == 16
        HelperFunctions.saveTiff(uint16(imFrame),fname, options);
    elseif bitDepth == 8
        HelperFunctions.saveTiff(uint8(imFrame),fname, options);
    else
        HelperFunctions.saveTiff(uint16(imFrame),fname, options);
    end       
end
end