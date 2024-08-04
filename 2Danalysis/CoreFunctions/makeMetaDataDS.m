function [expInfoDS] = makeMetaDataDS(paramTextFolder)
% makeExpInfo Stores all relevant information/data that will be used
% in the following analyses
% correct folder
expInfoFileName = append(paramTextFolder, filesep,'expInfo.txt');
expInfoFileChar = convertStringsToChars(expInfoFileName);
%where is the  text file containing the relevant parameters (should be in
%the analysis folder)
if exist(expInfoFileName, 'file')
    info = HelperFunctions.readtext(expInfoFileChar,'\t');
    userInfo = cell2struct(info(:,2),info(:,1),1);
else
    error('User must supply expInfo.txt in the analysis directory')
end
%where is the raw file stored (takes in stack of images)
if isfield(userInfo,'rawDataFilePath')
    expInfoDS.expInfo.rawDataFilePath = userInfo.rawDataFilePath;
    [rawFilePath,rawFileName,rawFileExt] = fileparts(userInfo.rawDataFilePath);
%     expInfoDS.expInfo.analysisFolder = analysisFolder;
%     expInfoDS.expInfo.rawImagesFolderName = strcat(analysisFolder, filesep, 'out');
%     expInfoDS.expInfo.procImagesFolderName = strcat(rawFilePath, filesep, 'DS_', rawFileName);
else
    error('expInfo.txt should include rawDataFilePath')
end

if isfield(userInfo,'DSFilePath')
    expInfoDS.expInfo.analysisFolder = userInfo.DSFilePath;
    expInfoDS.expInfo.rawImagesFolderName = strcat(expInfoDS.expInfo.analysisFolder, filesep, 'out');
    expInfoDS.expInfo.procImagesFolderName = strcat(expInfoDS.expInfo.analysisFolder, filesep, 'DS_', rawFileName);
    else
    error('expInfo.txt should include valid DSFilePath')
end

%what is the bit depth of the acquisition
if isfield(userInfo,'bitDepth')
    bitDepth = userInfo.bitDepth;
else
    error('image bit depth not found')
end
% what is the time resolution
if isfield(userInfo,'timeResolution')
    DeltaT = userInfo.timeResolution;
else
    error('image acquisition time not found')
end
% psf in the z direction
if isfield(userInfo,'Zpsf')
    Zpsf = userInfo.Zpsf;
else
    Zpsf = 500;
    warning('no feature size reference found, using 500 nanometers')
end
%psf alonf the horizontal plane
if isfield(userInfo,'XYpsf')
    XYpsf = userInfo.XYpsf;
else
    XYpsf = 200;
    warning('no feature size reference found, using 200 nanometers')
end
%how many z sectons are used
if isfield(userInfo,'Zslices')
    Zslices = userInfo.Zslices;
else
    Zslices = 1;
end
% Flag to determine which image to use for all analyses
if isfield(userInfo,'imageToUse')
    imageToUse = userInfo.imageToUse;
else
    warning('no use image type specified, using raw images for all analyses')
end
switch imageToUse 
    case 'raw'
        imUseType = 1;
    case 'smooth'
        imUseType = 2;
    case 'sharp'
        imUseType = 3;
    otherwise
        imUseType = 1;
end
% if the nuclei touching the borders need to be cleared
if isfield(userInfo,'imClearBorder')
    imClearBorder = userInfo.imClearBorder;
else
    warning('no valid input smoothing filter type found, not clearing border')
end
switch imClearBorder 
    case 'yes'
        clearBorder = 1;
    case 'no'
        clearBorder = 0;
    otherwise
        clearBorder = 0;
end
% what kind of smoothing filter should be used on the iamges
if isfield(userInfo,'smoothFilterType')
    smoothFilterType = userInfo.smoothFilterType;
else
    warning('no valid input smoothing filter type found, using wiener')
end
switch smoothFilterType 
    case 'gaussian'
        nucSmoothFilter = 1;
    case 'wiener'
        nucSmoothFilter = 2;    
    case 'mean'
        nucSmoothFilter = 3;
    case 'median'
        nucSmoothFilter = 4;
    case 'bilateral'
        nucSmoothFilter = 5;
    case 'nonlocal'
        nucSmoothFilter = 6;
    case 'dog'
        nucSmoothFilter = 7;    %Don't use
    otherwise
        nucSmoothFilter = 5;
end
% what kind of sharpening filter should be used on the iamges
if isfield(userInfo,'deconvolutionFilterType')
    deconvolutionFilterType = userInfo.deconvolutionFilterType;
else
    warning('no valid input deconvolution filter found, using wiener')
end
switch deconvolutionFilterType 
    case 'lucy'
        deconvolutionType = 1;
    case 'wiener'
        deconvolutionType = 2;    
    case 'regular'
        deconvolutionType = 3;    
    otherwise
        deconvolutionType = 2;
end
%how to detect spots (most common would be deviation from mean). devation
%from max is better paired with "imClearBorder" set to yes.
if isfield(userInfo,'spotDetectType')
    spotDetectType = userInfo.spotDetectType;
else
    warning('no valid input spot filter type found, deviationFromMax')
end
switch spotDetectType 
    case 'deviationFromMean'
        spotDetect = 1;
    case 'deviationFromMax'
        spotDetect = 2;    
    case 'deviationOfGradient'
        spotDetect = 3;
    otherwise
        spotDetect = 2;
end
%what kind of filtering to use in the spot dete analysis. (most common
%would be top-hat
if isfield(userInfo,'spotFilterType')
    spotFilterType = userInfo.spotFilterType;
else
    warning('no valid input spot filter type found, using top hat')
end
switch spotFilterType 
    case 'gaussian'
        spotFilter = 1;
    case 'bilateral'
       spotFilter = 2;    
    case 'tophat'
        spotFilter = 3;
    case 'median'
        spotFilter = 4;
    otherwise
        spotFilter = 4;
end
%options for detecting nucleus (most common would be watershed)
if isfield(userInfo,'nucleusDetectMethod')
    nucleusDetectMethod = userInfo.nucleusDetectMethod;
else
    nucleusDetectMethod = 'morphclose';
    warning('no nucleus detection type specified, using morphclose')
end
switch nucleusDetectMethod
    case 'reconstruct'
        nucleusDetectMode = 1;
    case 'morphclose'        
        nucleusDetectMode = 2;
    case 'morpherode'
        nucleusDetectMode = 3;
    case 'watershed'
        nucleusDetectMode = 4;
    otherwise
        nucleusDetectMode = 1;
end
% parameter for image thresholding for nucleus detection
if isfield(userInfo,'nucleusIntensityThreshold')
    nucThres = userInfo.nucleusIntensityThreshold;
else
    nucThres = 0.85;
    warning('no threshold value found, using 0.85')
end
% parameter for image smoothing
if isfield(userInfo,'smoothingParam')
    smoothingParam = userInfo.smoothingParam;
else
    smoothingParam = 5;
    warning('no smoothing parameter found, using 5')
end
% cutoff for detected spots in nm (0 = Auto)
if isfield(userInfo,'minSpotSize')
    minSpotSize = userInfo.minSpotSize;
else
    minSpotSize = 0;
    warning('no sharp parameter found, using auto mode = 0')
end
% desired size of the structural elements for image morphological
% transformation
if isfield(userInfo,'elementSize')
    elementSize = userInfo.elementSize;
else
    elementSize = 3;
    warning('no structural element size found, using 3')
end
% size of boundary pad on the edges to be blackened out
if isfield(userInfo,'padSize')
    padSize = userInfo.padSize;
    fprintf('padding of size %d is being used\n', padSize)
else
    padSize = 0;    
end
% how many time points to analyze
if isfield(userInfo,'overrideTimePoints')
    overrideTimePoints = userInfo.overrideTimePoints;
else
    overrideTimePoints = 0;
    warning('using all time points from the file')
end
%starting time frame of analysis
if isfield(userInfo,'startTimePoint')
    startingFrame = userInfo.startTimePoint;
else
    startingFrame = 1;
    warning('using starting frame = 1')
end

% if isfield(userInfo,'startingFrame')
%     startingFrame = userInfo.startingFrame;
% else
%     startingFrame = 1;
%     warning('using starting frame = 1')
% end


%rough cutoff for nucleus diameter in [um]
if isfield(userInfo,'nucleusFeatureSize')
    nucleusFeatureSize = userInfo.nucleusFeatureSize;
else
    nucleusFeatureSize = 3;
    warning('no feature size reference found, using 3 microns')
end
imDS = imfinfo(expInfoDS.expInfo.rawDataFilePath);
T = (length(imDS));
X = (imDS(1).Width);
Y = (imDS(1).Height);
expInfoDS.imagingInfo = struct('X',X,'Y',Y,'T',T,'Zslices',Zslices,'DeltaT',DeltaT);
expInfoDS.imagingInfo.XYpsf = (XYpsf);
expInfoDS.imagingInfo.Zpsf = (Zpsf);
expInfoDS.imagingInfo.bitDepth = (bitDepth);
expInfoDS.imagingInfo.imUseType = imUseType;
expInfoDS.imagingInfo.nucSmoothFilter = nucSmoothFilter;
expInfoDS.imagingInfo.deconvolutionType = deconvolutionType;
expInfoDS.imagingInfo.nucleusDetectMode = nucleusDetectMode;
expInfoDS.imagingInfo.nucThres = nucThres;
expInfoDS.imagingInfo.smoothingParam = smoothingParam;
expInfoDS.imagingInfo.minSpotSize = minSpotSize;
expInfoDS.imagingInfo.clearBorder = clearBorder;
expInfoDS.imagingInfo.elementSize = elementSize;
expInfoDS.imagingInfo.padSize = padSize;
expInfoDS.imagingInfo.nucleusFeatureSize = nucleusFeatureSize;
expInfoDS.imagingInfo.overrideTimePoints = (overrideTimePoints);
expInfoDS.imagingInfo.startingFrame = (startingFrame);
expInfoDS.imagingInfo.spotDetect = spotDetect;
expInfoDS.imagingInfo.spotFilter = spotFilter;
expInfoDS.imagingInfo.timeResolution = DeltaT;
if isfield(userInfo,'XpixelSize')&&isfield(userInfo,'YpixelSize')&&isfield(userInfo,'ZpixelSize')
    expInfoDS.imagingInfo.XpixelSize  = userInfo.XpixelSize;
    expInfoDS.imagingInfo.YpixelSize  = userInfo.YpixelSize;
    expInfoDS.imagingInfo.ZpixelSize  = userInfo.ZpixelSize;
else
    error('expInfo.txt should include pixel size [nm] Xpixel, Ypixel and Zpixel')
end
% by how much the centroid shift of each nucleus between consecutive frames
% be allowed (in pixels)
if isfield(userInfo,'centroidShiftConstraint')
    centroidShiftConstraint = userInfo.centroidShiftConstraint;
else
    centroidShiftConstraint = ceil(nucleusFeatureSize/(2*userInfo.XpixelSize));
    warning('no centroid shift constraint found, using calculated = %d',centroidShiftConstraint)
end
expInfoDS.imagingInfo.centroidShiftConstraint = centroidShiftConstraint;

% total time for which the hotspots are to be analyzed (in seconds)
if isfield(userInfo,'analysisTime')
    analysisTime = userInfo.analysisTime;
else
    analysisTime = floor(overrideTimePoints/DeltaT);    
    warning('no hotspot analysis time span found, using %d seconds', analysisTime);
end

% total time for which the hotspots are to be analyzed (in seconds)
if isfield(userInfo,'timeWindow')
    timeWindow = userInfo.timeWindow;
else
    timeWindow = floor(analysisTime/3);    
    warning('no moving time window for hotspot analysis found, using %d seconds', timeWindow);
end

% slide side of moving time window for hotspot analysis (in seconds)
if isfield(userInfo,'timeWindow')
    timeSlide = userInfo.timeSlide;
else
    timeSlide = '1';    
    warning('no time slide for hotspot analysis found, using %d seconds', timeSlide);
end
expInfoDS.analysisInfo.analysisTime = analysisTime;
expInfoDS.analysisInfo.timeWindow = timeWindow;
expInfoDS.analysisInfo.timeSlide = timeSlide;
end

