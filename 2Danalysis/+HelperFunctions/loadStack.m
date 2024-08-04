function [ stack ] = loadStack(folder, l)
% Load an entire folder of TIFFs into a single 4D matrix

files = dir([folder,filesep,'*.tif']);
if nargin < 2
    l = length(files);
end

for i = 1:l
    if i == 1
        im1 = loadTiff([folder, filesep, files(1).name]);
        stack = zeros([size(im1), l],'single');
        stack(:,:,:,1) = im1;
    else
        stack(:,:,:,i) = loadTiff([folder, filesep, files(i).name]);
    end
end

end