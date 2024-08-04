function [fim] = mexicanHatFilter2D(I)
    im = single(I);
    len = 2;
    dist = 2*len + 1;
    % First, enhance contrast by removing the background
    ox = ones(1, (2 * len) + 1); 
    ox = ox / sum(ox(:)); 
    ox2 = ox' * ox;
    im_mean = imfilter(im, ox2);
    im_std = sqrt(imfilter((im - im_mean).^2, ox2));
    goodim = (im - im_mean) ./ (im_std + eps);    
    % In my experience, the filter that works fine is
    % F=fspecial('gaussian',30,1)-fspecial('gaussian',30,6);
    % On my pictures a typical internuclear distance was 12 pixels. 
    % So, scale the parameters appropriately.
    %F=fspecial('gaussian',30,dist/12)-fspecial('gaussian',30,dist/2);
    F = fspecial('gaussian', 30, dist / 15) - fspecial('gaussian', 30, dist / 2);
    filtered = imfilter(goodim, F); 
    % The first pass tries to separate the nuclei that may get lumped together.
    % Then, remove all negative values from between the nuclei and run the filter
    % again, on a cleaner image this time; this should remove the false nuclei
    % -- tiny bright spots too close to the other nuclei.
    %F2 = fspecial('gaussian', 30, dist / 12) - fspecial('gaussian', 30, dist / 3);
    F2=fspecial('gaussian', 30, dist / 15)-fspecial('gaussian', 30, dist / 3);
    fim = filtered; 
    fim(fim < 0) = 0; 
    fim = imfilter(fim, F2);
end