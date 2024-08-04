function [fim, L] = watershedSpot(im, labelMat, mask, minSpotSize, spotFilter)
% mask = bwareaopen(mask,minSpotSize);
% im = rescale(im);
im(mask==0) = 0;
% el = 3;
% switch spotFilter
%     case 1% gaussian (recomended for global spots)
%         imBlur = imgaussfilt(im);        
%     case 2% bilateral (recomended for global spots)
%         imBlur = imbilatfilt(im);
%     case 3% tophat (not recomended for global spots)
%         se = strel('disk', el, 8);
%         imBlur = imtophat(im,se); 
%         imBlur = imadjust(imBlur);        
%     case 4% median (not recomended for global spots)
%         imBlur = medfilt2(im);        
%     otherwise% none
%         imBlur = im;
% end

imBlur = imgaussfilt(im, 3);   
bw = imBlur;%>0;
bw = bwareaopen(bw,minSpotSize);
D = bwdist(~bw);
D = -D;
D(~bw) = -Inf;
L = watershed(D);
bw2 = mask;
bw2(L==0) = 0;
% bw3 = bw2;
bw3 = bwareaopen(bw2, minSpotSize); % opened
bw3 = imclearborder(bw3);
spotLabel = labelMat.*bw3;
fim = bw3;
L = spotLabel;
end