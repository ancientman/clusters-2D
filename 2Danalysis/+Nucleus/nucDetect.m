function [fim] = nucDetect(im, thres, mode, el, minNucSize, clearBorder, t, stopMark)
if mode==1 % reconstruct
    maskSeed = Nucleus.seedNuclearMask(im, el, minNucSize, thres, t, stopMark);
elseif mode==2 % morphclose
    maskSeed = Nucleus.seedNuclearMask2(im, el, minNucSize, thres, t, stopMark);
elseif mode==3 % morpherode
    maskSeed = Nucleus.seedNuclearMask3(im, el, minNucSize, thres, t, stopMark);
elseif mode==4 % watershed
    [maskSeed, bgm] = Nucleus.seedNuclearMask4(im, el, minNucSize, thres, t, stopMark);
else % reconstruct
    maskSeed = Nucleus.seedNuclearMask(im, el, minNucSize, thres, t, stopMark);
end
nucBW = maskSeed;
if t~=0
%     [nucBW, thres] = Nucleus.noisyNuclearMask(im, maskSeed, thres, minNucSize, clearBorder, t, stopMark);
end
nucBW2 = bwareaopen(nucBW, minNucSize);
if clearBorder
    nucBW2 = imclearborder(nucBW2);
end
se = strel('disk', el);
Io = imopen(nucBW2, se);
If = imfill(Io, 'holes');
nucBW3 = bwareaopen(If, minNucSize);    
if t == stopMark
    fprintf("in nucDetect");
end
fim = nucBW3;
end