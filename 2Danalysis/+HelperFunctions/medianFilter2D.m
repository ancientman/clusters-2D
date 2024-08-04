function [fim] = medianFilter2D(im,sx, sy)
kernel = [sx, sy];
fim = medfilt2(im, kernel);
end