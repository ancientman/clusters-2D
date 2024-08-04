function [fim] = meanFilter2D(imConEn, sx, sy)
kernel = 1/(sx*sy)*ones(sx, sy);
fim = conv2(single(imConEn), kernel, 'same');
end