function [fim] = nonlocalFilter2D(im, sx, sy)
J = imnlmfilt(im, 'DegreeOfSmoothing',sx);
fim = J;
end