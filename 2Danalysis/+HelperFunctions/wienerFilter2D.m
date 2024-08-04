function [fim] = wienerFilter2D(I, sx, sy)
fim = wiener2(I, [sx, sy]);
end