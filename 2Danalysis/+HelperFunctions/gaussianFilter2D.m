function [fim] = gaussianFilter2D(im,sx,sy)
x = -2*round(sx):2*round(sx);
y = -2*round(sy):2*round(sy);
%% separable implementation
Fx  = normpdf(x,0,sx);
Fy  = normpdf(y,0,sy);
Fx = Fx/sum(Fx);
Fy = Fy/sum(Fy);
fim = imfilter(im,Fx','replicate','conv');
fim = imfilter(fim,Fy,'replicate','conv');
end

