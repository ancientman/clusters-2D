function [deconvIm] = sharpRawImage(im,deconvolutionType,sharpParam)
hSize = [13 13];
sigma = sharpParam;
psf = fspecial('gaussian',hSize,sigma)/sum(ones(hSize),'all');
I = im;
if deconvolutionType==1 % using "lucy"
    I = edgetaper(I,psf);
    J = deconvlucy(I, psf, 20);
elseif deconvolutionType==2 % using "wiener"
    I = edgetaper(I,psf);
    J = deconvreg(I,psf);
elseif deconvolutionType==3 %using "regular"
    I = edgetaper(I,psf);
    J = deconvwnr(I,psf);        
end
JJ = deconvblind({J}, {psf});
deconvIm = uint16(JJ{2});
end