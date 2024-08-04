function [fim] = dogFilter2D(im,Sxy)
sxyIn  = Sxy(1);
sxyOut = Sxy(2);
xy = -2*round(sxyOut):2*round(sxyOut);
FxyIn  = normpdf(xy,0,sxyIn);
FxyOut = normpdf(xy,0,sxyOut);
FxyIn  = FxyIn/sum(FxyIn);
FxyOut = FxyOut/sum(FxyOut);
imIn  = imfilter(im,FxyIn','replicate','conv');
imOut = imfilter(im,FxyOut','replicate','conv');
imIn  = imfilter(imIn,FxyIn,'replicate','conv');
imOut = imfilter(imOut,FxyOut,'replicate','conv');
fim = imIn - imOut;
% t1 = fim(:,:,:,1);
% t2 = fim2(:,:,:,1);
% figure(1)
% imagesc(t1(:,:,5))
% figure(2)
% imagesc(t2(:,:,5))
% figure(3)
% dog2 = FxyIn'*FxyIn - FxyOut'*FxyOut;
% imagesc(dog2)
% figure(4)
% [X,Y] = ndgrid(xy,xy);
% Fout = mvnpdf([X(:),Y(:)],zeros(1,2),diag([sxyOut,sxyOut])^2);
% Fin  = mvnpdf([X(:),Y(:)],zeros(1,2),diag([sxyIn,sxyIn])^2);
% Fout = reshape(Fout,numel(xy),[])/sum(Fout(:));
% Fin  = reshape(Fin,numel(xy),[])/sum(Fin(:));
% dog1 = Fin - Fout;
% imagesc(dog1)

end

