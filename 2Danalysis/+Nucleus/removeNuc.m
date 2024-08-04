function [imNucRemoved, mask] = removeNuc(im, iter, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allows to interactively mask unwanted nuclei one by one
% Roi masking happens only in the first frame
% in all other frames this function is called 
% but the same roi, drawn in
% the first frame is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iter==1
%     polyVert = {};
    prompt = "another nucleus?";
    i = 1;
    mask = zeros(size(im));
    fprintf("\ndraw polygons around nuclei to reject, one by one\n");
    while i>=1
        imshow(im,[]);
        p = drawpolygon('LineWidth',5,'Color','cyan');
        mask = mask + createMask(p, size(mask,2), size(mask, 1));
%         polyVert{i} = p.Position;        
        another = input(prompt);% 1 if another nucleus is needed to be masked, 0 otherwise
        if another == 1
            i = i+1;
        else
            i = 0;
        end
    end
%     J = im;
%     for j=1:length(polyVert)
%         J = regionfill(J,polyVert{j}(:,1), polyVert{j}(:,2));
%     end
end
J = im.*(~mask);
imNucRemoved = J;
end
        
  