function nucPropPlot(nucProp, spotProp, tMax)
figure;
for t = 1:tMax
    boundaries = bwboundaries(nucProp{t}.labelMat);
    numberOfBoundaries = size(boundaries, 1); 
    plot(nucProp{t}.centroid(:,1), nucProp{t}.centroid(:,2), 'r*'); 
    hold on;
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
        hold on;
        plot(spotProp{t}.spotCentroid{k}(:,1), spotProp{t}.spotCentroid{k}(:,2), 'black.');
        hold on;
    end
    fprintf('plot time = %d\n', t);
end
xlim([0, size(nucProp{t}.labelMat,2)])
ylim([0, size(nucProp{t}.labelMat,1)])
end
% % first run
% x1 = rand(100,1);
% y1 = rand(100,1);
% x2 = rand(100,1);
% y2 = rand(100,1);
% x3 = rand(100,1);
% y3 = rand(100,1);
% x4 = rand(100,1);
% y4 = rand(100,1);
% fOp.figName='my_title';
% hf = findobj('Type','figure','Name',fOp.figName);
% if ~isempty(hf)
%     % get the handles if the figure already exist
%     hf = figure(hf);
%     h1 = hf.Children(4); % subplot 221
%     h2 = hf.Children(3); % subplot 222
%     h3 = hf.Children(2); % subplot 223
%     h4 = hf.Children(1); % subplot 224
% else
%     % create figure and subplots if the figure does not exist
%     hf = figure('Name',fOp.figName);
%     h1=subplot(2,2,1);
%     h2=subplot(2,2,2);
%     h3=subplot(2,2,3);
%     h4=subplot(2,2,4);
% end
% hold(h1,'on');
% plot(h1,x1,y1);
% hold(h2,'on');
% plot(h2,x2,y2);
% hold(h3,'on');
% plot(h3,x3,y3);
% hold(h4,'on');
% plot(h4,x4,y4);
% % second run
% x1 = rand(100,1);
% y1 = rand(100,1);
% x2 = rand(100,1);
% y2 = rand(100,1);
% x3 = rand(100,1);
% y3 = rand(100,1);
% x4 = rand(100,1);
% y4 = rand(100,1);
% fOp.figName='my_title';
% hf = findobj('Type','figure','Name',fOp.figName);
% if ~isempty(hf)
%     % get the handles if the figure already exist
%     hf = figure(hf);
%     h1 = hf.Children(4); % subplot 221
%     h2 = hf.Children(3); % subplot 222
%     h3 = hf.Children(2); % subplot 223
%     h4 = hf.Children(1); % subplot 224
% else
%     % create figure and subplots if the figure does not exist
%     hf = figure('Name',fOp.figName);
%     h1=subplot(2,2,1);
%     h2=subplot(2,2,2);
%     h3=subplot(2,2,3);
%     h4=subplot(2,2,4);
% end
% hold(h1,'on');
% plot(h1,x1,y1);
% hold(h2,'on');
% plot(h2,x2,y2);
% hold(h3,'on');
% plot(h3,x3,y3);
% hold(h4,'on');
% plot(h4,x4,y4);