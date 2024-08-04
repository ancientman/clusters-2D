function [corrStruct] = autoCorr(imUse, registerNucPropStruct, globalHotSpotStruct, registeredFrames, breakingTimeFrame, metadataDS, plotOn)
im = imUse;
hsDS = globalHotSpotStruct;
nucDS = registerNucPropStruct;
frames = registeredFrames;
tFinal = breakingTimeFrame;
delT = metadataDS.imagingInfo.timeResolution; %in seconds

%~~~~~~~~~~~~~~~~~~~
inParam.wholeFrame = 'no'; %   'Options: 'yes', 'no'. Yes is for whole frame
inParam.wholeNuc = 'no'; % 'Options: 'yes', 'no'. Yes is for whole nuc vs background. use only imclearborder.
inParam.win = 12; % Hardcoded (about one micron): Actual win size is 2*win + 1
inParam.hotSpotBB = 'on'; % 'on' for hotspot bb, 'off' for fixes win around centroid.
inParam.filterType = 'detrend'; %   Oprions 'none', 'medfilt', 'detrend'
%~~~~~~~~~~~~~~~~~~~

nucCent = round(hsDS.globalNucCentroidNL);
corrCoeffFrame = zeros(numel(frames),1);
meanFrame = zeros(numel(frames),1);
%________________________________________________________________________________
if strcmp(inParam.wholeFrame, 'yes')    
    t = 1;
    while t <=numel(frames)
        ff = frames(t);
        meanFrame(t) = mean(im(:,:,t),'all');
        corrCoeffFrame(t) = corr2(im(:,:,1), im(:,:,ff)./(meanFrame(t)/meanFrame(1)));
        t = t+1;
    end
    timeFrame = frames'.*delT;
    fitFrame = fit(timeFrame, meanFrame, 'exp1');  
    corrStruct.frameMean = meanFrame;
    corrStruct.frameCorr = corrCoeffFrame;
    corrStruct.timeFrame = timeFrame;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (strcmp(plotOn,'yes'))
        figure('Color', 'w')
        plot(timeFrame, corrCoeffFrame, 'k-', 'Linewidth', 1.5);
        xlabel('Time (s)');
        ylabel('{\sigma}');
        title('whole frame correlation');
        figure('Color', 'w');
        plot(fitFrame, timeFrame, meanFrame);
        xlabel('Time (s)');
        ylabel('Fluorescence intensity (a. u)');
        title('whole frame bleaching');
    end
else
    %_______________________________________________________________________________
    if strcmp(inParam.wholeNuc, 'yes')
        labelMatGlobal = imclearborder(globalHotSpotStruct.globalNucLabel);
        labelMatGlobal = bwlabel(imbinarize(labelMatGlobal));        
        bgMask = zeros(size(labelMatGlobal));
        [~, bgMask] = Nucleus.removeNuc(im(:,:,1), 1,bgMask);
        
        nNuc = numel(unique(labelMatGlobal)) - 1;
        meanNucWhole = cell(1, (nNuc));
        corrCoeffNucWhole = cell(1, (nNuc));        
        mockIm1 = zeros(size(im, 2), size(im, 1));
        for i=1:nNuc
            t = 1;
            imWin = cell(1, length(frames));
            
            while t <=length(frames)
                ff = frames(t);    
                labelMat = (nucDS{ff}.labelMat == i);
                imWin{ff} = im(:,:,ff);
                imWin{ff}(~labelMat) = 0;
%                 imWin{t} = im(sub2ind(size(im), YNucCent, XNucCent, repmat(ff, size(XNucCent))));
                mean1 = mean(imWin{1},'all');
                meanT = mean(imWin{ff},'all');
                meanNucWhole{i} = vertcat(meanNucWhole{i}, mean(nonzeros(imWin{ff}),'all'));

                switch inParam.filterType
                    case 'none'
                        corrCoeffNucWhole{i} = vertcat(corrCoeffNucWhole{i}, corr2(imWin{1}, imWin{ff}));
                    case 'medfilt'
                        corrCoeffNucWhole{i} = vertcat(corrCoeffNucWhole{i}, corr2(medfilt2(imWin{1}), medfilt2(imWin{ff})));
                    case 'detrend'
                        corrCoeffNucWhole{i} = vertcat(corrCoeffNucWhole{i}, corr2(imWin{1}, imWin{ff}./(meanT/mean1)));
                    otherwise
                        corrCoeffNucWhole{i} = vertcat(corrCoeffNucWhole{i}, corr2(imWin{1}, imWin{ff}));
                end
                t = t+1;
            end
        end  
        
        meanNucWhole = meanNucWhole(~cellfun('isempty',meanNucWhole));
        valNucWholeMean = mean(horzcat(meanNucWhole{:}), 2, 'omitnan');
        valNucWholeSem = std(horzcat(meanNucWhole{:}), 0, 2, 'omitnan')./sqrt(nNuc);     
        
        corrCoeffNucWhole = corrCoeffNucWhole(~cellfun('isempty',corrCoeffNucWhole));
        corrNucWholeMean = mean(horzcat(corrCoeffNucWhole{:}), 2, 'omitnan');
        corrNucWholeSem = std(horzcat(corrCoeffNucWhole{:}), 0, 2, 'omitnan')./sqrt(nNuc);

        timeNucWhole = repmat({frames'.*delT}, [1, length(meanNucWhole)]);
        fitNucWhole = cellfun(@(x, y) fit(x , y, 'exp1'), timeNucWhole, meanNucWhole, 'un', 0);  
        
        %-------------------------------------------------------------------------
        
        t = 1;
        bgWin = cell(1, length(frames));
        meanBg = 0;
        corrCoeffBg = 0;
        while t <=length(frames)
            ff = frames(t);   
            bgWin{ff} = im(:,:,ff);
            bgWin{ff}(~bgMask) = 0;

            mean1 = mean(nonzeros(bgWin{1}),'all');
            meanT = mean(nonzeros(bgWin{ff}),'all');
            meanBg = vertcat(meanBg, meanT);

            switch inParam.filterType
                case 'none'
                    corrCoeffBg = vertcat(corrCoeffBg, corr2(bgWin{1}, bgWin{ff}));
                case 'medfilt'
                    corrCoeffBg = vertcat(corrCoeffBg, corr2(medfilt2(bgWin{1}), medfilt2(bgWin{ff})));
                case 'detrend'
                    corrCoeffBg = vertcat(corrCoeffBg, corr2(bgWin{1}, bgWin{ff}./(meanT/mean1)));
                otherwise
                    corrCoeffBg = vertcat(corrCoeffBg, corr2(bgWin{1}, bgWin{ff}));
            end
            t = t+1;
        end
        
        valBgMean = meanBg(2:end);
        valBgSem = zeros(size(valBgMean));
        timeBg = zeros(size(valBgMean));
        
        corrBgMean = corrCoeffBg(2:end);
        corrBgSem = zeros(size(corrBgMean));
        
        %-------------------------------------------------------------------------
        
        corrStruct.nucWholeMean = meanNucWhole;
        corrStruct.nucWholeCorr = corrCoeffNucWhole;
        corrStruct.nucWholeTime = timeNucWhole;
        
        corrStruct.bgMean = meanBg;
        corrStruct.bgCorr = corrCoeffBg;
        corrStruct.bgTime = timeBg;      
        
        if (strcmp(plotOn, 'yes'))
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            plot1(frames, delT, valNucWholeMean, valNucWholeSem, valBgMean, valBgSem, 'intensity', inParam);

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            plot1(frames, delT, corrNucWholeMean, corrNucWholeSem, corrBgMean, corrBgSem, 'corr', inParam);
        end
        %___________________________________________________________________________
    else
        
        meanNucCent = cell(1, size(nucCent, 1));
        corrCoeffNucCent = cell(1, size(nucCent, 1));
        mockIm1 = zeros(size(im, 2), size(im, 1));
        win = inParam.win;

        if strcmp(inParam.hotSpotBB,'on')
            hsProp = regionprops(logical(globalHotSpotStruct.hotspotBW), 'BoundingBox');
            hsProp = num2cell(hsProp)';
            bbHs = cellfun(@(x) round(x.BoundingBox), hsProp, 'un', 0);
        %     bbMeanArea = mean(cellfun(@(x) (x(3)*x(4)), bbHs));
            nucwin = round(sqrt(mean(cellfun(@(x) (x(3)*x(4)), bbHs))));
            win = nucwin;
        end

        for i = 1:size(nucCent, 1)
            xNucCent = (nucCent(i, 1) - win) : (nucCent(i, 1) + win);
            yNucCent = (nucCent(i, 2) - win) : (nucCent(i, 2) + win);
            [XNucCent,YNucCent] = meshgrid(xNucCent, yNucCent);
            if  ~isempty(setdiff(xNucCent,1:size(im, 1))) || ~isempty(setdiff(yNucCent,1:size(im, 2)))
                continue
            else  
                t = 1;
                mockIm1(sub2ind(size(mockIm1), YNucCent, XNucCent)) = 1;
                imWin = cell(1, numel(frames));
                while t <=numel(frames)
                    ff = frames(t);
                    imWin{t} = im(sub2ind(size(im), YNucCent, XNucCent, repmat(ff, size(XNucCent))));
                    mean1 = mean(imWin{1},'all');
                    meanT = mean(imWin{t},'all');
                    meanNucCent{i} = vertcat(meanNucCent{i}, meanT);

                    switch inParam.filterType
                        case 'none'
                            corrCoeffNucCent{i} = vertcat(corrCoeffNucCent{i}, corr2(imWin{1}, imWin{t}));
                        case 'medfilt'
                            corrCoeffNucCent{i} = vertcat(corrCoeffNucCent{i}, corr2(medfilt2(imWin{1}), medfilt2(imWin{t})));
                        case 'detrend'
                            corrCoeffNucCent{i} = vertcat(corrCoeffNucCent{i}, corr2(imWin{1}, imWin{t}./(meanT/mean1)));
                        otherwise
                            corrCoeffNucCent{i} = vertcat(corrCoeffNucCent{i}, corr2(imWin{1}, imWin{t}));
                    end
                    t = t+1;
                end
            end
        end

        hsCent = round(vertcat(hsDS.hotspotCentroidUniq{:}));
        meanHs = cell(1, size(hsCent, 1));
        corrCoeffHs = cell(1, size(hsCent, 1));

        mockIm2 = zeros(size(im, 2), size(im, 1));

        for i = 1:size(hsCent, 1)

            if strcmp(inParam.hotSpotBB, 'on')    
                %   Use bounding box of each global hotspot
                xHs = bbHs{i}(1)  : bbHs{i}(1) + bbHs{i}(3);
                yHs = bbHs{i}(2)  : bbHs{i}(2) + bbHs{i}(4);
                [XHs,YHs] = meshgrid(xHs, yHs);        
            else
                xHs = (hsCent(i, 1) - win) : (hsCent(i, 1) + win);
                yHs = (hsCent(i, 2) - win) : (hsCent(i, 2) + win);
                [XHs,YHs] = meshgrid(xHs, yHs);
            end       

            if  ~isempty(setdiff(xHs,1:size(im, 1))) || ~isempty(setdiff(yHs,1:size(im, 2)))
                continue
            else    
                mockIm2(sub2ind(size(mockIm2), YHs, XHs)) = 1;    
                t = 1;
                imWin = cell(1, numel(frames));
                while t <=numel(frames)
                    ff = frames(t);        
                    imWin{t} = im(sub2ind(size(im), YHs, XHs, repmat(ff, size(XHs))));
                    mean1 = mean(imWin{1},'all');
                    meanT = mean(imWin{t},'all');
                    meanHs{i} = vertcat(meanHs{i}, mean(imWin{t},'all'));

                    switch inParam.filterType
                        case 'none'
                            corrCoeffHs{i} = vertcat(corrCoeffHs{i}, corr2(imWin{1}, imWin{t}));
                        case 'medfilt'
                            corrCoeffHs{i} = vertcat(corrCoeffHs{i}, corr2(medfilt2(imWin{1}), medfilt2(imWin{t})));
                        case 'detrend'
                            corrCoeffHs{i} = vertcat(corrCoeffHs{i}, corr2(imWin{1}, imWin{t}./(meanT/mean1)));
                        otherwise
                            corrCoeffHs{i} = vertcat(corrCoeffHs{i}, corr2(imWin{1}, imWin{t}));
                    end
                    t = t+1;
                end
           end
        end
        %________________________________________________________________________________
        %
        %   Mean values of the areas under consideration
        %________________________________________________________________________________

        meanNucCent = meanNucCent(~cellfun('isempty',meanNucCent));
        meanHs = meanHs(~cellfun('isempty',meanHs));
        valHsMean = mean(horzcat(meanHs{:}), 2, 'omitnan');
        valHsSem = std(horzcat(meanHs{:}), 0, 2, 'omitnan')./sqrt(size(hsCent, 1));
        valNucMean = mean(horzcat(meanNucCent{:}), 2, 'omitnan');
        valNucSem = std(horzcat(meanNucCent{:}), 0, 2, 'omitnan')./sqrt(size(nucCent, 1));

        timeNuc = repmat({frames'.*delT}, [1, length(meanNucCent)]);
        fitNuc = cellfun(@(x, y) fit(x , y, 'exp1'), timeNuc, meanNucCent, 'un', 0);  

        timeHs = repmat({frames'.*delT}, [1, length(meanHs)]);
        fitHs = cellfun(@(x, y) fit(x , y, 'exp1'), timeHs, meanHs, 'un', 0);  
        
        corrHsMean = mean(horzcat(corrCoeffHs{:}), 2, 'omitnan');
        corrHsSem = std(horzcat(corrCoeffHs{:}), 0, 2, 'omitnan')./sqrt(size(hsCent, 1));
        corrNucMean = mean(horzcat(corrCoeffNucCent{:}), 2, 'omitnan');
        corrNucSem = std(horzcat(corrCoeffNucCent{:}), 0, 2, 'omitnan')./sqrt(size(nucCent, 1));
        
         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        corrStruct.nucCentMean = meanNucCent;
        corrStruct.nucCentCorr = meanNucCent;
        corrStruct.nucCentTime = meanNucCent;
        
        corrStruct.hsMean = meanHs;
        corrStruct.hsCorr = corrCoeffHs;
        corrStruct.hsTime = timeHs;          
        
        if (strcmp(plotOn, 'yes'))
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            plot1(frames, delT, valHsMean, valHsSem, valNucMean, valNucSem, 'intensity', inParam);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            plot1(frames, delT, corrHsMean, corrHsSem, corrNucMean, corrNucSem, 'corr', inParam);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
        
        figure('Color', 'w')
        imshow(bwperim(globalHotSpotStruct.globalNucLabel) | ...
        globalHotSpotStruct.hotspotBW  | ...
        bwperim(mockIm1) | bwperim(mockIm2)); 
        hold on; 
        plot(hsCent(:,1), hsCent(:,2), '*r');
        hold off;
    end
    %_________________________________________________________________________________
end
end

function plot1(frames, delT, mean1, sem1, mean2, sem2, tag, param)
figure('Color', 'w')
xx = frames'.*delT;
ys = sem1;
ym = smooth(mean1);
curve1 = ym + ys;
curve2 = ym - ys;
xx = [xx; flipud(xx)];
inBetween = [curve1; flipud(curve2)];
d1 = fill(xx, inBetween,'g');
d1.FaceColor = [0.5, 0.8, 0.9];
d1.FaceAlpha = 0.5;
d1.EdgeColor = 'none';
hold on;
%----------------------------------------------
p1 = plot(frames'.*delT, smooth(mean1));
p1.LineWidth = 2;
p1.Color = [0.5, 0.8, 0.9];
hold on;
%-----------------------------------------------
xx = frames'.*delT;
ys = sem2;
ym = smooth(mean2);
curve1 = ym + ys;
curve2 = ym - ys;
xx = [xx; flipud(xx)];
inBetween = [curve1; flipud(curve2)];
d2 = fill(xx, inBetween,'g');
d2.FaceColor = [0.9, 0.5, 0.7];
d2.FaceAlpha = 0.5;
d2.EdgeColor = 'none';
hold on;
%-----------------------------------------------
%-----------------------------------------------
p2 = plot(frames'.*delT, smooth(mean2));
p2.LineWidth = 2;
p2.Color = [0.9, 0.5, 0.7];
hold on;
%-----------------------------------------------
legend([p1 p2],{'rich zone','nuc center'})
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
grid ('on');
xlabel('Time (s)');

if (strcmp(tag, 'intensity'))
    if (strcmp(param.wholeNuc, 'yes'))
        title('Mean intensity of nuclei');
        legend([p1 p2],{'nuclei','background'})
    elseif (strcmp(param.wholeNuc, 'no'))
        if strcmp(param.hotSpotBB, 'on')
            title(['Mean intensity within areas','hotspot:bb']);
        else
            title(['Mean intensity within areas',' pixels']);
        end    
        legend([p1 p2],{'rich zone','nuc center'})        
    end
    ylabel('Intensity (a. u.)');
elseif (strcmp(tag, 'corr'))
    if (strcmp(param.wholeNuc, 'yes'))
        title('Aurocorrelation of nuclei pixels');
        legend([p1 p2],{'nuclei','background'});        
    elseif (strcmp(param.wholeNuc, 'no'))
        if strcmp(param.hotSpotBB, 'on')
            title(['Autocorrelation: nuc win ', ' pixels. hotspot: bb']);
        else
            title(['Autocorrelation: win ', ' pixels']);
        end
         legend([p1 p2],{'rich zone','nuc center'})        
    end
        ylabel('correlation Coefficient')
end
hold off;
%_________________________________________________________________________________
end
