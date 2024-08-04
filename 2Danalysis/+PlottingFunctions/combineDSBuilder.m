function combineDSBuilder(txtFilePath, type)
DS1 = PlottingFunctions.spotPropPlot4Combine(txtFilePath, type);
DS2 = PlottingFunctions.spotPropPlot7Combine(txtFilePath, type);
DS = cell2struct([struct2cell(DS1);struct2cell(DS2)],[fieldnames(DS1);fieldnames(DS2)]);
DS.type = type;
outDSName = (['inOutDS', type]);
procFolder = 'G:\Tyrone_analysis\bcd2x\2D_2xA\combine\combineDS';
save([procFolder, filesep, outDSName, '.mat'],'DS');
end