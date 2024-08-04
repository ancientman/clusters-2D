function saveSomeStuff(ProcFolder, ...
    registerNucPropStruct, ...
    globalRegSpotStruct, ...
    globalHotSpotStruct, ...
    globalSpotHotSpotStruct, ...
    imAddDS, ...
    hotSpotPropNL, ...
    crossCorrStruct, ...
    pairCorrStruct, ...
    conditionDS)
save([ProcFolder,'/nucDS.mat'],'registerNucPropStruct');
save([ProcFolder,'/spotDS.mat'],'globalRegSpotStruct');
save([ProcFolder,'/hotSpotDS.mat'],'globalHotSpotStruct');
save([ProcFolder,'/globalHotSpotDS.mat'],'globalSpotHotSpotStruct');
save([ProcFolder,'/imAddDS.mat'],'imAddDS');
save([ProcFolder,'/hotSpotNLDS.mat'],'hotSpotPropNL'); 
save([ProcFolder,'/xCorrDS.mat'],'crossCorrStruct'); 
save([ProcFolder,'/pairCorrDS.mat'],'pairCorrStruct'); 
save([ProcFolder,'/procConditionsDS.mat'],'conditionDS');
end