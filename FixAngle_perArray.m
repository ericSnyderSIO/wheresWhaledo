%Fix angle for all

% TO DEFINE: Set up data and paths
OrDetPath = 'D:\';
saveDetPath = 'D:\Corrected angle'; % set path to save tracking detections

FileMask = {'*detections*'};
SearchRecursiv = 0;

[PathFileListMat, FileListMat, PathListMat] = ...
    utFindFiles(FileMask, OrDetPath, SearchRecursiv);

for n = 1:length(PathFileListMat)
    load(PathFileListMat{n});
    
    DET = fixAngle(DET);
    
    detFile = split(FileListMat{n},'.');
    detFileName = [detFile{1},'_corrAngle.mat'];
     save(fullfile(saveDetPath,detFileName),'DET')
    fprintf('File %d/%d saved: %s\n',n,height(n),detFileName)
end