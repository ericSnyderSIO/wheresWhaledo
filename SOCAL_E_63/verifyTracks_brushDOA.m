filepath = 'M:\Corrected angle';
D = dir([filepath, '\*mod*.mat']);

charAftermodStr = 7; % how many characters after beginning of "mod" is the modification number

nd = 1;
while nd <=numel(D)
    % find all DET files associated with this encounter:
    Imod = strfind(D(nd).name, 'mod'); % index in file name with "mod"
    
    Denc = dir(fullfile(filepath, [D(nd).name(1:Imod-1), 'mod*'])); % all files associated with this encounter

    load(fullfile(filepath, Denc(end).name))
    
    DET = brushDOA(DET{1}, DET{2});
    
    newPath = ['D:\SOCAL_E_63\tracking\interns2022\ericEdits\', Denc(end).name(23:Imod-2)];

    mkdir(newPath)
    
    save(fullfile(newPath, Denc(end).name(1:Imod-2)), 'DET')

    nd = nd+numel(Denc); % skip to next encounter

end