function loadParams(filename)
% read in parameter file

global brushing
fid = fopen(filename, 'r');
while ~feof(fid)
    tline = fgets(fid);
    eval(tline)
end
fclose(fid);

end