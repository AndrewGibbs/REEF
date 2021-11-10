function addPaths()
    folder = fileparts(which('addPaths'));
%    addpath(folder);
    addpath(genpath(folder));
end