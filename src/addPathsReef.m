function addPathsReef()
    folder = fileparts(which('addPaths'));
%    addpath(folder);
    addpath(genpath(folder));
end