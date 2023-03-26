url = 'https://github.com/casadi/casadi/releases/download/3.5.5/casadi-windows-matlabR2016a-v3.5.5.zip'; 
websave('casadi.zip',url);
unzip('casadi.zip','casadi');

% get the current Matlab version number
prj = currentProject;
% Get the project root folder:
projectRoot = prj.RootFolder;

if exist(fullfile(projectRoot, 'casadi'), 'dir')
    addPath(prj, fullfile(projectRoot, filesep, 'casadi'));
end