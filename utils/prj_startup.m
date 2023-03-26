function prj_startup()
%set_up_project  Configure the environment for this project
%
%   Set up the environment for the current project. This function is set to
%   Run at Startup.

% get the current Matlab version number
prj = currentProject;
% Get the project root folder:
projectRoot = prj.RootFolder;

% if exist(fullfile(projectRoot, 'casadi'), 'dir')
%     addPath(prj, fullfile(projectRoot, 'casadi'));
% end


% Set the location of slprj to be the "work" folder of the current project:
myCacheFolder = fullfile(projectRoot, 'work');
myCacheGenFolder = fullfile(projectRoot, 'build');
% if ~exist(myCacheFolder, 'dir')
%     mkdir(myCacheFolder)
% end
Simulink.fileGenControl('set', 'CacheFolder', myCacheFolder, ...
   'CodeGenFolder', myCacheGenFolder,'createDir', true);