function prj_shutdown()
%prj_shutdown   Clean up local customizations of the environment
% 
%   Clean up the environment for the current project. This function undoes
%   the settings applied in "set_up_project". It is set to Run at Shutdown.

% Reset the location where generated code and other temporary files are
% created (slprj) to the default:
Simulink.fileGenControl('reset');

% Save and Close all Data Dictionaries
Simulink.data.dictionary.closeAll('-save');
