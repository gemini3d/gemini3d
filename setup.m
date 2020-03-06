function setup()
%% run this before running Gemini Matlab scripts

cwd = fileparts(mfilename('fullpath'));

for p = {'matlab', 'matlab/vis', 'matlab/vis/plotfunctions', 'matlab/setup', 'matlab/setup/gridgen'}
  addpath([cwd, filesep, p{:}])
end

end % function