function setup()
%% configure paths to work with Gemini Matlab

matgemini_git = "https://github.com/gemini3d/mat_gemini.git";
matgemini_dir_name = "mat_gemini";

cwd = fileparts(mfilename('fullpath'));
gemini_matlab = fullfile(cwd, "..", matgemini_dir_name);
if ~isfolder(gemini_matlab)
  cmd = "git -C " + fullfile(cwd, "..") + " clone --recursive " + matgemini_git + " " + matgemini_dir_name;
  ret = system(cmd);
  assert(ret==0, "problem downloading Gemini Matlab functions")
end
run(fullfile(gemini_matlab, "setup.m"))

end
