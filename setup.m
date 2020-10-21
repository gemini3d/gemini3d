function setup(version)
arguments
  version string = string.empty
end
%% configure paths to work with Gemini Matlab

matgemini_git = "https://github.com/gemini3d/mat_gemini.git";
matgemini_dir_name = "mat_gemini";

cwd = fileparts(mfilename('fullpath'));
gemini_matlab = fullfile(cwd, "..", matgemini_dir_name);
if ~isfolder(gemini_matlab)
  cmd = append("git -C ", fullfile(cwd, ".."), " clone --recursive ", matgemini_git, " ", matgemini_dir_name);

  ret = system(cmd);
  assert(ret==0, "problem downloading Gemini Matlab functions")
end

%% checkout MatGemini version, IF no MatGemini development changes are detected
if ~isempty(version) && strlength(version) > 0
  cmd = append("git -C ", gemini_matlab, " status --porcelain");
  [ret, stat] = system(cmd);
  if ret==0 && isempty(stat)
    cmd = append("git -C ", gemini_matlab, " checkout ", version);
    ret = system(cmd);
    assert(ret==0, "could not checkout MatGemini %s", version)
  end
end
%% enact
run(fullfile(gemini_matlab, "setup.m"))

end
