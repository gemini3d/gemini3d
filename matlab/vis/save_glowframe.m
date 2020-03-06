function save_glowframe(flagoutput, filename, saveplot_fmt, hf)
%% CREATES IMAGE FILES FROM PLOTS
narginchk(4,4)
validateattributes(flagoutput, {'numeric'}, {'scalar'}, mfilename, 'output flag', 1)
validateattributes(filename, {'char'}, {'vector'}, mfilename, 'aurora filename', 2)

res = '-r150';

if isempty(saveplot_fmt)
  return
elseif ischar(saveplot_fmt)
  saveplot_fmt = {saveplot_fmt};
end

[outdir, outname] = fileparts(filename);
outdir = [outdir, '/../plots'];

for i=1:length(saveplot_fmt)

  [flag, suffix] = printflag(saveplot_fmt{i});

  outfile = [outdir, '/aurora-', outname, suffix];
  if flagoutput~=3
    disp(['writing ', outfile])
    print(hf,flag,outfile, res)
  end
end

end % function
