function save_glowframe(flagoutput, filename, saveplot_fmt, hf)
%% CREATES IMAGE FILES FROM PLOTS
narginchk(4,4)
validateattributes(flagoutput, {'numeric'}, {'scalar'}, mfilename, 'output flag', 1)
validateattributes(filename, {'char'}, {'vector'}, mfilename, 'aurora filename', 2)
validateattributes(saveplot_fmt, {'cell'}, {'vector'}, mfilename, 'format to save', 3)

[outdir, outname] = fileparts(filename);
outdir = [outdir, '/../Aurplots'];

makedir(outdir)

for i=1:length(saveplot_fmt)

  switch saveplot_fmt{i}
    case 'png', flag = '-dpng'; suffix = '.png';
    case 'eps', flag = '-depsc2'; suffix = '.eps';
    otherwise, continue
  end

  outfile = [outdir, '/', outname, suffix];
  if flagoutput~=3
    disp(['writing ',outfile])
    print(hf,flag,outfile,'-r150')
  end
end

end % function
