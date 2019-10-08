function save_glowframe(flagoutput, direc, filename, saveplot_fmt, hf)
%% CREATES IMAGE FILES FROM PLOTS
narginchk(5,5)
validateattributes(flagoutput, {'numeric'}, {'scalar'}, mfilename)
validateattributes(direc, {'char'}, {'vector'}, mfilename)
validateattributes(filename, {'char'}, {'vector'}, mfilename)
validateattributes(saveplot_fmt, {'cell'}, {'vector'}, mfilename)

outdir = [direc, filesep, 'Aurplots'];

if ~is_folder(outdir)
  mkdir(outdir);
end

for i=1:length(saveplot_fmt)

  switch saveplot_fmt{i}
    case 'png', flag = '-dpng'; suffix = '.png';
    case 'eps', flag = '-depsc2'; suffix = '.eps';
    otherwise, continue
  end

  outfile = [outdir,filesep,filename,suffix];
  if flagoutput~=3
    disp(['writing ',outfile])
    print(hf,flag,outfile,'-r150')
  end
end

end % function
