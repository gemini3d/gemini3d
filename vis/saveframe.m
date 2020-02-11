function saveframe(flagoutput, direc, filename, saveplot_fmt, h)
%% CREATES IMAGE FILES FROM PLOTS
narginchk(5,5)
validateattributes(flagoutput, {'numeric'}, {'scalar'}, mfilename)
validateattributes(direc, {'char'}, {'vector'}, mfilename)
validateattributes(filename, {'char'}, {'vector'}, mfilename)
validateattributes(h, {'struct'}, {'vector'}, mfilename, 'figure handles', 5)

if ischar(saveplot_fmt), saveplot_fmt = {saveplot_fmt}; end

assert(is_folder(direc), [direc, ' is not a directory.'])

% filename has the suffix, let's ditch the suffix.
[~, stem] = fileparts(filename);

plotdir = [direc, '/plots'];

makedir(plotdir)

disp(['writing plots to ', plotdir])

for i=1:length(saveplot_fmt)

  switch saveplot_fmt{i}
    case 'png', flag = '-dpng'; suffix = '.png';
    case 'eps', flag = '-depsc2'; suffix = '.eps';
    otherwise, continue
  end

  if flagoutput~=3
    print(h.f1,flag,[plotdir, '/v1-', stem, suffix],'-r300')
    print(h.f2,flag,[plotdir, '/Ti-', stem, suffix],'-r300')
    print(h.f3,flag,[plotdir, '/Te-', stem, suffix],'-r300')
    print(h.f4,flag,[plotdir, '/J1-', stem, suffix],'-r300')
    print(h.f5,flag,[plotdir, '/v2-', stem, suffix],'-r300')
    print(h.f6,flag,[plotdir, '/v3-', stem, suffix],'-r300')
    print(h.f7,flag,[plotdir, '/J2-', stem, suffix],'-r300')
    print(h.f8,flag,[plotdir, '/J3-', stem, suffix],'-r300')
    if ~isempty(h.f9)
      print(h.f9,flag,[plotdir, '/Phitop-', stem, suffix],'-r300')
    end
  end
  print(h.f10,flag,[plotdir, '/ne-', stem, suffix],'-r300')
end % for

end % function
