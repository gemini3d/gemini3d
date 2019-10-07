function saveframe(flagoutput, direc, filename, saveplot_fmt, h)
%% CREATES IMAGE FILES FROM PLOTS
narginchk(5,5)
validateattributes(flagoutput, {'numeric'}, {'scalar'}, mfilename)
validateattributes(direc, {'char'}, {'vector'}, mfilename)
validateattributes(filename, {'char'}, {'vector'}, mfilename)
validateattributes(saveplot_fmt, {'cell'}, {'vector'}, mfilename)
validateattributes(h, {'struct'}, {'vector'}, mfilename)

dirs = {'v1plots', 'Tiplots', 'Teplots', 'J1plots', 'v2plots', 'v3plots', 'J2plots', 'J3plots', 'Phiplots', 'nplots'};
for i = 1:length(dirs)
  dirs{i} = [direc, filesep, dirs{i}];
  if ~is_folder(dirs{i})
    mkdir(dirs{i});
  end
end

disp(['writing plots to ',direc])

for i=1:length(saveplot_fmt)

  switch saveplot_fmt{i}
    case 'png', flag = '-dpng'; suffix = '.png';
    case 'eps', flag = '-depsc2'; suffix = '.eps';
    otherwise, continue
  end

  if flagoutput~=3
    print(h.f1,flag,[dirs{1}, filesep, filename,suffix],'-r300')
    print(h.f2,flag,[dirs{2}, filesep, filename,suffix],'-r300')
    print(h.f3,flag,[dirs{3}, filesep, filename,suffix],'-r300')
    print(h.f4,flag,[dirs{4}, filesep, filename,suffix],'-r300')
    print(h.f5,flag,[dirs{5}, filesep, filename,suffix],'-r300')
    print(h.f6,flag,[dirs{6}, filesep, filename,suffix],'-r300')
    print(h.f7,flag,[dirs{7}, filesep, filename,suffix],'-r300')
    print(h.f8,flag,[dirs{8}, filesep, filename,suffix],'-r300')
    if ~isempty(h.f9)
      print(h.f9,flag,[dirs{9}, filesep, filename,suffix],'-r300')
    end
  end
  print(h.f10,flag,[dirs{10}, filesep, filename,suffix],'-r300')
end % for

end % function
