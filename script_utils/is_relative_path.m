function isrel = is_relative_path(path)
%% detect if a path is relative

isrel = strcmp(path(1), '.');
if isrel, return, end

isrel = ~contains(path, filesep);
if isrel, return, end

isrel =  isletter(path(1)) && ~strcmp(path(2), ':');

end
