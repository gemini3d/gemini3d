function hdir = homepath()

persistent h;

if isempty(h)
    if ispc % windows
        h = [getenv('HOMEDRIVE'), getenv('HOMEPATH')];
    else %linux,mac
        h = getenv('HOME');
    end
end

hdir = h;  % for Matlab

end