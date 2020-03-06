function hdir = homepath()

persistent h;

if isempty(h)
    if ispc % windows
        h = getenv('USERPROFILE');
    else %linux,mac
        h = getenv('HOME');
    end
end

hdir = h;  % for Matlab

end
