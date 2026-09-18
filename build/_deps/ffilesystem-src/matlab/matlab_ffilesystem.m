%% Example of using Ffilesystem Matlab class
% we haven't added all the Ffilesystem function as methods yet.
% However as seen in Ffilesystem.m, it is straightforward to add function methods.
%
% First, build Ffilesystem shared library with CMake
%   cmake --workflow shared-nofortran-build

function matlab_ffilesystem()

libdir = ('../build-matlab');

% create an instance of Ffilesystem class
lib = get_lib(libdir);
inc = '../include/ffilesystem.h';

fs = Ffilesystem(lib, inc);

cpu_arch = fs.cpu_arch();
disp("CPU arch: " + cpu_arch)

backend = fs.backend();
disp("Ffilesystem backend: " + backend)

compiler = fs.compiler();
disp("Compiler: " + compiler)

s = fs.file_size('Readme.md');
disp("Readme.md file size (bytes): " + string(s))

fs_type = fs.filesystem_type(pwd());
disp("Filesystem type: " + fs_type)

e = fs.equivalent('.', pwd());
disp("Equivalent . and pwd: " + e)

% close the Ffilesystem instance
% this unloads the C library
delete(fs);

end


function lib = get_lib(libdir)
arguments (Input)
  libdir (1,1) string {mustBeFolder}
end
arguments (Output)
  lib (1,1) string {mustBeFile}
end

if ispc()
  names = fullfile(libdir, ["", "Release"], "ffilesystem.dll");
  i = find(isfile(names), 1);
  lib = names(i);
elseif ismac()
  lib = fullfile(libdir, "libffilesystem.dylib");
else
  lib = fullfile(libdir, "libffilesystem.so");
end

end
