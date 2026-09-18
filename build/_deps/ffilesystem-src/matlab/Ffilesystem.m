%% uses "ffilesystem.{dll,so,dylib}" shared library created by Matlab-compatible compiler by:
%
%   cmake --workflow shared-nofortran-build

classdef Ffilesystem < handle
properties
  functions
  lib
  max_path
end

methods
  function obj = Ffilesystem(lib, hdr)
    arguments
      lib {mustBeTextScalar,mustBeFile}
      hdr {mustBeTextScalar,mustBeFile}
    end

    [~, obj.lib] = fileparts(lib);

    if ~libisloaded(obj.lib)
      loadlibrary(lib, hdr);
    end
    obj.max_path = calllib(obj.lib, "fs_get_max_path");
    obj.functions = libfunctions(obj.lib);
  end

    % https://www.mathworks.com/help/matlab/matlab_oop/handle-class-destructors.html
  function delete(obj)
    if ~isempty(obj.lib)
      unloadlibrary(obj.lib);
    end
  end

  function pid = getpid(obj)
    pid = calllib(obj.lib, 'fs_getpid');
  end

  function o = is_optimized(obj)
    o = calllib(obj.lib, 'fs_is_optimized');
  end

  function t = backend(obj)
    t = obj.string_call("fs_backend");
  end

  function t = compiler(obj)
    t = obj.string_call("fs_compiler");
  end

  function t = cpu_arch(obj)
    t = obj.string_call("fs_cpu_arch");
  end

  function t = exe_path(obj)
    t = obj.string_call("fs_exe_path");
  end

  function t = lib_path(obj)
    t = obj.string_call("fs_lib_path");
  end

  function t = shell(obj)
    t = obj.string_call("fs_get_shell");
  end

  function t = terminal(obj)
    t = obj.string_call("fs_get_terminal");
  end

  function s = file_size(obj, file)
    s = calllib(obj.lib, "fs_file_size", file);
  end

  function t = filesystem_type(obj, file)
    t = obj.string_call("fs_filesystem_type", file);
  end

  function y = equivalent(obj, file, file2)
    y = calllib(obj.lib, "fs_equivalent", file, file2);
  end
end

methods (Access = private)

function t = string_call(obj, fun, path1, path2)
arguments
  obj
  fun {mustBeTextScalar}
  path1 {mustBeTextScalar} = ''
  path2 {mustBeTextScalar} = ''
end

  s = string(pad("", obj.max_path));
  buf = libpointer('string', s);

  if strlength(path2)
    L = calllib(obj.lib, fun, char(path1), char(path2), obj.max_path);
  elseif strlength(path1)
    L = calllib(obj.lib, fun, char(path1), buf, obj.max_path);
  else
    L = calllib(obj.lib, fun, buf, obj.max_path);
  end

  if L > 0
    t = buf.value;
  else
    t = string.empty;
  end
end

end
end
