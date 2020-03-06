clear,clc;
fid=fopen('MHD1D.dat');
lz=fscanf(fid,'%f',1);
z=fscanf(fid,'%f \n',lz+4);

while ~feof(fid)
  t=fscanf(fid,'%f',1);
  if feof(fid)
    break;
  end

  n=fscanf(fid,'%f',lz+4);
  u1=fscanf(fid,'%f',lz+4);
  u2=fscanf(fid,'%f',lz+4);
  u3=fscanf(fid,'%f',lz+4);
  B1=fscanf(fid,'%f',lz+4);
  B2=fscanf(fid,'%f',lz+4);
  B3=fscanf(fid,'%f',lz+4);
  p=fscanf(fid,'%f',lz+4);

  plotMHD(n,u2,u3,u1,B2,B3,B1,p,z,t)
  pause%(0.05)
end

fclose(fid);
