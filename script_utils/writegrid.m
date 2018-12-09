function writegrid(xg,outdir)

%--------------------------------------------------------
%-----NOTE THAT WE DO STORE *ALL*
%-----GRID DATA (INCLUDING STUFF NOT NEEDED BY
%-----FORTRAN CODE BUT POSSIBLY USEFUL FOR PLOTTING)
%--------------------------------------------------------
 
  %% MAKE THE OUTPUT DIRECTORY IF IT DOESN'T EXIST AND NOTIFY USER
  if (~(exist(outdir,'dir')==7))
    mkdir(outdir);
    disp(['Note that directory:  ',outdir,' needed to be created for grid file...']);
  end


  filename=[outdir,'/simsize.dat'];
  fid=fopen(filename,'w');
  fwrite(fid,xg.lx,'integer*4');
  fclose(fid);
  
  fid=fopen([outdir,'/simgrid.dat'],'w');
  
  fwrite(fid,xg.x1,'real*8');    %coordinate values
  fwrite(fid,xg.x1i,'real*8');
  fwrite(fid,xg.dx1b,'real*8');
  fwrite(fid,xg.dx1h,'real*8');
  
  fwrite(fid,xg.x2,'real*8');
  fwrite(fid,xg.x2i,'real*8');
  fwrite(fid,xg.dx2b,'real*8');
  fwrite(fid,xg.dx2h,'real*8');
  
  fwrite(fid,xg.x3,'real*8');
  fwrite(fid,xg.x3i,'real*8');
  fwrite(fid,xg.dx3b,'real*8');
  fwrite(fid,xg.dx3h,'real*8');
  
  fwrite(fid,xg.h1,'real*8');   %cell-centered metric coefficients
  fwrite(fid,xg.h2,'real*8');
  fwrite(fid,xg.h3,'real*8');

  fwrite(fid,xg.h1x1i,'real*8');    %interface metric coefficients
  fwrite(fid,xg.h2x1i,'real*8');
  fwrite(fid,xg.h3x1i,'real*8');
  
  fwrite(fid,xg.h1x2i,'real*8');
  fwrite(fid,xg.h2x2i,'real*8');
  fwrite(fid,xg.h3x2i,'real*8');
  
  fwrite(fid,xg.h1x3i,'real*8');
  fwrite(fid,xg.h2x3i,'real*8');
  fwrite(fid,xg.h3x3i,'real*8');
  
  %gravity, geographic coordinates, magnetic field strength? unit vectors?
  fwrite(fid,xg.gx1,'real*8');    %gravitational field components
  fwrite(fid,xg.gx2,'real*8');
  fwrite(fid,xg.gx3,'real*8');
  
  fwrite(fid,xg.alt,'real*8');    %geographic coordinates
  fwrite(fid,xg.glat,'real*8');
  fwrite(fid,xg.glon,'real*8');
  
  fwrite(fid,xg.Bmag,'real*8');    %magnetic field strength
  
  fwrite(fid,xg.I,'real*8');    %magnetic field inclination
  
  fwrite(fid,xg.nullpts,'real*8');    %points not to be solved
  
  
  %NOT ALL OF THE REMAIN INFO IS USED IN THE FORTRAN CODE, BUT IT INCLUDED FOR COMPLETENESS
  fwrite(fid,xg.e1,'real*8');   %4D unit vectors (in cartesian components)
  fwrite(fid,xg.e2,'real*8');
  fwrite(fid,xg.e3,'real*8');
  
  fwrite(fid,xg.er,'real*8');    %spherical unit vectors
  fwrite(fid,xg.etheta,'real*8');
  fwrite(fid,xg.ephi,'real*8');
  
  fwrite(fid,xg.r,'real*8');    %spherical coordinates
  fwrite(fid,xg.theta,'real*8');
  fwrite(fid,xg.phi,'real*8');
  
  fwrite(fid,xg.x,'real*8');     %cartesian coordinates
  fwrite(fid,xg.y,'real*8');
  fwrite(fid,xg.z,'real*8');  
 
  fclose(fid);  
end
