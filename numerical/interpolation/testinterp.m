function testinterp(filename)
  fid=fopen(filename,'r');

  lx1=fread(fid,1,'integer*4');
  lx2=fread(fid,1,'integer*4');
  x1=fread(fid,lx1,'real*8');
  x2=fread(fid,lx2,'real*8');
  f=fread(fid,lx1*lx2,'real*8');
  f=reshape(f,[lx1, lx2]);

  fclose(fid);


  %% PLOT
  figure

  if (lx2==1)
    plot(x1,f);
    xlabel('x_1')
    ylabel('f')
    title('1-D interp')
  else
    imagesc(x2,x1,f);
    axis xy;
    xlabel('x_2')
    ylabel('x_1')
    c=colorbar;
    ylabel(c,'f')
    title('2-D interp')
  end

  %print -dpng -r300 ~/testinterp.png;

end % function
