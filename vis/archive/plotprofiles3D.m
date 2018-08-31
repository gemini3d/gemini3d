function h=plotprofiles3D(t,xg,ns,vsx1,Ts)
clf;
h=gcf;   %return the handle to whichever figure is to be used for plots

dmy=[0,0,0];


%PARAMETERS TO PLOT
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
lsp=size(ns,4);


%FIELD ALIGNED 'OUTFLOW'
%gx1s=repmat(xg.gx1,[1, 1, lsp]);
%vsfa=-1*vsx1.*sign(gx1s);
if (xg.alt(2,1,1)>xg.alt(1,1,1))    %non-inverted grid type
  vsfa=vsx1;
else
  vsfa=-1*vsx1;
  fprintf('\n  PLOTPROFILES3D --> Detected an inverted-type grid...\n')
end


%SUBSET OF DATA
ix2=floor(lx2/4);
ix3=floor(lx3/4);
%if abs(xg.r(1,1)-xg.r(lx1,1))<1         %closed dipole grid
if (xg.alt(1,1,1)-xg.alt(lx1,1,1) < 1)    %closed dipole grid
    inds=1:floor(lx1/2);    %southern hemisphere on a closed grid
    z=xg.alt(inds,ix2,ix3)/1e3;
   fprintf('\n  PLOTPROFILES3D --> Detected an closed, dipole grid...\n')
else
    z=xg.alt(:,ix2,ix3)/1e3;
end
%maxz=max(z);
maxz=600;
izs=find(z>=80);


%PLOT
lege={'O^+','NO^+','N_2^+','O_2^+','N^+','H^+','e^-'};

subplot(131);
semilogx(squeeze(ns(izs,ix2,ix3,:)),z(izs));
axis([1e0 1e12 80 maxz]);
ax=axis;
xlabel('density [m^{-3}]');
ylabel('altitude [km]');
legend(lege,'Location','SouthWest');
text(ax(1),ax(4)-(ax(4)-ax(3))/20,[num2str(dmy(1)),'/',num2str(dmy(2)),'/',num2str(dmy(3)), ...
    ' ',num2str(t),' UT'],'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');

subplot(132);
plot(squeeze(Ts(izs,ix2,ix3,:)),z(izs));
axis([100 6000 80 maxz]);
xlabel('temperature [K]');
ylabel('altitude [km]');
%legend(lege);

subplot(133);
plot(squeeze(vsfa(izs,ix2,ix3,:)),z(izs));
axis([-500 500 80 maxz]);
xlabel('drift velocity [m/s]');
ylabel('altitude [km]');
%legend(lege,'Location','SouthEast');

end
