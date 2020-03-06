function plotMHD(n,ux,uy,uz,Bx,By,Bz,p,z,t)

subplot(411);
plot(z/1e3,n);
maxn=max(n);
minn=min(n);
dn=max(maxn-minn,1);
meann=mean(n);
axis([min(z)/1e3 max(z)/1e3 meann-0.6*dn meann+0.6*dn]);
xlabel('z [km]')
ylabel('n [m^{-3}]');
title(sprintf('t=%.1f',t));

subplot(412);
plot(z/1e3,[ux/1e3,uy/1e3,uz/1e3]);
maxu=max(abs([ux(:);uy(:);uz(:)]))/1e3;
maxu=max(maxu,1e-3);
axis([min(z)/1e3 max(z)/1e3 -maxu-0.2*maxu maxu+0.2*maxu]);
xlabel('z [km]')
ylabel('u [km/s]');

subplot(413);
maxB=max([Bx(:);By(:)])/1e-9;
maxB=max(maxB,1e-3);
plot(z/1e3,[Bx/1e-9,By/1e-9,Bz/1e-9]);
axis([min(z)/1e3 max(z)/1e3 -maxB-0.2*maxB maxB+0.2*maxB]);
xlabel('z [km]')
ylabel('B [nT]');
%legend('x','y','z')

subplot(414);
plot(z/1e3,p);
maxn=max(p);
minn=min(p);
dn=max(maxn-minn,1e-20);
meann=mean(p);
axis([min(z)/1e3 max(z)/1e3 meann-0.6*dn meann+0.6*dn]);
xlabel('z [km]')
ylabel('p [J*m^{-3}]');

end
