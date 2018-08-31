addpath ../script_utils/;
addpath ../vis/;


%READ IN THE SIMULATION INFORMATION
ID='~/simulations/input/KHI_periodic_highres/';
%ID='~/simulations/input/GDI_periodic/';
xg=readgrid(ID);


%LOAD THE FRAME OF THE SIMULATION THAT WE WANT TO PERTURB
direc=ID;
filebase='KHI_periodic_highres';
filename=[filebase,'_ICs.dat'];
loadframe3Dcurvnoelec;


%DEFINE A PERTURBATION AND CHANGE THE INITIAL CONDITIONS
%{
%%GDI nonperiodic
sigx2=30e3;
meanx3=0e3;
sigx3=30e3;
meanx2=-30e3;

for isp=1:lsp
  for ix3=1:xg.lx(3)
    for ix2=1:xg.lx(2)
      amplitude=rand(xg.lx(1),1);
      amplitude=0.1*amplitude;
      nsperturb(:,ix2,ix3,isp)=ns(:,ix2,ix3,isp)+ ...                                           %original data
                amplitude.*ns(:,ix2,ix3,isp)+ ...                                    %noise
                7.5d0*ns(:,ix2,ix3,isp).*exp(-1d0*(xg.x2(2+ix2)-meanx2).^18/2d0/sigx2.^18).* ...
                exp(-1d0*(xg.x3(2+ix3)-meanx3).^18/2d0/sigx3.^18);    %patch, note offset in the x2 index!!!!
    end
  end
end
%}


%%GDI EXAMPLE (PERIODIC)
%{
sigx2=20e3;
meanx3=0e3;
sigx3=20e3;
meanx2=-50e3;

for isp=1:lsp
  for ix2=1:xg.lx(2)
    amplitude=rand(xg.lx(1),1,xg.lx(3));
    amplitude=0.1*amplitude;
    nsperturb(:,ix2,:,isp)=ns(:,ix2,:,isp)+...                                           %original data
                4d0*ns(:,ix2,:,isp).*exp(-1d0*(xg.x2(2+ix2)-meanx2).^18/2d0/sigx2.^18)+...    %patch, note offset in the x2 index!!!!
                amplitude.*ns(:,ix2,:,isp);                                    %noise
  end
end
%nsperturb=max(nsperturb,1e4);
%}

%%KHI EXAMPLE

v0=500d0;
vn=500d0;
voffset=100d0;

sigx2=770e0;      %from Keskinen, 1988 growth rate formulas
meanx3=0e3;
sigx3=20e3;
meanx2=0e0;

for isp=1:lsp
  for ix2=1:xg.lx(2)
    amplitude=rand(xg.lx(1),1,xg.lx(3));
    amplitude=0.025*amplitude;
    nsperturb(:,ix2,:,isp)=ns(:,ix2,:,isp).*(v0+vn+voffset)./(-v0*tanh((xg.x2(2+ix2))/sigx2)+vn+voffset)+ ...
                          amplitude.*10e0.*ns(:,ix2,:,isp);
  end
end


%DIAGNOSTIC PLOT
imagesc(xg.x2(2:end-2),xg.x3(2:end-2),squeeze(nsperturb(50,:,:,7)));
axis xy;
xlabel('x_3');
ylabel('x_2');
title('Perturbed initial condition');
colorbar;


%WRITE OUT THE RESULTS TO A NEW FILE
dmy=[simdate(3),simdate(2),simdate(1)];
UTsec=simdate(4)*3600;
writedata(dmy,UTsec,nsperturb,vs1,Ts,[filebase,'_perturb']);


%RESET PATH
rmpath ../script_utils/;
rmpath ../vis/
