function xg = plotframe(direc,ymd,UTsec,saveplots,plotfun,xg)

narginchk(3,6)
if nargin<4
  saveplots=false;
end


%SET THE CAXIS LIMITS FOR THE PLOTS
nelim =  [0 1.9e12];
v1lim = [-400 400];
Tilim = [100 3000];
Telim = [100 3000];
J1lim = [-0.1 0.1];
v2lim = [-10 10];
v3lim = [-10 10];
J2lim = [-0.1 0.1];
J3lim=[-0.25 0.25];


%MAKE DIRECTORIES FOR OUTPUT FILES
if saveplots
  mkdir([direc,'/nplots']);    %store output plots with the simulation data
  mkdir([direc,'/v1plots']);    
  mkdir([direc,'/v2plots']);    
  mkdir([direc,'/v3plots']);   
  mkdir([direc,'/J1plots']);   
  mkdir([direc,'/Tiplots']);   
  mkdir([direc,'/Teplots']);   
  mkdir([direc,'/J2plots']);  
  mkdir([direc,'/J3plots']);  
  mkdir([direc,'/Phiplots']); 
  mkdir([direc,'/nplots_eps']);    %store output plots with the simulation data
  mkdir([direc,'/v1plots_eps']);
  mkdir([direc,'/v2plots_eps']);
  mkdir([direc,'/v3plots_eps']);
  mkdir([direc,'/J1plots_eps']);
  mkdir([direc,'/Tiplots_eps']);
  mkdir([direc,'/Teplots_eps']);
  mkdir([direc,'/J2plots_eps']);
  mkdir([direc,'/J3plots_eps']);
  mkdir([direc,'/Phiplots_eps']);
end

cwd = fileparts(mfilename('fullpath'));

addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

%READ IN THE SIMULATION INFORMATION (this is low cost so reread no matter what)
fprintf('Reloading input file...\n');
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs/config.ini']);


%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if (~exist('xg','var'))
  disp('Reloading grid...')
  xg = readgrid([direc,filesep,'inputs',filesep]);
end


% DEFINE THE PLOTTING FUNCTION BASED ON THE TYPE OF GRID USED
if (~exist('plotfun','var') | isempty(plotfun))
  minh1=min(xg.h1(:));
  maxh1=max(xg.h1(:));
  if (abs(minh1-1)>1e-4 | abs(maxh1-1)>1e-4)    %curvilinear grid
    if (xg.lx(2)>1 & xg.lx(3)>1)
      plotfun=@plot3D_curv_frames_long;
    else
      plotfun=@plot2D_curv;
    end
  else     %cartesian grid
    if (xg.lx(2)>1 & xg.lx(3)>1)
      plotfun=@plot3D_cart_frames_long_ENU;
    else
      plotfun=@plot2D_cart;
    end 
  end
end


% COMPUTE SOURUCE LOCATION IN MCOORDS
if ~isempty(mloc)
 mlat=mloc(1);
 mlon=mloc(2);
else
 mlat=[];
 mlon=[];
end


% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
lt=numel(times);

Nsp = ceil(sqrt(length(times)));    %number of subplots (if used)


%CREATE HANDLES FOR FIGURE SET
h1=figure('name','V1');
h2=figure('name','Ti');
h3=figure('name','Te');
h4=figure('name','J1');
h5=figure('name','V2');
h6=figure('name','V3');
h7=figure('name','J2');
h8=figure('name','J3');
if (xg.lx(2)~=1 & xg.lx(3)~=1)
  h9=figure('name','phiTop');
end
h10=figure('name','Ne');


%DETERMINE OUTPUT TIME NEAREST TO THAT REQUESTED
tnow=UTsec0+tdur;    %not date corrected, but doesn't really matter
[ymdend,UTsecend]=dateinc(tdur,ymd0,UTsec0);


%LOCATE TIME NEAREST TO THE REQUESTED DATE
ymdprev=ymd0;
UTsecprev=UTsec0;
ymdnext=ymd0;
UTsecnext=UTsec0;
it=1;
while(datenum([ymdnext,UTsecnext/3600,0,0])<datenum([ymd,UTsec/3600,0,0]) & datenum([ymdnext,UTsecnext/3600,0,0])<datenum([ymdend,UTsecend/3600,0,0]))
  ymdprev=ymdnext;
  UTsecprev=UTsecnext;
  [ymdnext,UTsecnext]=dateinc(dtout,ymdprev,UTsecprev);
  it=it+1;
end
if(UTsecnext-UTsec<=UTsec-UTsecprev)   %we are closer to next frame
  ymdnow=ymdnext;
  UTsecnow=UTsecnext;
else    %closer to previous frame
  ymdnow=ymdprev;
  UTsecnow=UTsecprev;
end


%LOAD THE FRAME NEAREST TO THE REQUESTED TIME
[ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,ymd0,UTsec0,tdur,dtout,flagoutput,mloc,xg);


%MAKE THE PLOTS (WHERE APPROPRIATE)
if (xg.lx(2)~=1 & xg.lx(3)~=1 | lt>16)    %3D simulation or a very long 2D simulation - do separate plots for each time frame
    fprintf('Detected a long 2D or 3D simulation...\n');
    h10=figure(h10);
%        plotfun(ymd,UTsec,xg,log10(ne(:,:,:)),'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h10);
    plotfun(ymd,UTsec,xg,ne(:,:,:),'n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h10);
    
    if (flagoutput~=3)
        h1=figure(h1);
        plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',v1lim,[mlatsrc,mlonsrc],h1);
        
        h2=figure(h2);
        plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',Tilim,[mlatsrc,mlonsrc],h2);
        
        h3=figure(h3);
        plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',Telim,[mlatsrc,mlonsrc],h3);
        
        h4=figure(h4);
        plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',J1lim,[mlatsrc,mlonsrc],h4);
        
        h5=figure(h5);
        plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],h5);
        
        h6=figure(h6);
        plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],h6);
        
        h7=figure(h7);
        plotfun(ymd,UTsec,xg,J2(:,:,:)*1e6,'J_2 (uA/m^2)',J2lim,[mlatsrc,mlonsrc],h7);
        
        h8=figure(h8);
        plotfun(ymd,UTsec,xg,J3(:,:,:)*1e6,'J_3 (uA/m^2)',J3lim,[mlatsrc,mlonsrc],h8);
        
        if (xg.lx(2)~=1 & xg.lx(3)~=1)
            h9=figure(h9);
            h9a=gca;
            imagesc(h9a,Phitop)
            colorbar;
        end
    end
    if (saveplots)    %for 3D or long 2D plots print and output file every time step
        if (flagoutput~=3)
            print(h1,'-dpng',[direc,'/v1plots/',filename,'.png'],'-r300')
            print(h2,'-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300')
            print(h3,'-dpng',[direc,'/Teplots/',filename,'.png'],'-r300')
            print(h4,'-dpng',[direc,'/J1plots/',filename,'.png'],'-r300')
            print(h5,'-dpng',[direc,'/v2plots/',filename,'.png'],'-r300')
            print(h6,'-dpng',[direc,'/v3plots/',filename,'.png'],'-r300')
            print(h7,'-dpng',[direc,'/J2plots/',filename,'.png'],'-r300')
            print(h8,'-dpng',[direc,'/J3plots/',filename,'.png'],'-r300')
            if (xg.lx(2)~=1 & xg.lx(3)~=1)
                print(h9,'-dpng',[direc,'/Phiplots/',filename,'.png'],'-r300')
            end
        end
        print(h10,'-dpng',[direc,'/nplots/',filename,'.png'],'-r300')

        if (flagoutput~=3)     %now make .eps prints of the plots
            print(h1,'-depsc2',[direc,'/v1plots_eps/',filename,'.eps'])
            print(h2,'-depsc2',[direc,'/Tiplots_eps/',filename,'.eps'])
            print(h3,'-depsc2',[direc,'/Teplots_eps/',filename,'.eps'])
            print(h4,'-depsc2',[direc,'/J1plots_eps/',filename,'.eps'])
            print(h5,'-depsc2',[direc,'/v2plots_eps/',filename,'.eps'])
            print(h6,'-depsc2',[direc,'/v3plots_eps/',filename,'.eps'])
            print(h7,'-depsc2',[direc,'/J2plots_eps/',filename,'.eps'])
            print(h8,'-depsc2',[direc,'/J3plots_eps/',filename,'.eps'])
            if (xg.lx(2)~=1 & xg.lx(3)~=1)
                print(h9,'-depsc2',[direc,'/Phiplots_eps/',filename,'.eps'])
            end
        end
        print(h10,'-depsc2',[direc,'/nplots_eps/',filename,'.eps'])
    end
else    %short 2D simulation - put the entire time series in a single plot
    fprintf('Detected a short 2D simulations...\n');
    h10=figure(h10);
    ha = subplot(Nsp,Nsp,it,'parent',h10);
    nelim =  [9 11.3];
    plotfun(ymd,UTsec,xg,log10(ne(:,:,:)),'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],ha);
    
    if (flagoutput~=3)
        h1=figure(h1);
        ha = subplot(Nsp,Nsp,it,'parent',h1);
        plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',v1lim,[mlatsrc,mlonsrc],ha);
        
        h2=figure(h2);
        ha = subplot(Nsp,Nsp,it,'parent',h2);
        plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',Tilim,[mlatsrc,mlonsrc],ha);
        
        h3=figure(h3);
        ha = subplot(Nsp,Nsp,it,'parent',h3);
        plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',Telim,[mlatsrc,mlonsrc],ha);
        
        h4=figure(h4);
        ha = subplot(Nsp,Nsp,it,'parent',h4);
        plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',J1lim,[mlatsrc,mlonsrc],ha);
        
        h5=figure(h5);
        ha = subplot(Nsp,Nsp,it,'parent',h5);
        plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],ha);
        
        h6=figure(h6);
        ha = subplot(Nsp,Nsp,it,'parent',h6);
        plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],ha);
        
        h7=figure(h7);
        ha = subplot(Nsp,Nsp,it,'parent',h7);
        plotfun(ymd,UTsec,xg,J2(:,:,:)*1e6,'J_2 (uA/m^2)',J2lim,[mlatsrc,mlonsrc],ha);
        
        h8=figure(h8);
        ha = subplot(Nsp,Nsp,it,'parent',h8);
        plotfun(ymd,UTsec,xg,J3(:,:,:)*1e6,'J_3 (uA/m^2)',J3lim,[mlatsrc,mlonsrc],ha);
        
        if (xg.lx(2)~=1 & xg.lx(3)~=1)
            h9=figure(h9);
            ha = subplot(Nsp,Nsp,it,'parent',h9);
            imagesc(ha,Phitop)
            colorbar;
        end
    end
end


if (saveplots & (xg.lx(2)==1 & xg.lx(3)==1 | lt<=16) )    %save the short 2D sim plots
    if (flagoutput~=3)
        print(h1,'-dpng',[direc,'/v1plots/',filename,'.png'],'-r300')
        print(h2,'-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300')
        print(h3,'-dpng',[direc,'/Teplots/',filename,'.png'],'-r300')
        print(h4,'-dpng',[direc,'/J1plots/',filename,'.png'],'-r300')
        print(h5,'-dpng',[direc,'/v2plots/',filename,'.png'],'-r300')
        print(h6,'-dpng',[direc,'/v3plots/',filename,'.png'],'-r300')
        print(h7,'-dpng',[direc,'/J2plots/',filename,'.png'],'-r300')
        print(h8,'-dpng',[direc,'/J3plots/',filename,'.png'],'-r300')
        if (xg.lx(2)~=1 & xg.lx(3)~=1)
            print(h9,'-dpng',[direc,'/Phiplots/',filename,'.png'],'-r300')
        end
    end
    print(h10,'-dpng',[direc,'/nplots/',filename,'.png'],'-r300')

    if (flagoutput~=3)     %now make .eps prints of the plots
        print(h1,'-depsc2',[direc,'/v1plots_eps/',filename,'.eps'])
        print(h2,'-depsc2',[direc,'/Tiplots_eps/',filename,'.eps'])
        print(h3,'-depsc2',[direc,'/Teplots_eps/',filename,'.eps'])
        print(h4,'-depsc2',[direc,'/J1plots_eps/',filename,'.eps'])
        print(h5,'-depsc2',[direc,'/v2plots_eps/',filename,'.eps'])
        print(h6,'-depsc2',[direc,'/v3plots_eps/',filename,'.eps'])
        print(h7,'-depsc2',[direc,'/J2plots_eps/',filename,'.eps'])
        print(h8,'-depsc2',[direc,'/J3plots_eps/',filename,'.eps'])
        if (xg.lx(2)~=1 & xg.lx(3)~=1)
            print(h9,'-depsc2',[direc,'/Phiplots_eps/',filename,'.eps'])
        end
    end
    print(h10,'-depsc2',[direc,'/nplots_eps/',filename,'.eps'])
end
    
end % function
