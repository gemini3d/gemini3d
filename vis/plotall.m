function xg = plotall(direc, saveplots, plotfun, xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(1,3)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)
if nargin<2, saveplots=false; end
validateattr(saveplots, {'logical'}, {'scalar'}, mfilename)
if nargin<3
  plotfun=[]; 
else
  validateattr(plotfun, {'function_handle'}, {'scalar'}, mfilename, 'plotting function handle', 3)
end
if nargin<4
  xg=[]; 
else
  validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 4)
end

%SET THE CAXIS LIMITS FOR THE PLOTS
%nelim =  [9 11.3];
nelim =  [0 8e11];
v1lim = [-300 300];
Tilim = [100 3000];
Telim = [100 3000];
%J1lim = [-0.15 0.15];
J1lim = [-40 40];
v2lim = [-1500 1500];
v3lim = [-1500 1500];
%J2lim = [-0.25 0.25];
J2lim = [-10 10];
%J3lim = [-0.5 0.5];
J3lim=[-10 10];


%% MAKE DIRECTORIES FOR OUTPUT FILES
% store output plots with the simulation data
if saveplots
  dlist = {'nplots', 'v1plots', 'v2plots', 'v3plots', 'J1plots',...
           'Tiplots', 'Teplots', 'J2plots', 'J3plots', 'Phiplots', ...
           'nplots_eps', 'v1plots_eps', 'v2plots_eps', 'v3plots_eps', ...
           'J1plots_eps', 'Tiplots_eps', 'Teplots_eps', 'J2plots_eps', 'J3plots_eps', ...
           'Phiplots_eps'};
  for i=1:length(dlist)
    mkdir([direc, filesep, dlist{i}]);
  end
end

%%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs', filesep, 'config.ini']);


%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if isempty(xg)
  xg = readgrid([direc,filesep,'inputs',filesep]);
end


%CUSTOM FUNCTINO TO PLOT WITH
%%plotfun=@plot2D_curv;
%%plotfun=@plot2D_cart;
%%plotfun=@plot2D_curv_north;
%%plotfun=@plot2D_curv_south;
%%plotfun=@plot3D_curv_frames;
%%plotfun=@plot3D_cart_frames;
%plotfun=@plot3D_curv_frames_long;
%%plotfun=@plot3D_cart_frames_long_ENU;


%% DEFINE THE PLOTTING FUNCTION BASED ON THE TYPE OF GRID USED

minh1=min(xg.h1(:));
maxh1=max(xg.h1(:));
if (abs(minh1-1)>1e-4 || abs(maxh1-1)>1e-4)    %curvilinear grid
  if xg.lx(2)>1 && xg.lx(3)>1
    plotfun=@plot3D_curv_frames_long;
  else
    plotfun=@plot2D_curv;
  end
else     %cartesian grid
  if xg.lx(2)>1 && xg.lx(3)>1
    plotfun=@plot3D_cart_frames_long_ENU;
  else
    plotfun=@plot2D_cart;
  end 
end


%% COMPUTE SOURUCE LOCATION IN MCOORDS
if ~isempty(mloc)
 mlat=mloc(1);
 mlon=mloc(2);
else
 mlat=[];
 mlon=[];
end


%% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
lt=numel(times);

Nsp = ceil(sqrt(length(times)));    %number of subplots (if used)
ymd=ymd0;
UTsec=UTsec0;
%% initialize figure sets
h1=figure('name','V1');
h2=figure('name','Ti');
h3=figure('name','Te');
h4=figure('name','J1');
h5=figure('name','V2');
h6=figure('name','V3');
h7=figure('name','J2');
h8=figure('name','J3');
if (xg.lx(2)~=1 && xg.lx(3)~=1)
  h9=figure('name','phiTop');
end
h10=figure('name','Ne');
%% main figure loop
for it=1:lt
    [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop] = loadframe(direc, UTsec, ymd, UTsec0, ymd0, mloc, xg);
    disp([filename, ' => ', func2str(plotfun)])
    

    
    %% Electron number density
    if xg.lx(2)~=1 && xg.lx(3)~=1 || lt>16    %3D simulation or a very long 2D simulation - do separate plots for each time frame
        disp('long 2D or 3D simulation...');
        h10=figure(h10);
%        plotfun(ymd,UTsec,xg,log10(ne(:,:,:)),'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h10);
        plotfun(ymd,UTsec,xg, ne(:,:,:),'n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h10);
        
        if flagoutput~=3
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
            
            if xg.lx(2)~=1 && xg.lx(3)~=1
                h9=figure(h9);
                h9a=gca;
                imagesc(h9a,Phitop)
                colorbar;
            end
        end
        if saveplots   % for 3D or long 2D plots print and output file every time step
            if flagoutput~=3
                print(h1,'-dpng',[direc,'/v1plots/',filename,'.png'],'-r300')
                print(h2,'-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300')
                print(h3,'-dpng',[direc,'/Teplots/',filename,'.png'],'-r300')
                print(h4,'-dpng',[direc,'/J1plots/',filename,'.png'],'-r300')
                print(h5,'-dpng',[direc,'/v2plots/',filename,'.png'],'-r300')
                print(h6,'-dpng',[direc,'/v3plots/',filename,'.png'],'-r300')
                print(h7,'-dpng',[direc,'/J2plots/',filename,'.png'],'-r300')
                print(h8,'-dpng',[direc,'/J3plots/',filename,'.png'],'-r300')
                if xg.lx(2)~=1 && xg.lx(3)~=1
                    print(h9,'-dpng',[direc,'/Phiplots/',filename,'.png'],'-r300')
                end
            end
            print(h10,'-dpng',[direc,'/nplots/',filename,'.png'],'-r300')

            if flagoutput~=3     %now make .eps prints of the plots
                print(h1,'-depsc2',[direc,'/v1plots_eps/',filename,'.eps'])
                print(h2,'-depsc2',[direc,'/Tiplots_eps/',filename,'.eps'])
                print(h3,'-depsc2',[direc,'/Teplots_eps/',filename,'.eps'])
                print(h4,'-depsc2',[direc,'/J1plots_eps/',filename,'.eps'])
                print(h5,'-depsc2',[direc,'/v2plots_eps/',filename,'.eps'])
                print(h6,'-depsc2',[direc,'/v3plots_eps/',filename,'.eps'])
                print(h7,'-depsc2',[direc,'/J2plots_eps/',filename,'.eps'])
                print(h8,'-depsc2',[direc,'/J3plots_eps/',filename,'.eps'])
                if xg.lx(2)~=1 && xg.lx(3)~=1
                    print(h9,'-depsc2',[direc,'/Phiplots_eps/',filename,'.eps'])
                end
            end
            print(h10,'-depsc2',[direc,'/nplots_eps/',filename,'.eps'])
        end
    else    %short 2D simulation - put the entire time series in a single plot
        disp('short 2D simulations...')
        h10=figure(h10);
        ha = subplot(Nsp,Nsp,it,'parent',h10);
        nelim =  [9 11.3];
        plotfun(ymd,UTsec,xg,log10(ne(:,:,:)),'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],ha);
        
        if flagoutput~=3
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
            
            if xg.lx(2)~=1 && xg.lx(3)~=1
                h9=figure(h9);
                ha = subplot(Nsp,Nsp,it,'parent',h9);
                imagesc(ha,Phitop)
                colorbar;
            end
        end
    end
    
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end % for
    
    if saveplots && (xg.lx(2)==1 && xg.lx(3)==1 || lt<=16)    %save the short 2D sim plots
        if flagoutput~=3
            print(h1,'-dpng',[direc,'/v1plots/',filename,'.png'],'-r300')
            print(h2,'-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300')
            print(h3,'-dpng',[direc,'/Teplots/',filename,'.png'],'-r300')
            print(h4,'-dpng',[direc,'/J1plots/',filename,'.png'],'-r300')
            print(h5,'-dpng',[direc,'/v2plots/',filename,'.png'],'-r300')
            print(h6,'-dpng',[direc,'/v3plots/',filename,'.png'],'-r300')
            print(h7,'-dpng',[direc,'/J2plots/',filename,'.png'],'-r300')
            print(h8,'-dpng',[direc,'/J3plots/',filename,'.png'],'-r300')
            if xg.lx(2)~=1 && xg.lx(3)~=1
                print(h9,'-dpng',[direc,'/Phiplots/',filename,'.png'],'-r300')
            end
        end
        print(h10,'-dpng',[direc,'/nplots/',filename,'.png'],'-r300')

        if flagoutput~=3     %now make .eps prints of the plots
            print(h1,'-depsc2',[direc,'/v1plots_eps/',filename,'.eps'])
            print(h2,'-depsc2',[direc,'/Tiplots_eps/',filename,'.eps'])
            print(h3,'-depsc2',[direc,'/Teplots_eps/',filename,'.eps'])
            print(h4,'-depsc2',[direc,'/J1plots_eps/',filename,'.eps'])
            print(h5,'-depsc2',[direc,'/v2plots_eps/',filename,'.eps'])
            print(h6,'-depsc2',[direc,'/v3plots_eps/',filename,'.eps'])
            print(h7,'-depsc2',[direc,'/J2plots_eps/',filename,'.eps'])
            print(h8,'-depsc2',[direc,'/J3plots_eps/',filename,'.eps'])
            if xg.lx(2)~=1 && xg.lx(3)~=1
                print(h9,'-depsc2',[direc,'/Phiplots_eps/',filename,'.eps'])
            end
        end
        print(h10,'-depsc2',[direc,'/nplots_eps/',filename,'.eps'])
    end
    
if nargout==0, clear('xg'), end
    
end % function
