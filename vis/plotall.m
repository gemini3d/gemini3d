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
           'Tiplots', 'Teplots', 'J2plots', 'J3plots', 'Phiplots'};
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
h.f1=figure('name','V1');
h.f2=figure('name','Ti');
h.f3=figure('name','Te');
h.f4=figure('name','J1');
h.f5=figure('name','V2');
h.f6=figure('name','V3');
h.f7=figure('name','J2');
h.f8=figure('name','J3');
if xg.lx(2)>1 && xg.lx(3)>1 % a 3-D simulation
  h.f9=figure('name','phiTop');
else
  h.f9 = [];
end
h.f10=figure('name','Ne');%, 'position', [.1, .1, .5, .5], 'units', 'normalized');
%% main figure loop
for it=1:lt
    [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop] = loadframe(direc, UTsec, ymd, UTsec0, ymd0, mloc, xg);
    disp([filename, ' => ', func2str(plotfun)])

    %% Electron number density, 'position', [.1, .1, .5, .5], 'units', 'normalized'
    if ~isempty(h.f9) || lt>16    % 3D simulation or a very long 2D simulation - do separate plots for each time frame
        if it==1, disp('long 2D or 3D simulation...'), end
%        plotfun(ymd,UTsec,xg,log10(ne(:,:,:)),'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h10);
        plotfun(ymd,UTsec,xg, ne, 'n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],h.f10);
        
        if flagoutput~=3
            plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',v1lim,[mlatsrc,mlonsrc],h.f1);
            
            plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',Tilim,[mlatsrc,mlonsrc],h.f2);
            
            plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',Telim,[mlatsrc,mlonsrc],h.f3);
            
            plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',J1lim,[mlatsrc,mlonsrc],h.f4);
            
            plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],h.f5);
            
            plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],h.f6);
            
            plotfun(ymd,UTsec,xg,J2(:,:,:)*1e6,'J_2 (uA/m^2)',J2lim,[mlatsrc,mlonsrc],h.f7);
            
            plotfun(ymd,UTsec,xg,J3(:,:,:)*1e6,'J_3 (uA/m^2)',J3lim,[mlatsrc,mlonsrc],h.f8);
            
            if ~isempty(h.f9)
                h9a = axes('parent', h.f9);
                imagesc(Phitop, 'parent', h9a)
                colorbar;
            end
        end
        if saveplots   % for 3D or long 2D plots print and output file every time step
            dosave(flagoutput, direc, filename, h)
        end
    else    %short 2D simulation - put the entire time series in a single plot
        if it==1, disp('short 2D simulations...'), end
        ha = subplot(Nsp,Nsp,it,'parent',h.f10);
        nelim =  [9 11.3];
        plotfun(ymd,UTsec,xg,log10(ne), 'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],ha);
        
        if flagoutput~=3
            ha = subplot(Nsp,Nsp,it,'parent',h.f1);
            plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',v1lim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f2);
            plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',Tilim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f3);
            plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',Telim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f4);
            plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',J1lim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f5);
            plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f6);
            plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f7);
            plotfun(ymd,UTsec,xg,J2(:,:,:)*1e6,'J_2 (uA/m^2)',J2lim,[mlatsrc,mlonsrc],ha);
            
            ha = subplot(Nsp,Nsp,it,'parent',h.f8);
            plotfun(ymd,UTsec,xg,J3(:,:,:)*1e6,'J_3 (uA/m^2)',J3lim,[mlatsrc,mlonsrc],ha);
            
            if ~isempty(h.f9)
                ha = subplot(Nsp,Nsp,it,'parent',h.f9);
                imagesc(Phitop, 'parent', ha)
                colorbar;
            end
        end
    end
    
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end % for
    
    if saveplots && (xg.lx(2)==1 && xg.lx(3)==1 || lt<=16)    %save the short 2D sim plots
        dosave(flagoutput, direc, filename, h)
    end
    
if nargout==0, clear('xg'), end
    
end % function


function dosave(flagoutput, direc, filename, h)

if flagoutput~=3
    print(h.f1,'-dpng',[direc,'/v1plots/',filename,'.png'],'-r300')
    print(h.f2,'-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300')
    print(h.f3,'-dpng',[direc,'/Teplots/',filename,'.png'],'-r300')
    print(h.f4,'-dpng',[direc,'/J1plots/',filename,'.png'],'-r300')
    print(h.f5,'-dpng',[direc,'/v2plots/',filename,'.png'],'-r300')
    print(h.f6,'-dpng',[direc,'/v3plots/',filename,'.png'],'-r300')
    print(h.f7,'-dpng',[direc,'/J2plots/',filename,'.png'],'-r300')
    print(h.f8,'-dpng',[direc,'/J3plots/',filename,'.png'],'-r300')
    if ~isempty(h.f9)
        print(h.f9,'-dpng',[direc,'/Phiplots/',filename,'.png'],'-r300')
    end
end
print(h.f10,'-dpng',[direc,'/nplots/',filename,'.png'],'-r300')

if flagoutput~=3     %now make .eps prints of the plots
    print(h.f1,'-depsc2',[direc,'/v1plots/',filename,'.eps'])
    print(h.f2,'-depsc2',[direc,'/Tiplots/',filename,'.eps'])
    print(h.f3,'-depsc2',[direc,'/Teplots/',filename,'.eps'])
    print(h.f4,'-depsc2',[direc,'/J1plots/',filename,'.eps'])
    print(h.f5,'-depsc2',[direc,'/v2plots/',filename,'.eps'])
    print(h.f6,'-depsc2',[direc,'/v3plots/',filename,'.eps'])
    print(h.f7,'-depsc2',[direc,'/J2plots/',filename,'.eps'])
    print(h.f8,'-depsc2',[direc,'/J3plots/',filename,'.eps'])
    if ~isempty(h.f9)
        print(h.f9,'-depsc2',[direc,'/Phiplots/',filename,'.eps'])
    end
end
print(h.f10,'-depsc2',[direc,'/nplots/',filename,'.eps'])

end
