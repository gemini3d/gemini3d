function xg = plotframe(direc,ymd,UTsec,saveplots,plotfun,xg,h,visible)


cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

%% CHECK ARGS AND SET SOME DEFAULT VALUES ON OPTIONAL ARGS
narginchk(3,8)
validateattributes(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)
validateattributes(ymd, {'numeric'}, {'vector','numel',3}, mfilename, 'year, month, day', 2)
validateattributes(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 3)

if nargin<4
  saveplots={};  % 'png', 'eps' or {'png', 'eps'}
end
if nargin<5
  plotfun=[];    %code will choose a plotting function if not given
end
if nargin<6
  xg=[];         %code will attempt to reload the grid
end

if nargin<8, visible = 'on'; end

Ncmap = parula(256);
Tcmap = parula(256);
Vcmap = bwr();
Jcmap = bwr();

lotsplots=true;   %@scivision may want to fix this...

direc = expanduser(direc);  % for mkdir on Octave.

%% NEED TO ADD SOME CODE TO VALIDATE THE INPUT DATA...


%% SET THE CAXIS LIMITS FOR THE PLOTS - really needs to be user provided somehow...
flagcaxlims=false;
%{
nelim =  [0 6e11];
v1lim = [-400 400];
Tilim = [100 3000];
Telim = [100 3000];
J1lim = [-25 25];
v2lim = [-10 10];
v3lim = [-10 10];
J2lim = [-10 10];
J3lim=[-10 10];
%}

%% MAKE DIRECTORIES FOR OUTPUT FILES
if ~isempty(saveplots)
  dlist = {'nplots', 'v1plots', 'v2plots', 'v3plots', 'J1plots',...
           'Tiplots', 'Teplots', 'J2plots', 'J3plots', 'Phiplots'};
  for i=1:length(dlist)
    mkdir([direc, filesep, dlist{i}]);
  end
end


%% READ IN THE SIMULATION INFORMATION (this is low cost so reread no matter what)
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc] = readconfig(direc, filesep, 'inputs');


%% CHECK WHETHER WE NEED TO RELOAD THE GRID (check if one is given because this can take a long time)
if isempty(xg)
  disp('Reloading grid...  Provide one as input if you do not want this to happen.')
  xg = readgrid([direc,filesep,'inputs',filesep]);
end

if nargin<7 || isempty(h)
  h = plotinit(xg, visible);
end

plotfun = grid2plotfun(plotfun, xg);

%% DETERMINE OUTPUT TIME NEAREST TO THAT REQUESTED
%tnow = UTsec0+tdur;    %not date corrected, but doesn't really matter
[ymdend,UTsecend]=dateinc(tdur,ymd0,UTsec0);


%% LOCATE TIME NEAREST TO THE REQUESTED DATE
%ymdprev=ymd0;
%UTsecprev=UTsec0;
ymdnext=ymd0;
UTsecnext=UTsec0;
it=1;
while(datenum([ymdnext,UTsecnext/3600,0,0])<datenum([ymd,UTsec/3600,0,0]) && datenum([ymdnext,UTsecnext/3600,0,0])<datenum([ymdend,UTsecend/3600,0,0]))
  ymdprev=ymdnext;
  UTsecprev=UTsecnext;
  [ymdnext,UTsecnext]=dateinc(dtout,ymdprev,UTsecprev);
  it=it+1;
end
%{
if(UTsecnext-UTsec<=UTsec-UTsecprev)   %we are closer to next frame
  ymdnow=ymdnext;
  UTsecnow=UTsecnext;
else    %closer to previous frame
  ymdnow=ymdprev;
  UTsecnow=UTsecprev;
end
%}

%% LOAD THE FRAME NEAREST TO THE REQUESTED TIME
[ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = loadframe(direc,ymd,UTsec,ymd0,UTsec0,tdur,dtout,flagoutput,mloc,xg);
disp([filename, ' => ', func2str(plotfun)])


%% UNTIL WE PROVIDE A WAY FOR THE USER TO SPECIFY COLOR AXES, JUST TRY TO SET THEM AUTOMATICALLY
if (~flagcaxlims)
 nelim =  [min(ne(:)), max(ne(:))];
 v1mod=max(abs(v1(:)));
 v1lim = [-v1mod, v1mod];
 Tilim = [0, max(Ti(:))];
 Telim = [0, max(Te(:))];
 J1mod=max(abs(J1(:)));
 J1lim = [-J1mod, J1mod];
 v2mod=max(abs(v2(:)));
 v2lim = [-v2mod, v2mod];
 v3mod=max(abs(v3(:)));
 v3lim = [-v3mod, v3mod];
 J2mod=max(abs(J2(:)));
 J2lim = [-J2mod, J2mod];
 J3mod=max(abs(J3(:)));
 J3lim=[-J3mod, J3mod];
%   nelim =  [min(ne(:)), max(ne(:))];
%   v1lim = [NaN,NaN];
%   Tilim = [0, max(Ti(:))];
%   Telim = [0, max(Te(:))];
%   J1lim = [NaN,NaN];
%   v2lim = [NaN,NaN];
%   v3lim = [NaN,NaN];
%   J2lim = [NaN,NaN];
%   J3lim=[NaN,NaN];
end


%% MAKE THE PLOTS (WHERE APPROPRIATE)
%Electron number density, 'position', [.1, .1, .5, .5], 'units', 'normalized'
if lotsplots   % 3D simulation or a very long 2D simulation - do separate plots for each time frame

    clf(h.f10)
    plotfun(ymd,UTsec,xg, ne, 'n_e (m^{-3})', nelim,[mlatsrc,mlonsrc], h.f10, Ncmap);

    if flagoutput~=3
        clf(h.f1)
        plotfun(ymd,UTsec,xg,v1, 'v_1 (m/s)', v1lim,[mlatsrc,mlonsrc], h.f1, Vcmap);
        clf(h.f2)
        plotfun(ymd,UTsec,xg,Ti,'T_i (K)',Tilim,[mlatsrc,mlonsrc], h.f2, Tcmap);
        clf(h.f3)
        plotfun(ymd,UTsec,xg,Te,'T_e (K)',Telim,[mlatsrc,mlonsrc], h.f3, Tcmap);
        clf(h.f4)
        plotfun(ymd,UTsec,xg,J1,'J_1 (A/m^2)',J1lim,[mlatsrc,mlonsrc],h.f4, Jcmap);
        clf(h.f5)
        plotfun(ymd,UTsec,xg,v2,'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],h.f5, Vcmap);
        clf(h.f6)
        plotfun(ymd,UTsec,xg,v3,'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],h.f6, Vcmap);
        clf(h.f7)
        plotfun(ymd,UTsec,xg,J2,'J_2 (A/m^2)',J2lim,[mlatsrc,mlonsrc],h.f7, Jcmap);
        clf(h.f8)
        plotfun(ymd,UTsec,xg,J3,'J_3 (A/m^2)',J3lim,[mlatsrc,mlonsrc],h.f8, Jcmap);

        if ~isempty(h.f9)
            clf(h.f9)
            h9a = axes('parent', h.f9);
            imagesc(Phitop, 'parent', h9a)
            colorbar('peer', h9a)
        end
    end

    % for 3D or long 2D plots print and output file every time step
    dosave(flagoutput, direc, filename, saveplots, h)

else    %short 2D simulation - put the entire time series in a single plot

    figure(h.f10)
    ha = subplot(Rsp, Csp,it,'parent',h.f10);
    nelim =  [9 11.3];
    plotfun(ymd,UTsec,xg,log10(ne), 'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],ha);

    if flagoutput~=3
        ha = subplot(Rsp, Csp,it,'parent',h.f1);
        plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',v1lim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f2);
        plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',Tilim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f3);
        plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',Telim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f4);
        plotfun(ymd,UTsec,xg,J1(:,:,:),'J_1 (A/m^2)',J1lim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f5);
        plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',v2lim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f6);
        plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',v3lim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f7);
        plotfun(ymd,UTsec,xg,J2(:,:,:),'J_2 (A/m^2)',J2lim,[mlatsrc,mlonsrc],ha);

        ha = subplot(Rsp, Csp,it,'parent',h.f8);
        plotfun(ymd,UTsec,xg,J3(:,:,:),'J_3 (A/m^2)',J3lim,[mlatsrc,mlonsrc],ha);

        if ~isempty(h.f9)

            ha = subplot(Rsp, Csp,it,'parent',h.f9);
            imagesc(Phitop, 'parent', ha)
            colorbar('peer', ha)
        end
    end

end

if ~lotsplots    %save the short 2D sim plots
    dosave(flagoutput, direc, filename, saveplots, h)
end

%% Don't print
if nargout==0, clear('xg'), end

end % function plotframe



%%FUNCTION THAT CREATES IMAGE FILES FROM PLOTS
function dosave(flagoutput, direc, filename, fmt, h)

narginchk(5,5)
validateattr(flagoutput, {'numeric'}, {'scalar'}, mfilename)
validateattr(direc, {'char'}, {'vector'}, mfilename)
validateattr(filename, {'char'}, {'vector'}, mfilename)

if any(strcmpi(fmt, 'png'))
    disp(['writing png plots to ',direc])

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
end

if any(strcmpi(fmt, 'eps'))
    disp(['writing eps plots to ',direc])

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

end %function dosave
