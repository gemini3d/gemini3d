function xg = plotframe(direc,ymd,UTsec,saveplot_fmt,plotfun,xg,h,visible)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])

%% CHECK ARGS AND SET SOME DEFAULT VALUES ON OPTIONAL ARGS
narginchk(3,8)
validateattributes(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)
validateattributes(ymd, {'numeric'}, {'vector','numel',3}, mfilename, 'year, month, day', 2)
validateattributes(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 3)

if nargin<4
  saveplot_fmt={};  % {'png'} or {'png', 'eps'}
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
Phi_cmap = bwr();
Vcmap = bwr();
Jcmap = bwr();

lotsplots=true;   %@scivision may want to fix this...

direc = expanduser(direc);  % for mkdir on Octave.

% FIXME: VALIDATE THE INPUT DATA...

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

%% READ IN THE SIMULATION INFORMATION (this is low cost so reread no matter what)
params = read_config([direc, filesep, 'inputs']);

%% CHECK WHETHER WE NEED TO RELOAD THE GRID (check if one is given because this can take a long time)
if isempty(xg)
  disp('Reloading grid...  Provide one as input if you do not want this to happen.')
  xg = readgrid([direc,filesep,'inputs']);
end

if nargin<7 || isempty(h)
  h = plotinit(xg, visible);
end

plotfun = grid2plotfun(plotfun, xg);

%% DETERMINE OUTPUT TIME NEAREST TO THAT REQUESTED
%tnow = UTsec0+tdur;    %not date corrected, but doesn't really matter
[ymdend,UTsecend] = dateinc(params.tdur, params.ymd, params.UTsec0);

%% LOCATE TIME NEAREST TO THE REQUESTED DATE
%ymdprev=ymd0;
%UTsecprev=UTsec0;
ymdnext = params.ymd(:)';
UTsecnext = params.UTsec0;
it=1;
while(datenum([ymdnext,UTsecnext/3600,0,0]) < datenum([ymd,UTsec/3600,0,0]) && datenum([ymdnext,UTsecnext/3600,0,0]) < datenum([ymdend,UTsecend/3600,0,0]))
  ymdprev=ymdnext;
  UTsecprev=UTsecnext;
  [ymdnext,UTsecnext]=dateinc(params.dtout,ymdprev,UTsecprev);
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
[ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop] = loadframe(direc,ymd,UTsec, params.flagoutput, params.mloc,xg);
disp([filename, ' => ', func2str(plotfun)])

%% UNTIL WE PROVIDE A WAY FOR THE USER TO SPECIFY COLOR AXES, JUST TRY TO SET THEM AUTOMATICALLY
if (~flagcaxlims)
 nelim =  [min(ne(:)), max(ne(:))];
% v1mod=max(abs(v1(:)));
 v1mod=80;
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
 Phitop_lim = [min(Phitop(:)), max(Phitop(:))];
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

  if params.flagoutput ~= 3
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
    clf(h.f9)
%    % TODO: check units
%    plotfun(ymd,UTsec,xg,Phitop,'topside potential \Phi_{top} (V)', Phitop_lim, [mlatsrc, mlonsrc], h.f9, Phi_cmap)
  end

else    %short 2D simulation - put the entire time series in a single plot

  figure(h.f10) %#ok<*UNRCH>
  ha = subplot(Rsp, Csp, it, 'parent',h.f10);
  nelim =  [9 11.3];
  plotfun(ymd,UTsec,xg,log10(ne), 'log_{10} n_e (m^{-3})',nelim,[mlatsrc,mlonsrc],ha);

  if params.flagoutput ~= 3
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

    ha = subplot(Rsp, Csp,it,'parent',h.f9);
%    % TODO: check units
%    plotfun(ymd,UTsec,xg,Phitop,'topside potential \Phi_{top} (V)', Phitop_lim, [mlatsrc, mlonsrc], h.f9, Phi_cmap)
  end

end


if lotsplots
  % for 3D or long 2D plots print and output file every time step
  saveframe(params.flagoutput, direc, filename, saveplot_fmt, h)
end

%% Don't print
if nargout==0, clear('xg'), end

end % function plotframe

