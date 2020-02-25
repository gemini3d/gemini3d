function particles_BCs_2d(p)
% create particle precip for 2D
narginchk(1,1)
validateattributes(p, {'struct'}, {'scalar'})

%% read desired simulation config
params = read_config(p.nml);
%% GRID ALREADY NEEDS TO BE MADE, AS WELL
xg = readgrid(p.simdir, p.format, p.realbits);
%% CREATE PRECIPITATION CHARACTERISTICS data

% number of grid cells. This will be interpolated to grid, so 100x100 is arbitrary
llon=100;
llat=100;
if xg.lx(2) == 1    % cartesian
  llon=1;
elseif xg.lx(3) == 1
  llat=1;
end
thetamin=min(xg.theta(:));
thetamax=max(xg.theta(:));
mlatmin=90-thetamax*180/pi;
mlatmax=90-thetamin*180/pi;
mlonmin=min(xg.phi(:))*180/pi;
mlonmax=max(xg.phi(:))*180/pi;
%mlat=linspace(mlatmin,mlatmax,llat);
%mlon=linspace(mlonmin,mlonmax,llon);
latbuf=1/100*(mlatmax-mlatmin);
lonbuf=1/100*(mlonmax-mlonmin);
mlat=linspace(mlatmin-latbuf,mlatmax+latbuf,llat);
mlon=linspace(mlonmin-lonbuf,mlonmax+lonbuf,llon);
[MLON,MLAT]=ndgrid(mlon,mlat);
mlonmean=mean(mlon);
mlatmean=mean(mlat);

%% disturbance width
mlatsig = p.precip_latwidth*(mlatmax-mlatmin);
mlatsig = max(mlatsig,0.01);    % can't let this go to zero...
mlonsig = p.precip_lonwidth*(mlonmax-mlonmin);

%% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=params.tdur;
%Nt=tdur+1;
%time=linspace(tmin,tmax,Nt)';
time=tmin:5:tmax;
Nt=numel(time);

%% COMPUTE THE TOTAL ENERGY FLUX AND CHAR. ENERGY
%{
Q=zeros(Nt,llon,llat);
E0=zeros(Nt,llon,llat);
for it=1:Nt
  Qtmp=squeeze(Qdat(it,:,:));
  Qtmp=interp2(glat,gloncorrected,Qtmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  Q(it,:,:)=reshape(Qtmp,[1 llon llat]);
  E0tmp=squeeze(E0dat(it,:,:));
  E0tmp=interp2(glat,gloncorrected,E0tmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  E0(it,:,:)=reshape(E0tmp,[1 llon llat]);
end
%}

%% time
ymd=params.ymd;
UTsec = params.UTsec0+time;     % seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd(:)',[Nt,1]),UThrs(:),zeros(Nt,1),zeros(Nt,1));
% t=datenum(expdate);

%{
%INTERPOLATE/DECIMATE TO 1 SECOND RESOLUTION
samplerate=1;    %sampling rate in seconds
startdate=[ymd0,UTsec0/3600,0,0];
startt=datenum(startdate);
%tmin=ceil(min(t)*86400)/86400;
tmin=ceil(startt*86400)/86400;
tmax=floor(max(t)*86400)/86400;
outputt=tmin:samplerate/86400:tmax;
outputdate=datevec(outputt);
ltout=numel(outputt);
Qit=zeros(llon,llat,ltout);
E0it=zeros(llon,llat,ltout);
for ilat=1:llat
   for ilon=1:llon
       Qhere=squeeze(Qsmooth(:,ilon,ilat));
       E0here=squeeze(E0smooth(:,ilon,ilat));
       Qi=interp1(t,Qhere,outputt);
       E0i=interp1(t,E0here,outputt);
       Qit(ilon,ilat,:)=reshape(Qi,[1 1 ltout]);
       E0it(ilon,ilat,:)=reshape(E0i,[1 1 ltout]);
   end
end
%}

%% CREATE PRECIPITATION INPUT DATA
% Qit: [mW m^-2]
% E0it: [eV]
Qit=zeros(llon,llat,Nt);
E0it=zeros(llon,llat,Nt);
for it=1:Nt
   Qit(:,:,it)=10*exp(-(MLON-mlonmean).^2/(2*mlonsig^2)).*exp(-(MLAT-mlatmean).^2/(2*mlatsig^2));
%  Qit(:,:,it)=5;
  E0it(:,:,it)=5e3;%*ones(llon,llat);
end

%% CONVERT THE ENERGY TO EV
%E0it=max(E0it,0.100);
%E0it=E0it*1e3;

%% SAVE to files
% LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
% FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
% THE EFIELD DATA DO NOT NEED TO BE SMOOTHED.

outdir = absolute_path([p.simdir, '/inputs/prec_inputs/']);
makedir(outdir)

switch p.format
  case {'h5','hdf5'}, write_hdf5(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it, p.realbits)
  case {'dat','raw'}, write_raw(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it, p.realbits)
  otherwise, error(['unknown file format ', p.format])
end

end % function


function write_hdf5(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it, realbits)
narginchk(10,10)

fn = [outdir, '/simsize.h5'];
disp(['write ', fn])
h5save(fn, '/llon', int32(llon))
h5save(fn, '/llat', int32(llat))

freal = ['float',int2str(realbits)];

fn = [outdir, '/simgrid.h5'];
disp(['write ', fn])
h5save(fn, '/mlon', mlon, [], freal)
h5save(fn, '/mlat', mlat, [], freal)

for i = 1:Nt
  UTsec = expdate(i,4)*3600 + expdate(i,5)*60 + expdate(i,6);
  ymd = expdate(i, 1:3);

  fn = [outdir, filesep, datelab(ymd,UTsec), '.h5'];
  disp(['writing ', fn])

  h5save(fn, '/Qp', Qit(:,:,i), [], freal)
  h5save(fn, '/E0p', E0it(:,:,i), [], freal)
end

end % function


function write_raw(outdir, llon, llat, mlon, mlat, expdate, Nt, Qit, E0it, realbits)
narginchk(10,10)

filename=[outdir, '/simsize.dat'];
disp(['write ', filename])
fid=fopen(filename, 'w');
fwrite(fid,llon,'integer*4');
fwrite(fid,llat,'integer*4');
fclose(fid);

freal = ['float', int2str(realbits)];

filename=[outdir, '/simgrid.dat'];
disp(['write ', filename])

fid=fopen(filename,'w');
fwrite(fid,mlon, freal);
fwrite(fid,mlat, freal);
fclose(fid);

for i = 1:Nt
  UTsec = expdate(i,4)*3600 + expdate(i,5)*60 + expdate(i,6);
  ymd = expdate(i, 1:3);

  filename = [outdir, filesep, datelab(ymd,UTsec), '.dat'];
  disp(['writing ', filename])

  fid = fopen(filename,'w');
  fwrite(fid,Qit(:,:,i), freal);
  fwrite(fid,E0it(:,:,i), freal);
  fclose(fid);
end

end % function
