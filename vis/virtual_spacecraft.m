%FLIES A VIRTUAL SPACECRAFT THROUGH THE MODEL
% only works with a cartesian grid for now.


%INPUT DATA
direc='~/zettergmdata/simulations/Aether_discrete'

addpath ../script_utils;


%%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs', filesep, 'config.ini']);
ymd=ymd0; UTsec=UTsec0;


%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if ~exist('xg','var')
  fprintf('Loading grid...\n');
  xg = readgrid([direc,filesep,'inputs',filesep]);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  x1=xg.x1(3:end-2); x2=xg.x2(3:end-2); x3=xg.x3(3:end-2);
  [X2,X1,X3]=meshgrid(x2,x1,x3);

  %In case the simulation is cartesian we need to know the center geomgagetic coordinates of the grid
  thetactr=mean(xg.theta(:));
  phictr=mean(xg.phi(:));
  [glatctr,glonctr]=geomag2geog(thetactr,phictr);
end


%% COMPUTE SOURUCE LOCATION IN MCOORDS
if ~isempty(mloc)
 mlatctr=mloc(1);
 mlonctr=mloc(2);
else
 mlatctr=[];
 mlonctr=[];
end


%TIMES WHERE WE HAVE MODEL OUTPUT
times=UTsec0:dtout:UTsec0+tdur;
lt=numel(times);
datemod0=datenum([ymd0,UTsec0/3600,0,0]);
datemod=datemod0:dtout/86400:datemod0+tdur/86400;


%DEFINE SOME SORT OF SATELLITE ORB.
lorb=250;    %number of times series for satellite orbit
UTsat=linspace(min(times),max(times),lorb);
ymdsat=repmat(ymd0,[lorb,1]);
datevecsat=[ymdsat,UTsat(:)/3600,zeros(lorb,1),zeros(lorb,1)];
datesat=datenum(datevecsat);

thetasat=linspace(min(xg.theta(:)),max(xg.theta(:)),lorb);
phisat=linspace(min(xg.phi(:)),max(xg.phi(:)),lorb);
[glatsat,glonsat]=geomag2geog(thetasat,phisat);
%glatsat=linspace(min(xg.glat(:)),max(xg.glat(:)),lorb);
%glonsat=linspace(min(xg.glon(:)),max(xg.glon(:)),lorb);
altsat=linspace(100e3,600e3,lorb);


%MAIN LOOP OVER ORBIT SEGMENTS
datebufprev=datemod0;
datebufnext=datemod0;
neprev=[]; nenext=[];
viprev=[]; vinext=[];
Tiprev=[]; Tinext=[];
Teprev=[]; Tenext=[];
nesat=zeros(1,lorb);
for iorb=1:lorb
  datenow=datesat(iorb);


  %FIND THE TWO FRAMES THAT BRACKET THIS ORBIT TIME
  datemodnext=datemod0;
  datemodprev=datemod0;
  while(datemodnext<datenow & datemodnext<=datemod(end))
    datemodprev=datemodnext;
    datemodnext=datemodnext+dtout/86400;    %matlab datenums are in units of days from 0000
  end
  datestr(datemodprev),datestr(datenow),datestr(datemodnext)


  if (datemodnext==datemodprev | datemodnext>datemod(end))    %set everything to NaNs if outside model time domain
    fprintf('Requested time is out of bounds...\n');
    neprev=zeros(lx1,lx2,lx3); nenext=neprev;
    viprev=neprev; vinext=neprev;
    Tiprev=neprev; Tinext=neprev;
    Teprev=neprev; Tenext=neprev;
    datestr(datenow)
  else     %go ahead and read in data and set up the interpolations
    %DATA BUFFER UPDATES
    if (datebufprev~=datemodprev | isempty(neprev))    %need to reload the previous output frame data buffers
      fprintf('Loading previous buffer...\n');
      datevecmodprev=datevec(datemodprev);
      ymd=datevecmodprev(1:3);
      UTsec=datevecmodprev(4)*3600+datevecmodprev(5)*60+datevecmodprev(6);
      UTsec=round(UTsec);    %some accuracy problems...  this is fishy...
      [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop]=loadframe(direc,UTsec,ymd,UTsec0,ymd0,mloc,xg);
      neprev=ne; viprev=v1; Tiprev=Ti; Teprev=Te;
      clear ne v1 Ti Te J1 v2 v3 J2 J3 Phitop;    %avoid keeping extra copies of large data
      datebufprev=datemodprev; 
%      fprintf('Creating previous interpolant...\n');
%      neintpprev=scatteredInterpolant(xg.alt(:),xg.glon(:),xg.glat(:),neprev(:));    %way too slow
    end
    if (datebufnext~=datemodnext | isempty(nenext))    %need to reload the next output frame data buffers
      fprintf('Loading next buffer...\n');
      datevecmodnext=datevec(datemodnext);
      ymd=datevecmodnext(1:3);
      UTsec=datevecmodnext(4)*3600+datevecmodnext(5)*60+datevecmodnext(6);
      UTsec=round(UTsec);
      [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop]=loadframe(direc,UTsec,ymd,UTsec0,ymd0,mloc,xg);
      nenext=ne; vinext=v1; Tinext=Ti; Tenext=Te;
      clear ne v1 Ti T3 J1 v2 v3 J2 J3 Phitop;    %avoid keeping extra copies of large data
      datebufnext=datemodnext;
%      fprintf('Creating next interpolant...\n');
%      neintpnext=scatteredInterpolant(xg.alt(:),xg.glon(:),xg.glat(:),nenext(:));    %way too slow
    end
  end

    %INTERPOLATIONS
%    nesatprev=neintpprev(altsat(iorb),glonsat(iorb),glatsat(iorb));
%    nesatnext=neintpnext(altsat(iorb),glonsat(iorb),glatsat(iorb));
    [x1sat,x2sat,x3sat]=geog2UEN(altsat(iorb),glonsat(iorb),glatsat(iorb),thetactr,phictr);
    fprintf('Starting interpolations...\n')
    nesatprev=interp3(X2,X1,X3,neprev,x2sat,x1sat,x3sat)
    nesatnext=interp3(X2,X1,X3,nenext,x2sat,x1sat,x3sat)
    nesat(iorb)=nesatprev+(nesatnext-nesatprev)/(datemodnext-datemodprev)*(datenow-datemodprev);
end

rmpath ../script_utils;

