direc='~/simulations/ICI2_glow0/'
flagplot=0;    %whether or not to create plots


%array of volume emission rates at each altitude; cm-3 s-1:
%           3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A,
%           7774A, 8446A, 3726A, LBH, 1356, 1493, 1304; cm-3 s-1


%SET PATHS
cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])


%READ IN SIMULATION INFO
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs/config.ini']);


%READ IN THE GRID
disp('Reloading grid...  Provide one as input if you do not want this to happen.')
xg = readgrid([direc,filesep,'inputs',filesep]);


%GET THE SYSTEM SIZE
lwave=15;
lx2=xg.lx(2);
lx3=xg.lx(3);
x2=xg.x2(3:end-2);
x3=xg.x3(3:end-2);


%PLOT PARAMETERS
FS=14;
caxlims=[0.0,5.0];


%LOOP OVER FRAMES, SHOULDN'T BE A PROBLEM TO LOAD ALL AT ONCE
file_list=dir([direc,'*.dat']);
lt=numel(file_list);
bFrame = zeros(lx2,lx3,lwave,lt);
for ind = 1:lt
    filename=file_list(ind).name;
	fid=fopen([direc,filesep,'aurmaps/',filename],'r');
	cAur = fread(fid,lx2*lx3*lwave,'real*8');
	cAur = reshape(cAur,[lx2,lx3,lwave]);
	bFrame(:,:,:,ind)=cAur(:,:,:);
    fclose(fid);
	%max(max(bFrame))
    
    if (flagplot)
        h=imagesc(x2,x3,log10(bFrame)');
        axis xy;
        axis tight;
        caxis(caxlims);
        cb=colorbar;
        cby=get(cb,'ytick');
        set(cb,'yticklabel',sprintf('10^{%g}|',cby),'fontsize',FS);
        ylabel(cb,'427.8 nm Intensity (R)','fontsize',FS);
        xlabel('Eastward Distance (km)');
        ylabel('Northward Distance (km)');
        title(file_list{ind});
        set(gca,'fontsize',FS)%,'fontweight','bold');
        fname=[file_list{ind} '.png'];
        print('-dpng', fname);
        clf
    end
end
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libxvid -r 30 -q:v 0 isinglass_geminiglow_4278.avi
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
