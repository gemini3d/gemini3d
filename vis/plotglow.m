function plotglow(direc)
narginchk(1,1)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

assert(is_folder(direc), [direc, ' is not a directory.'])

%array of volume emission rates at each altitude; cm-3 s-1:
%           3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A,
%           7774A, 8446A, 3726A, LBH, 1356, 1493, 1304; cm-3 s-1





%READ IN SIMULATION INFO
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc, filesep, 'inputs']);


%READ IN THE GRID
xg = readgrid([direc,filesep,'inputs']);


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
file_list=dir([direc,filesep,'*.dat']);
lt=numel(file_list);
bFrame = zeros(lx2,lx3,lwave,lt);
for ind = 1:lt
    figure
    filename = [direc,filesep,'aurmaps/',file_list(ind).name];
    bFrame(:,:,:,ind) = loadglow_aurmap(filename, lx2, lx3, lwave);
    
    if lx3 > 1
        plot1wl(x2, x3, bFrame(:,:,:,ind))
    else
        imagesc(squeeze(bFrame(:,:,:,ind)))    %there's a singleton dimension here since this corresponds to a 2D simulation
        colorbar
    end
end

end % function


function plot1wl(x2, x3, bFrame)

h=imagesc(x2,x3,log10(squeeze(bFrame(:,:,2)))');
axis xy;
axis tight;
%caxis(caxlims);
cb=colorbar;
cby=get(cb,'ytick');
%set(cb,'yticklabel',sprintf('10^{%g}|',cby))
ylabel(cb,'427.8 nm Intensity (R)')
xlabel('Eastward Distance (km)')
ylabel('Northward Distance (km)')
title('427.8 nm intensity')

end

%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libxvid -r 30 -q:v 0 isinglass_geminiglow_4278.avi
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
