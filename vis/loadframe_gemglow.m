flagplot=0;


%GET THE SYSTEM SIZE
lwave=15;
lx2=128;
lx3=128;


%PLOT PARAMETERS
FS=14;
colormap("ocean");
caxlims=[0.0,5.0];


%LOOP OVER FRAMES, SHOULDN'T BE A PROBLEM TO LOAD ALL AT ONCE
file_list=glob('*.dat');
bFrame = zeros(lx2,lx3);
for ind = 1:size(file_list,1)
	fid=fopen(file_list{ind},'r');
	cAur = fread(fid,lx2*lx3*lwave,'real*8');
	cAur = reshape(cAur,[lx2,lx3,lwave]);
	bFrame(:,:)=cAur(:,:,2);
	%max(max(bFrame))
    
    if (flagplot)
        h=imagesc([-220,220],[-220,220],transpose(log10(bFrame)),climit=caxlims);
        axis xy;
        axis tight;
        caxis(caxlims);
        cb=colorbar;
        cby=get(cb,"ytick");
        set(cb,"yticklabel",sprintf("10^{%g}|",cby),"fontsize",FS);
        ylabel(cb,"427.8 nm Intensity (R)","fontsize",FS);
        xlabel("Eastward Distance (km)");
        ylabel("Northward Distance (km)");
        title(file_list{ind});
        set(gca,"fontsize",FS)%,"fontweight","bold");
        fname=[file_list{ind} ".png"];
        print("-dpng", fname);
        clf
    end if
endfor
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libxvid -r 30 -q:v 0 isinglass_geminiglow_4278.avi
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
