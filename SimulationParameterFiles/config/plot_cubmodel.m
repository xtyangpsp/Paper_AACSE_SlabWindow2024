define_latlon;

[outlat, outlon, outdep, outvs]=read_CUB_revised(minlon,maxlon,minlat,maxlat);
%%
%21 -> 80 km; 40 -> 156 km; 76-> 300km
zidx=5;
figure;pcolor(squeeze(outlon)-360,squeeze(outlat),squeeze(outvs(:,:,zidx))); 
shading interp; colorbar;
colormap(flip(colormap('jet')));