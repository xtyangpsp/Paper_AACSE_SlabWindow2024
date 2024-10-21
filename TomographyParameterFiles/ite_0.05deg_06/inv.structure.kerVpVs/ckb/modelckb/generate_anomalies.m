%generate anomaly bodies
clear all;
mx=46;
my=79;
mz=43;
fnm_blk=strcat('../../block.',num2str(mx),'x',num2str(my),'x',num2str(mz),'.1x1x1.1x1x1.nc');

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
xc=xc/1e3; yc=yc/1e3; zc=zc/1e3; %in km

minlon=min(y(1,:,1))-360; maxlon=max(y(1,:,1))-360; minlat=min(x(:,1,1));maxlat=max(x(:,1,1));

%%
outmodlefile='layer40km.txt';
vout=inv_make_anomaly(outmodlefile,mx,my,mz,...
    [1 mx 1 my mz-9 mz -0.1],...
    [1 mx 1 my mz-13 mz-10 0.1],...
    [1 mx 1 my mz-16 mz-14 -0.1],...
    [1 mx 1 my mz-19 mz-17 0.1],...
    [1 mx 1 my mz-22 mz-20 -0.1],...
    [1 mx 1 my mz-25 mz-23 0.1]);
%     [1 mx 1 my mz-19 mz-18 -0.1],...
%     [1 mx 1 my mz-21 mz-20 0.1],...
%     [1 mx 1 my mz-23 mz-22 -0.1]); % 5 layer
%     [1 mx 1 my mz-8 mz -0.1],...
%     [1 mx 1 my mz-11 mz-9 0.1],...
%     [1 mx 1 my mz-14 mz-12 -0.1],...
%     [1 mx 1 my mz-17 mz-15 0.1],...
%     [1 mx 1 my mz-20 mz-18 -0.1],...
%     [1 mx 1 my mz-23 mz-21 0.1]); % 5 layer

%
depmax=200;
figure('Position',[400 400 700 650]);
vout=load(outmodlefile);
vck=reshape(vout,[mx my mz]);
nx=25; ny=25;nz=[41,32]; %41,33
clear vplot;
vplot=squeeze(vck(nx,:,:)*100);
subplot(3,2,[1 2])
pcolor(squeeze(y(1,:,1))'-360,squeeze(z(1,1,:))',vplot');
colormap('jetwr');
set(gca,'CLim',[-10 10])

set(gca,'YDir','reverse')
colorbar;
ylim([0 depmax]);
grid on;
xlim([minlon maxlon])
title(['latitude: ' num2str(x(nx,1,1))]);

clear vplot;
vplot=squeeze(vck(:,ny,:)*100);
subplot(3,2,[3 4])
pcolor(squeeze(x(:,1,1))',squeeze(z(1,1,:))',vplot');
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','reverse')
colorbar;
ylim([0 depmax]);
grid on;
xlim([minlat maxlat])
title(['longitude: ' num2str(y(1,ny,1)-360)]);
%
clear vplot;
vplot=squeeze(vck(:,:,nz(1))*100);
subplot(3,2,5); 
pcolor(squeeze(y(1,:,1))'-360,squeeze(x(:,1,1))',vplot);
hold on;
% for sb=1:length(state)
%     plot(state(sb).polygon(:,1)*0.997-0.2, state(sb).polygon(:,2)*0.993+0.25,'color',[.3 .3 .3],'LineWidth',1);
% end
hold off;
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','normal')
colorbar;
grid on;
%xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);
daspect([1,cosd((minlat+maxlat)/2),1]);

title(['depth: ' num2str(z(1,1,nz(1)))]);
shading flat;

clear vplot;
% nz=36;
vplot=squeeze(vck(:,:,nz(2))*100);
subplot(3,2,6)
pcolor(squeeze(y(1,:,1))'-360,squeeze(x(:,1,1))',vplot);
hold on;
% for sb=1:length(state)
%     plot(state(sb).polygon(:,1)*0.997-0.2, state(sb).polygon(:,2)*0.993+0.25,'color',[.3 .3 .3],'LineWidth',1);
% end
hold off;
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','normal')
colorbar;
grid on;
axis([minlon maxlon minlat maxlat]);
daspect([1,cosd((minlat+maxlat)/2),1]);
title(['depth: ' num2str(z(1,1,nz(2)))]);
shading flat;