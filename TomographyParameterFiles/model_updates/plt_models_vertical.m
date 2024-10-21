%
clear all; close all;
% /opt/matlab/2007b/bin/matlab
MFILE_ROOT='../mfiles';
CODES_ROOT='../Codes';
NC_ROOT='../SeismicDataProcessing';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);
addpath([MFILE_ROOT]);
addpath([CODES_ROOT]);
addpath([NC_ROOT]);
addpath('~/DEPOTVINCE/AACSE_tomo_SURF/SeismicDataProcessing/MyMatlabToolboxPublic/mexcdf/snctools')

ite_nm = ['ite_0.05deg_01'];
fnm_conf=['./SeisFD3D.conf_' ite_nm];

iplot = 0; % 0, absolute velocity; 1, velocity perturbation
% 
% rd=[(0:31)/31,ones(1,32)];
% gn=[(0:31)/31,(31:-1:0)/31];
% bl=[ones(1,32),(31:-1:0)/31];
% rwb=[rd',gn',bl'];
% rwb=flipud(rwb);

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm];
dir_coord=['./input_' ite_nm];

disp(['Read model... ',dir_media]);

id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=1;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
[XSIM{n},YSIM{n},ZSIM{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
% convert from radian to degrees
XSIM{n}=90-XSIM{n}*180/pi; %latitude
YSIM{n}=YSIM{n}*180/pi;
%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM{1}(end-npml,1,end);maxlat=XSIM{1}(1+npml,1,end);
minlon=YSIM{1}(1,1+npml,end);maxlon=YSIM{1}(1,end-npml,end);

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;
mvs{n}=smooth3(mvs{n},'box',[5 5 1]);

% load Topography, ../misc/CascadiaTopo.dat
% origgrid=0.05;
% [lon,lat,elev]=textread('../misc/CascadiaTopo.dat','%f %f %f');
% Topo=-elev/1000;
%%
maxdepth=100;
xlabelunit='deg'; %deg (for degrees) or km (for km distance)
xlabeldegtype='auto'; %auto (automatically decides which of lat or lon is used;
                      %lon: use lontitude
                      %lat: use latitude.
nmarkpt=10; %degree to km coefficient.
plotmoho=0;
mohofile='depth_mat_cubiclinterp.txt';

% set starting and ending points for plotting vertical profile images
% #-75.76 46.98
% #-68.66 44.56
% #<
% #-75.95 45.83
% #-70.87 43.02
% #<
% #-75.8363 43.00
% #-70.78   41.90
% #<
% -70.047 47.25843
% -74.9728 39.95628
ptLon1 = [198.43968, 202.441051, 195.579244, 203.973767];
ptLon2 = [212.731443, 215.369599, 201.603772, 209.653452];
ptLat1 = [57.088515, 59.467408, 56.218923, 58.470721];
ptLat2 = [51.917168, 53.904338, 52.133488, 56.511018];
subplotpar=[3,1];
plongridnum = length(ptLon1);
platgridnum = length(ptLat1);
figlabel={'A-A'' ','B-B'' ','C-C'' ','D-D'' '};
vert_exaggeration=2;

for np=1:plongridnum,
% start point (slat slon)  and end point (elat elon) for plotting vertical image profile
    slat = ptLat1(np); 
	slon = ptLon1(np); 
    elat = ptLat2(np); 
	elon = ptLon2(np);
				    
	proflength = deg2km(distance(slat, slon, elat, elon));
    %nmarkpt = round(proflength/deg2kmcoe);
	markptdeltadist = proflength/(nmarkpt-1);
	%markptdist = zeros(1, nmarkpt);
    markptdist = 0:markptdeltadist:proflength;
% 	markptlon = zeros(1, nmarkpt);
% 	markptlat = zeros(1, nmarkpt);
												    
	% number of pts on this profile
	profptnum = 5*round((proflength/5));
	% final pts (flat flon) that are shown on the profile
	flat = slat + (0:(profptnum-1))*(elat - slat)/(profptnum - 1);
	flon = slon + (0:(profptnum-1))*(elon - slon)/(profptnum - 1);
																		       
	dist = zeros(1, profptnum);
    for i = 1:(profptnum-1)
        ptdist = deg2km(distance(flat(i), flon(i), flat(i+1), flon(i+1)));
        dist(i+1) = ptdist + dist(i);
    end
	markptlon = interp1(dist, flon, markptdist, 'pchip');
	markptlat = interp1(dist, flat, markptdist, 'pchip');
																												       
	intpdist = 0:(proflength/(profptnum-1)):proflength;
	proflength = dist(profptnum);
	depth=6371-abs(ZSIM{1}(npml,npml,:)/1000);
	[distmesh,depthmesh] = meshgrid(dist, depth);
    idepth=5:1:100;
	[intpdistmesh,intpdepthmesh] = meshgrid(intpdist, idepth);
																																	       
    intpvVert = nan(size(mvs{1},3),profptnum);

	for k = 1:size(mvs{1},3)
        v=squeeze(mvs{1}(:,:,k));
        intpvVert(k, 1:profptnum) = interp2(YSIM{1}(1,:,end), XSIM{1}(:,1,end), ...
            v, flon, flat, 'spline')/1000; 
    end
    
    intpvVert2 = interp2(distmesh,depthmesh, intpvVert, intpdistmesh,intpdepthmesh, 'spline');
    figure('Position',[400 400 200*max(intpdist)/(vert_exaggeration*maxdepth) 200]), hold on, box on, axis on
    %
    if plotmoho
        clear mx my mz md;
        [mx,my,mz,md]=surface2profile(mohofile,slon-360,slat,elon-360,elat,0.05,0.05,200,'yes','other'); 
%         mx=360+mx;
    end
    if strcmp(xlabelunit,'km')
        imagesc(intpdist, idepth, intpvVert2);
        axis([0 ceil(max(intpdist)) 0 maxdepth]);
        daspect([1 1/vert_exaggeration 1]);
        if plotmoho
           plot(md,mz,'w-','linewidth',2); 
        end
    elseif strcmp(xlabelunit,'deg')
        if strcmp(xlabeldegtype,'auto')
            if range(flon) >= range(flat)
                imagesc(flon-360, idepth, intpvVert2);
                axis([min(flon-360) max(flon-360) 0 maxdepth]);
                daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
                if plotmoho
                   plot(mx,mz,'w-','linewidth',2); 
                end
            else
                imagesc(flat, idepth, intpvVert2);
                axis([min(flat) max(flat) 0 maxdepth]);
                daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
                if plotmoho
                   plot(my,mz,'w-','linewidth',2); 
                end
            end
        elseif strcmp(xlabeldegtype,'lon')
            imagesc(flon-360, idepth, intpvVert2);
            axis([min(flon-360) max(flon-360) 0 maxdepth]);
            daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
            if plotmoho
               plot(mx,mz,'w-','linewidth',2); 
            end
        elseif strcmp(xlabeldegtype,'lat')
            imagesc(flat, idepth, intpvVert2);
            axis([min(flat) max(flat) 0 maxdepth]);
            daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
            if plotmoho
               plot(my,mz,'w-','linewidth',2); 
            end
        else
            error('***ERROR: wrong value for [xlabeldegtype]');
        end
    else
        error('***ERROR: wrong value for [xlabelunit]');
    end

    title(strcat(figlabel{np}, ' (',num2str(slon-360),', ',num2str(slat),') to (',num2str(elon-360),...
        ', ',num2str(elat),')'),'FontSize',14);
    
    %set(gca,'XTick',[0:100:max(intpdist)],'XTickLabel',[0:100:max(intpdist)]);
    set(gca,'YTick',0:20:maxdepth,'YTickLabel',0:20:maxdepth);
	set(gca,'TickDir','out');
    set(gca,'FontSize',14);
    set(gca,'YDir','reverse');
    %xlabel('Longitude (deg.)','FontSize',17);
    if strcmp(xlabelunit,'km')
        xlabel('Distance (km)','FontSize',14);
    elseif strcmp(xlabelunit,'deg')
        xlabel('Degree','FontSize',14);
    else
        error('***ERROR: wrong value for [xlabelunit]');
    end
    ylabel('Depth (km)','FontSize',14);
    shading interp;
    colormap('jetwr'); 
    %colormap(rwb);
    if iplot==0
    %caxis([3.2 5.0]);
    caxis([3.7 4.7]);
    elseif iplot==1
    caxis([-8 8])
    end
    colorbar;
%     set(h,'fontsize',13);
   
 
%%%% save figure
% set(gcf,'PaperPositionMode','auto');   
% switch iplot
%     case 0
%         figname = ['VelModel_EW_profile_OBS_' num2str(ptLat1(np)) '.eps']
%     case 1
%         figname = ['VelModel_Percent_EW_profile_OBS_' num2str(ptLat1(np)) '.eps']
% end
% eval(['print -depsc ' figname])
%pause
% close all
    %drawnow;
    hold off
    
%     if np==1
%         pos=get(gca,'Position');
%         pos4=pos(4);
%     end
%     pos=get(gca,'Position');
%     pos(4)=pos4;
%     set(gca,'Position',pos);
    saveas(gca,strcat(ite_nm,'-',num2str(slon),'_',num2str(slat),'_to_',num2str(elon),...
        '_',num2str(elat),'.png'),'png');
    
end

