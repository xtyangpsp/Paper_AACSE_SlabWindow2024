function plt_raypath(delayfile,stationfile,fband,ite_tag,ray_min,cmax)
% plot paths coverage at certain period
% delayfile: phase delay file. This needs to be the one after running
%           get_delay_latlon().
% stationfile: text file containing station information
% fband: frequency band [f1min, f1max; f2min, f2max]
% ite_tag: iteration tag.
% ray_min: minimum ray counts to plot and to make mask.
% cmax: maximum colorbar value.

% close all; %clear all;
% ite_tag = 'ite_0.05deg_01';
% delayfile=[ite_nm,'_measure_result_craton_flat15.dat_latlon'];
% PROJHOME = '/Users/xtyang/Work/Research/Projects/Craton';
% stationfile = [PROJHOME '/STinfo/craton_station_withdata.txt'];
[stationlon,stationlat]=textread(stationfile,'%*s %*s %f %f %*f\n');
stationlon=stationlon-360;
% pband=[7.5 15; 10 20; 15 30; 20 40; 30 60; 40 75; 60 100; 75 125];
% fband=flip(flip(1./pband),2);
%%
% set plotting area [lon1 lon2] x [lat1 lat2]
minlon=floor(min(stationlon)); maxlon=ceil(max(stationlon)); 
minlat=floor(min(stationlat));maxlat=ceil(max(stationlat));
magextend=0.5; %extend the map area from the station coverage.
lon1 = minlon+360; 
lon2 = maxlon+360;  
lat1 = minlat;  
lat2 = maxlat;
maparea.lon=[minlon-magextend,maxlon+magextend];
maparea.lat=[minlat-magextend, maxlat+magextend];
% imagearea=maparea;
origgrid = .5;  % longitude grid interval
lonpt = lon1:origgrid:lon2;   % data points lon
latpt = lat1:origgrid:lat2;   % data points lat 
% [gridLON,gridLAT] = meshgrid(lonpt, latpt);

nlonpt=nan(length(lonpt)-1,1); nlatpt=nan(length(latpt)-1,1);
for nn=1:length(lonpt)-1
  nlonpt(nn) = (lonpt(nn)+lonpt(nn+1))/2;   % data points lon
end
for nn=1:length(latpt)-1
  nlatpt(nn) = (latpt(nn)+latpt(nn+1))/2;   % data points lon
end
[ngridLON,ngridLAT] = meshgrid(nlonpt, nlatpt);

%load('us_states.mat');%[]; state(1).polygon(:,1)=Alaska.lon;state(1).polygon(:,2)=Alaska.lat;

fidtemp=fopen(delayfile,'r');
tempdelaydata=textscan(fidtemp,'%*s %*s %f %f %f %f %*f %*f %s');
fclose(fidtemp);
[slat,slon,rlat,rlon,fb] = tempdelaydata{1:5};
clear tempdelaydata;
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ','(j) ','(k) ','(l) '};
% figlabel={'a','b','c','d','e','f','g','h'};

% raymask=zeros(length(nlatpt), length(nlonpt));
figure('Position',[400 450 1300 850]);
for k=1:size(fband,1)
    ftag=['f' num2str(k)];
    % clear maskgcp;
    maskgcp=ones(length(nlatpt), length(nlonpt));
    cc=zeros(length(nlatpt), length(nlonpt));    
    subplot(2,4,k), hold on, box on, axis on    
    set(gca,'ydir','normal');
    daspect([1 cosd(mean([lat1 lat2])) 1]);
    
    clear fm;
    fm=strmatch(char(ftag),fb);
    
    count = length(fm);
    
    ss1lon=slon(fm);
    ss1lat=slat(fm);
    ss2lon=rlon(fm);
    ss2lat=rlat(fm);
    
    for idx=1:count
        deltalon = ss2lon(idx)-ss1lon(idx);
        deltalat = ss2lat(idx)-ss1lat(idx);
    
        if deltalon>=0&&deltalat>=0, index=1;
        elseif deltalon<0&&deltalat>0, index=2;
        elseif deltalon<0&&deltalat<0, index=3;
        elseif deltalon>0&&deltalat<0, index=4;
        end
            
        nx=ceil(abs(deltalon)/origgrid);
        ny=ceil(abs(deltalat)/origgrid);
        slope = abs(deltalat/deltalon);
    %     index
        switch index
        case 1
            if nx>=ny
            dx = ss1lon(idx)+origgrid*(0:nx);
            dy = ss1lat(idx)+origgrid*(0:nx)*slope;
            elseif nx<ny
            dy = ss1lat(idx)+origgrid*(0:ny);
            dx = ss1lon(idx)+origgrid*(0:ny)./slope;  
            end
            
            for pp=1:max(nx,ny)
    %          ii=find(nlonpt>=dx(pp)-origgrid&nlonpt<=dx(pp)+origgrid);
    %          jj=find(nlatpt>=dy(pp)-origgrid&nlatpt<=dy(pp)+origgrid);
             ii = find( abs(nlonpt-dx(pp))==min(abs(nlonpt-dx(pp))) );
             jj = find( abs(nlatpt-dy(pp))==min(abs(nlatpt-dy(pp))) );
             cc(jj(1),ii(1))=cc(jj(1),ii(1))+1;
            end
        
        case 2
            if nx>=ny
            dx = ss1lon(idx)-origgrid*(0:nx);
            dy = ss1lat(idx)+origgrid*(0:nx)*slope;
            elseif nx<ny
            dy = ss1lat(idx)+origgrid*(0:ny);
            dx = ss1lon(idx)-origgrid*(0:ny)./slope;  
            end
            
            for pp=1:max(nx,ny)
    %          ii=find(nlonpt>=dx(pp)-origgrid&nlonpt<=dx(pp)+origgrid);
    %          jj=find(nlatpt>=dy(pp)-origgrid&nlatpt<=dy(pp)+origgrid);
             ii = find( abs(nlonpt-dx(pp))==min(abs(nlonpt-dx(pp))) );
             jj = find( abs(nlatpt-dy(pp))==min(abs(nlatpt-dy(pp))) );
             cc(jj(1),ii(1))=cc(jj(1),ii(1))+1;
            end
              
        case 3      
            if nx>=ny
            dx = ss1lon(idx)-origgrid*(0:nx);
            dy = ss1lat(idx)-origgrid*(0:nx)*slope;
            elseif nx<ny
            dy = ss1lat(idx)-origgrid*(0:ny);
            dx = ss1lon(idx)-origgrid*(0:ny)./slope;  
            end
            
            for pp=1:max(nx,ny)
    %          ii=find(nlonpt>=dx(pp)-origgrid&nlonpt<=dx(pp)+origgrid);
    %          jj=find(nlatpt>=dy(pp)-origgrid&nlatpt<=dy(pp)+origgrid);
             ii = find( abs(nlonpt-dx(pp))==min(abs(nlonpt-dx(pp))) );
             jj = find( abs(nlatpt-dy(pp))==min(abs(nlatpt-dy(pp))) );
             cc(jj(1),ii(1))=cc(jj(1),ii(1))+1;
            end
            
        case 4
            if nx>=ny
            dx = ss1lon(idx)+origgrid*(0:nx);
            dy = ss1lat(idx)-origgrid*(0:nx)*slope;
            elseif nx<ny
            dy = ss1lat(idx)-origgrid*(0:ny);
            dx = ss1lon(idx)+origgrid*(0:ny)./slope;  
            end
            
            for pp=1:max(nx,ny)
    %          ii=find(nlonpt>=dx(pp)-origgrid&nlonpt<=dx(pp)+origgrid);
    %          jj=find(nlatpt>=dy(pp)-origgrid&nlatpt<=dy(pp)+origgrid);
             ii = find( abs(nlonpt-dx(pp))==min(abs(nlonpt-dx(pp))) );
             jj = find( abs(nlatpt-dy(pp))==min(abs(nlatpt-dy(pp))) );
             cc(jj(1),ii(1))=cc(jj(1),ii(1))+1;
            end
        end
    
    end            
    
    % max(max(cc))
    % count
    % if ite==1
    %     cmin=10; cmax=100;cint=20;
    % elseif ite==2
    %     cmin=10; cmax=250;cint=50;
    % elseif ite==3
    %     cmin=10; cmax=400;cint=50;
    % elseif ite==4
    %     cmin=10; cmax=400;cint=50;
    % elseif ite==5
    %     cmin=10; cmax=400;cint=50;
    % elseif ite==6
    %     cmin=10; cmax=400;cint=50;
    % elseif ite==7
    %     cmin=10; cmax=200;cint=50;
    % elseif ite==8
    %     cmin=10; cmax=100;cint=10;
    % end
    cmin=ray_min;%cint=20;
    % maskgcp=cc;
    % maskgcp(cc>=5)=1;
    maskgcp(cc<ray_min)=0;
    % raymask=raymask+maskgcp;
    clear raymask;
    raymask=maskgcp;
    
    % stores all areas with ray coverage.
    clear maskpoints
    raymask(raymask>=1)=1;
    nn=1;
    maskpoints=[];
    for j=1:size(raymask,1)
        for i=1:size(raymask,2)
            if raymask(j,i) == 1
                maskpoints(nn,1)=ngridLON(1,i)-360;
                maskpoints(nn,2)=ngridLAT(j,1);
                nn=nn+1;
            end
        end
    end
    
    %
    clear kk raycover;
    kk=boundary(maskpoints(:,1),maskpoints(:,2),1);
    raycover.data(:,1)=maskpoints(kk,1);
    raycover.data(:,2)=maskpoints(kk,2);
    raycover.cutoff=ray_min;
    
    hi=image(ngridLON(1,:)-360, ngridLAT(:,1), cc);
    hi.CDataMapping='scaled';
    hi.AlphaData=maskgcp;
    set(gca,'CLim',[cmin cmax],'YDir','normal');
    
%     if 1
%         for sb=1:length(state)
%                 plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color','k','LineWidth',.5);
%         end
%     end

    colormap(flip(colormap('jet')))
    hc=colorbar;
    hc.Location='southoutside';
    hc.Label.String=['count >= ' num2str(ray_min)];
    set(hc,'TickDirection','out');
    
    %plot stations.
    % if ite==8
    plot(stationlon,stationlat,'^','color',[.5 .5 .5],'markersize',4,'linewidth',0.5);
    % end
    
    % xlabel('Longitude','FontSize',13)  
    % ylabel('Latitude','FontSize',13)  
%     plot([imagearea.lon(1) imagearea.lon(2) imagearea.lon(2) imagearea.lon(1) imagearea.lon(1)],...
%         [imagearea.lat(1) imagearea.lat(1) imagearea.lat(2) imagearea.lat(2) imagearea.lat(1)],...
%         'k-','linewidth',2); % image area
%     text(mean(maparea.lon), maparea.lat(1)+1.5,[num2str(1/fband(k,2)) '-' num2str(1/fband(k,1)) ' s'],...
%         'fontsize',14,'HorizontalAlignment','center');
    axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
    daspect([1 cosd(mean(maparea.lat)) 1]);
    set(gca,'XTick',maparea.lon(1):6:maparea.lon(2), 'YTick',maparea.lat(1):4:maparea.lat(2));
    set(gca,'TickDir','out','Fontsize',14);
    title([figlabel{k} num2str(1/fband(k,2)) '-' num2str(1/fband(k,1)) ' s'],'FontSize',14);
    
    % plot(raycover.data(:,1),raycover.data(:,2),'m-','linewidth',2); hold off;
    
    hold off;
    drawnow
    
    clear ss1lon ss2lon ss1lat ss2lat
    
    %
    save(['RayCoverOutline_', ite_tag,'_', num2str(1/fband(k,2)),'-',num2str(1/fband(k,1)),...
        's_cutoff',num2str(ray_min),'.mat'],'raycover');
end
% cbarrow('right');
set(gcf,'PaperPositionMode','auto');
eval(['print -dpng -painters ' delayfile '_pathcoverage' num2str(ray_min) '.png']);

end
%% stores all areas with ray coverage.
% clear maskpoints
% % raymask(raymask<10)=0;
% raymask(raymask>=1)=1;
% nn=1;
% maskpoints=[];
% for j=1:size(raymask,1)
%     for i=1:size(raymask,2)
%         if raymask(j,i) > 0.1
%             maskpoints(nn,1)=ngridLON(1,i)-360;
%             maskpoints(nn,2)=ngridLAT(j,1);
%             nn=nn+1;
%         end
%     end
% end
% 
% %
% figure;
% clear kk raycover;
% kk=boundary(maskpoints(:,1),maskpoints(:,2),1);
% raycover.data(:,1)=maskpoints(kk,1);
% raycover.data(:,2)=maskpoints(kk,2);
% raycover.cutoff=5;
% plot(maskpoints(:,1),maskpoints(:,2),'.b'); hold on;
% plot(raycover.data(:,1),raycover.data(:,2),'r-'); hold off;

%%
% save('AlaskaRayCoverOutline_ite05.mat','raycover');
