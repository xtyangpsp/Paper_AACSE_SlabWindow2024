% Script combining CUB2.0 (0-396 km depth) and AK135 at great depth to generate a 
% media volume for finite-difference simulation in the spherical coordinate system
% The lat and long dimensions are defined in read_CUB.m
% Version: modified from the script in ./SE.Tibet/config/structure, Y.S. 03/2010
%
% Xiaotao Yang: extend the model to the depth of 1057 from ak135.
% added reference for the conversion from Vp to density for the crust.
% reference for mantle conversion is still not clear.

clear all

%------------------- read AK135     -----------
ak135=load('./vmodels/ak135.tvel');
[mak, nak]=size(ak135);

%
define_latlon;
%------------------ read CUB ------------------
% study area is defined in read_CUB.m, which calls define_latlon.m
[lat, lon, dep, vs]=read_CUB_revised(minlon,maxlon,minlat,maxlat);
lat=double(lat);
lon=double(lon);
dep=double(dep);
vs=double(vs);
% In CUB the surface sediment has a vs of 1.197 km/s and vs<1.197 means water.
% Since the model is discretized at 4 km in the radial direction, a composite model
% using  a bathymetry dataset would be more accurate. 
%vs(find(vs<1.197))=1.197; % replace the shallow water in the Bay of Bengal with soft sed. 
%vs(:,:,1)=vs(:,:,2);	% On the plateau, the soft sediment are thin and 
			% localized. So we replace this ubiquitous low velocity layer in CUB
			% with the velocity at 4 km depth.  This modification may make the
			% surface velocity in the basins too large.
figure;
testidx=2;
disp(dep(testidx))
imagesc(lon,lat,vs(:,:,testidx));set(gca,'ydir','normal');
size(vs)
%%
colat=90-lat;
nlat=numel(lat);
nlon=numel(lon);
ndep=numel(dep);

moho=zeros(nlat,nlon);
%------------------- convert Vs to Vp  -----------
for i=1:nlat
  for j=1:nlon
    for k=1:ndep
      depth=(k-1)*4;
      if(vs(i,j,k)<4.0);  % proxy for crust
        vp(i,j,k)=vs(i,j,k)*1.74;
        moho(i,j)=k;
      else
	% use vp/vs of ak135
        for m=1:mak
          if(ak135(m,1)>depth), break; end;
        end
        ratu=ak135(m-1,2)/ak135(m-1,3);
        ratd=ak135(m,2)/ak135(m,3);
        vp(i,j,k)=vs(i,j,k)*(ratu+(ratd-ratu)/(ak135(m,1)-ak135(m-1,1))*(depth-ak135(m-1,1)));
      end
      % water layer
      if(vs(i,j,k)<=0.0001)
% 	    vp(i,j,k)=1.5;
        if dep(k)>4
            vs(i,j,k)=2.5;
        else
            vs(i,j,k)=1.;
        end
        
        vp(i,j,k)=vs(i,j,k)*1.74;
        moho(i,j)=k;
      end
    end
  end
end

%-------------------calculate density using vp-density relation  -----------
for i=1:nlat
  for j=1:nlon
    % for the mantle (reference?)
    for k=ndep:-1:moho(i,j)+1
      if(k==ndep)
        rho(i,j,k)=3.5068+(vp(i,j,k)-9.0302)*0.384;
      else
        rho(i,j,k)=rho(i,j,k+1)+(vp(i,j,k)-vp(i,j,k+1))*0.384;
      end
    end
    % for the crust (reference: Christensen and Mooney, JGR1995, Table 7)
    for k=moho(i,j):-1:1
      if(dep(k)<=10),rho(i,j,k)=0.9893+0.2891*vp(i,j,k);
      elseif(dep(k)>10 && dep(k)<=20),rho(i,j,k)=0.9473+0.2966*vp(i,j,k);
      elseif(dep(k)>20 && dep(k)<=30),rho(i,j,k)=0.9466+0.2997*vp(i,j,k);
      elseif(dep(k)>30 && dep(k)<=40),rho(i,j,k)=0.9645+0.3005*vp(i,j,k);
      elseif(dep(k)>40),rho(i,j,k)=1.0783+0.2990*vp(i,j,k); 
      end;
      % water layer
      if(vs(i,j,k)==0); rho(i,j,k)=1.0; end
    end
  end
end

ndep=ndep+1;
dep(ndep)=409;
vp(:,:,ndep)=9.0302;
vs(:,:,ndep)=4.8702;
rho(:,:,ndep)=3.5068;

ndep=ndep+1;
dep(ndep)=410;
vp(:,:,ndep)=9.3601;
vs(:,:,ndep)=5.0806;
rho(:,:,ndep)=3.9317;

ndep=ndep+1;
dep(ndep)=460;
vp(:,:,ndep)=9.5208;
vs(:,:,ndep)=5.1864;
rho(:,:,ndep)=3.9273;

ndep=ndep+1;
dep(ndep)=510;
vp(:,:,ndep)=9.6962;
vs(:,:,ndep)=5.2922;
rho(:,:,ndep)=3.9233;

ndep=ndep+1;
dep(ndep)=560;
vp(:,:,ndep)=9.8640;
vs(:,:,ndep)=5.3989;
rho(:,:,ndep)=3.9218;

ndep=ndep+1;
dep(ndep)=610;
vp(:,:,ndep)=10.0320;
vs(:,:,ndep)=5.5047;
rho(:,:,ndep)=3.9206;

ndep=ndep+1;
dep(ndep)=659;
vp(:,:,ndep)=10.2000;
vs(:,:,ndep)=5.6104;
rho(:,:,ndep)=3.9201;

ndep=ndep+1;
dep(ndep)=660;
vp(:,:,ndep)=10.7909;
vs(:,:,ndep)=5.9607;
rho(:,:,ndep)=4.2387;

ndep=ndep+1;
dep(ndep)=710;
vp(:,:,ndep)=10.9222;
vs(:,:,ndep)=6.0898;
rho(:,:,ndep)=4.2986;

ndep=ndep+1;
dep(ndep)=760;
vp(:,:,ndep)=11.0558;
vs(:,:,ndep)=6.2095;
rho(:,:,ndep)=4.4305;

ndep=ndep+1;
dep(ndep)=809.5;
vp(:,:,ndep)=11.1353;
vs(:,:,ndep)=6.2426;
rho(:,:,ndep)=4.4596;

ndep=ndep+1;
dep(ndep)=859;
vp(:,:,ndep)=11.2221;
vs(:,:,ndep)=6.2798;
rho(:,:,ndep)=4.4885;

ndep=ndep+1;
dep(ndep)=908.5;
vp(:,:,ndep)=11.3068;
vs(:,:,ndep)=6.3160;
rho(:,:,ndep)=4.5173;

ndep=ndep+1;
dep(ndep)=958;
vp(:,:,ndep)=11.3896;
vs(:,:,ndep)=6.3512;
rho(:,:,ndep)=4.5459;

ndep=ndep+1;
dep(ndep)=1007.5;
vp(:,:,ndep)=11.4705;
vs(:,:,ndep)=6.3854;
rho(:,:,ndep)=4.5744;

ndep=ndep+1;
dep(ndep)=1057;
vp(:,:,ndep)=11.5495;
vs(:,:,ndep)=6.4187;
rho(:,:,ndep)=4.6028;

% convert to SI units
vp=vp*1000;
vs=vs*1000;
rho=rho*1000;
dep=dep*1000;


lon_test=min(lon):0.2:max(lon);
lat_test=min(lat):0.1:max(lat);
[xmesh,ymesh,zmesh]=meshgrid(lon_test,lat_test,dep);
vs=interp3(lon,lat,dep,vs,xmesh,ymesh,zmesh,'linear'); 
vp=interp3(lon,lat,dep,vp,xmesh,ymesh,zmesh,'linear');  
rho=interp3(lon,lat,dep,rho,xmesh,ymesh,zmesh,'linear');  
vs=smooth3(vs,'box',[5 5 3]);
vp=smooth3(vp,'box',[5 5 3]);
rho=smooth3(rho,'box',[5 5 3]);

lon=lon_test;
lat=lat_test;
nlat=numel(lat);
nlon=numel(lon);
ndep=numel(dep);
colat=90-lat;

size(vp)
vp=permute(vp,[3 2 1]); % depth lon lat
size(vp)
vp=flipdim(vp,3);
size(vp)
vs=permute(vs,[3 2 1]); % depth lon lat
vs=flipdim(vs,3);
rho=permute(rho,[3 2 1]); % depth lon lat
rho=flipdim(rho,3);
colat=fliplr(colat);
%------------------- create output nc file -----------------------
fnm_out='SeisMedia.volume.AACSE_area.CUB2_AK135.nc';
if 1
my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
nc_create_empty ( fnm_out, my_mode );
% nc_create_empty(fnm_out);
nc_add_dimension(fnm_out,'phi',nlon);
nc_add_dimension(fnm_out,'theta',nlat);
nc_add_dimension(fnm_out,'depth',ndep);

nc_attput(fnm_out,nc_global,'sealevel',6371*1e3);

var.Nctype='float';var.Attribute=[];
var.Name='theta';var.Dimension={'theta'};nc_addvar(fnm_out,var);
var.Name='phi';var.Dimension={'phi'};nc_addvar(fnm_out,var);
var.Name='depth';var.Dimension={'depth'};nc_addvar(fnm_out,var);
var.Name='depth2sealevel';var.Dimension={'depth'};nc_addvar(fnm_out,var);

var.Name='Vp'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
var.Name='Vs'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
var.Name='rho';var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);

nc_varput(fnm_out,'theta',colat);
nc_varput(fnm_out,'phi', lon);
nc_varput(fnm_out,'depth2sealevel',dep);

nc_varput(fnm_out,'Vp',vp);
nc_varput(fnm_out,'Vs',vs);
nc_varput(fnm_out,'rho',rho);

disp('finished creating')
end

%%
testidx=2
vs2=ncread(fnm_out,'Vs');
lat2=90-ncread(fnm_out,'theta');
lon2=ncread(fnm_out,'phi');
size(vs2)
figure
imagesc(lon2-360,lat2,squeeze(vs2(:,:,testidx)/1000));set(gca,'ydir','normal');

