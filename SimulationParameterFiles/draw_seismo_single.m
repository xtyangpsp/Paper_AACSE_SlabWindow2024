% draw_seismo: Plot seismography from inline or recv output.

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Initial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set_mfiles_path
% [fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
%   dir_station,dir_out]=get_simul_path
fnm_conf='./SeisFD3D.conf';
dir_coord='./input';  
dir_metric='./input'; 
dir_media='./input'; 
dir_source='./input'; 
dir_station='./input';
dir_out='./output';
% ----------------------- parameter -----------------------
spec_line='r';

%N=4; Fc=0.1;
N=2; Fc=[0.05 0.1] ; %Fc= 0.2;

cmp=[];cmp0=[];cmp1=[];
cmp{end+1}= 'Ux'; cmp0{end+1}='Ux0'; cmp1{end+1}='( 1)';
cmp{end+1}= 'Uy'; cmp0{end+1}='Uy0'; cmp1{end+1}='( 1)';
cmp{end+1}= 'Uz'; cmp0{end+1}='Uz0'; cmp1{end+1}='( 1)';
% cmp{end+1}= 'Ux'; cmp0{end+1}='1'; cmp1{end+1}='( 5)';
% cmp{end+1}= 'Uy'; cmp0{end+1}='1'; cmp1{end+1}='( 5)';
% cmp{end+1}= 'Uz'; cmp0{end+1}='1'; cmp1{end+1}='( 5)';
%cmp{end+1}= 'Vx'; cmp0{end+1}='1'; cmp1{end+1}='( 1)';
%cmp{end+1}= 'Vy'; cmp0{end+1}='1'; cmp1{end+1}='( 1)';
%cmp{end+1}= 'Vz'; cmp0{end+1}='1'; cmp1{end+1}='( 1)';
ncmp=numel(cmp);

% -------------------- load data --------------------------
ryid=[];stid=[];
%ryid(end+1)= 0; stid(end+1)= 1;
%ryid(end+1)= 0; stid(end+1)= 2;
%ryid(end+1)= 0; stid(end+1)= 3;
%ryid(end+1)= 0; stid(end+1)= 4;
%ryid(end+1)= 0; stid(end+1)= 5;
%ryid(end+1)= 0; stid(end+1)= 6;
%ryid(end+1)= 0; stid(end+1)= 7;
%ryid(end+1)= 0; stid(end+1)= 8;

%ryid(end+1)= 1; stid(end+1)= 25*(5+10)+1;spec_line='b';
ryid(end+1)=2;
stid(end+1)=1;

spec_line='b';

%for n=5:5:750
%    ryid(end+1)=2; stid(end+1)=n;
%end
%ryid=fliplr(ryid); stid=fliplr(stid);

nrecv=numel(ryid);

for n=1:nrecv
    seismoinfo=locate_seismo(fnm_conf,ryid(n),stid(n),dir_station,@get_fnm_station);
    coord_grid(n,:)=retrieve_station(seismoinfo,'grid','stationdir',dir_station);
     Sx    = retrieve_seismo(seismoinfo,'Vx','outdir',dir_out);
     Sy    = retrieve_seismo(seismoinfo,'Vy','outdir',dir_out);
    [Sz,T] = retrieve_seismo(seismoinfo,'Vz','outdir',dir_out);
    Vx(n,:)=Sx; Vy(n,:)=Sy; Vz(n,:)=Sz; Vt(n,:)=T;
    stept=T(2)-T(1);
    Ax(n,:)=gradient(Sx,stept); Ay(n,:)=gradient(Sy,stept); Az(n,:)=gradient(Sz,stept);
    Ux(n,:)=cumtrapz(Sx)*stept; Uy(n,:)=cumtrapz(Sy)*stept; Uz(n,:)=cumtrapz(Sz)*stept;
    At=Vt; Ut=Vt;
end

% -- filter --
if exist('Fc','var')
for n=1:nrecv
   stept=Vt(n,2)-Vt(n,1); Feff=1/stept/2; Wn=Fc/Feff; [b,a]=butter(N,Wn);
   Vx(n,:)=filter(b,a,Vx(n,:)); Vy(n,:)=filter(b,a,Vy(n,:)); Vz(n,:)=filter(b,a,Vz(n,:));
   stept=At(n,2)-At(n,1); Feff=1/stept/2; Wn=Fc/Feff; [b,a]=butter(N,Wn);
   Ax(n,:)=filter(b,a,Ax(n,:)); Ay(n,:)=filter(b,a,Ay(n,:)); Az(n,:)=filter(b,a,Az(n,:));
   stept=Ut(n,2)-Ut(n,1); Feff=1/stept/2; Wn=Fc/Feff; [b,a]=butter(N,Wn);
   Ux(n,:)=filter(b,a,Ux(n,:)); Uy(n,:)=filter(b,a,Uy(n,:)); Uz(n,:)=filter(b,a,Uz(n,:));
end
end

Vx0=max(max(abs(Vx)));Vy0=max(max(abs(Vy)));Vz0=max(max(abs(Vz)));
Ax0=max(max(abs(Ax)));Ay0=max(max(abs(Ay)));Az0=max(max(abs(Az)));
Ux0=max(max(abs(Ux)));Uy0=max(max(abs(Uy)));Uz0=max(max(abs(Uz)));
V0=max([Vx0,Vy0,Vz0]);A0=max([Ax0,Ay0,Az0]);U0=max([Ux0,Uy0,Uz0]);

% -------------------- plot figures ------------------------
figure('Position',[500 400 800 700])
for m=1:ncmp
    %figure;
    subplot(3,1,m)
    for n=1:nrecv
        T=eval([cmp{m}(1) 't(' num2str(n) ',:)']);
        W=eval([cmp{m} '(' num2str(n) ',:)'])/eval(cmp0{m})*eval(cmp1{m});
        plot(T,W,spec_line);
        hold on;
    end
    title(cmp{m});
end

saveas(gca,'seismogram_a1.5t4.5.jpg','jpg');