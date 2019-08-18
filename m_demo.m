function m_demo(num);
% M_DEMO  Demonstration program showing various maps in M_Map package
%         Dig into this to look for examples of things you want to do.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 7/May/1997
% (thanks to Art Newhall for putting these examples into an m-file).
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

clf;
m_proj('ortho','lat',48','long',-123');
m_coast('patch','r');
m_grid('linest','-','xticklabels',[],'yticklabels',[],'ytick',[-80:40:80]);
xlabel('Orthographic Projection','visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
m_proj('lambert','long',[-160 -40],'lat',[30 80]);
m_coast('patch',[1 .85 .7]);
m_elev('contourf',[500:500:6000]);
m_grid;
colormap(flipud(copper));
xlabel('Conic Projection of North America with elevations','visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clf;
m_proj('stereographic','lat',90,'long',30,'radius',25);
m_elev('contour',[-3500:1000:-500],'edgecolor','b');
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','r');
xlabel('Polar Stereographic Projection with bathymetry','visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clf;
Slongs=[-100 0;-75 25;-5 45; 25 145;45 100;145 295;100 290];
Slats= [  8 80;-80  8; 8 80;-80   8; 8  80;-80   0;  0  80];
for l=1:7,
 m_proj('sinusoidal','long',Slongs(l,:),'lat',Slats(l,:));
 m_coast('patch','g');
 m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],...
        'ytick',[-80:20:80],'yticklabels',[],'linest','-','color',[.9 .9 .9]);
end;
xlabel('Interrupted Sinusoidal Projection of World Oceans','visible','on');

% The multiple maps trick is useful only with this projection. In order to
% see all the maps we must undo the axis limits set by m_grid calls:

set(gca,'xlimmode','auto','ylimmode','auto');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clf
%% Nice looking data
[lon,lat]=meshgrid([-136:2:-114],[36:2:54]);
u=sin(lat/6);
v=sin(lon/6);

m_proj('oblique','lat',[56 30],'lon',[-132 -120],'aspect',.8);

subplot(121);
m_coast('patch',[.9 .9 .9],'edgecolor','none');
m_grid('tickdir','out','yaxislocation','right',...
       'xaxislocation','top','xlabeldir','end','ticklen',.02);
hold on;
m_quiver(lon,lat,u,v);
xlabel('Simulated surface winds');

subplot(122);
m_coast('patch',[.9 .9 .9],'edgecolor','none');
m_grid('tickdir','out','yticklabels',[],...
       'xticklabels',[],'linestyle','none','ticklen',.02);
hold on;
[cs,h]=m_contour(lon,lat,sqrt(u.*u+v.*v));
clabel(cs,h,'fontsize',8);
xlabel('Simulated something else');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clf
% Plot a circular orbit
lon=[-180:180];
lat=atan(tan(60*pi/180)*cos((lon-30)*pi/180))*180/pi;

m_proj('miller','lat',82);
m_coast('color',[0 .6 0]);
m_grid('linest','none');

m_line(lon,lat,'linewi',3,'color','r');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  hit return to continue');
pause
disp('        ...drawing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clf
m_proj('lambert','lon',[-10 20],'lat',[33 48]);
m_tbase('contourf');
m_grid('linestyle','none');
colormap(jet);
