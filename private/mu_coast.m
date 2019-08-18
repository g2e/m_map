function [ncst,Area,k]=mu_coast(optn,varargin);
% MU_COAST Add a coastline to a given map.
%         MU_COAST draw a coastline as either filled patches (slow) or
%         lines (fast) on a given coastline. It uses a coastline database with
%         a resolution of about 1/4 degree. 
%
%         MU_COAST( (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         MU_COAST('patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%
%         There are a number of options relating to the use of different
%         coastline datasets. If you have coastline data in the M_Map format
%         stored in the .mat file FILENAME (default ./private/m_coasts.mat) 
%         then
%              MU_COAST('user',FILENAME,...) 
%         uses this data in your map. This is not particularly useful unless
%         you make your own coastlines! The primary purpose of this feature is 
%         to enable use of subsampled high-resolution coastline databases
%         for repeated plotting of the same map.
%
%         High-resolution coastlines in the GSHHS format, stored in FILNAME can 
%         be used through the 'gshhs' property:
%
%             MU_COAST(SIZE,FILENAME,...
%
%         However, be warned that the loading this data can be very 
%         timeconsuming. It is perhaps better to subsample the desired data
%         once using
%
%             [ncst,Area,k]=MU_COAST('f','FILENAME);
%             save MMAPFILE ncst Area k
%            
%             MU_COAST('user',MMAPFILE,...)
%
%
%    
%         See also M_PROJ, M_GRID     

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% Notes: 15/June/98 - fixed some obscure problems that occured sometimes
%                     when using conic projections with large extents that
%                     crossed the 180-deg longitude line.
%        31/Aug/98  - added "f" gshhs support (Thanks to RMO)
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;


% m_coasts.mat contains 3 variables:
% ncst: a Nx2 matrix of [LONG LAT] line segments, each of which form
%       a closed contour, separated by NaN
%  k=[find(isnan(ncst(:,1)))];
% Area: a vector giving the areas of the different patches. Both ncst
%     and Area should be ordered with biggest regions first. Area can
%     be computed as follows:
%
%  Area=zeros(length(k)-1,1);
%  for i=1:length(k)-1,
%    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
%    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
%    nl=length(x);
%    Area(i)=sum( diff(x).*(y(1:nl-1)+y(2:nl))/2 );
%  end;
%
%     Area should be >0 for land, and <0 for lakes and inland seas.

switch optn(1),
  case {'c','l','i','h','f'},  
    [ncst,k,Area]=get_coasts(optn,varargin{1});
    varargin(1)=[];
  case  'u',                   
    eval(['load ' varargin{1} ' -mat']);
    varargin(1)=[];
  otherwise
    load m_coasts
end;

% If all we wanted was to extract a sub-coastline, return.
if nargout>0,
  return;
end;

% Handle wrap-arounds (not needed for azimuthal and oblique projections)

switch MAP_PROJECTION.routine,
 case {'mp_cyl','mp_conic','mp_tmerc'}
  if MAP_VAR_LIST.longs(2)<-180,
   ncst(:,1)=ncst(:,1)-360;
  elseif MAP_VAR_LIST.longs(1)>180,
   ncst(:,1)=ncst(:,1)+360;
  elseif MAP_VAR_LIST.longs(1)<-180,
   Area=[Area;Area];
   k=[k;k(2:end)+k(end)-1];
   ncst=[ncst;[ncst(2:end,1)-360 ncst(2:end,2)]];
   % This is kinda kludgey - but sometimes adding all these extra points
   % causes wrap-around in the conic projection, so we want to limit the
   % longitudes to the range needed. However, we don't just clip them to
   % min long because that can cause problems in trying to decide which way
   % curves are oriented when doing the fill algorithm below. So instead
   % I sort of crunch the scale, preserving topology.
   nn=ncst(:,1)<MAP_VAR_LIST.longs(1);
   ncst(nn,1)=(ncst(nn,1)-MAP_VAR_LIST.longs(1))/100+MAP_VAR_LIST.longs(1);
  elseif MAP_VAR_LIST.longs(2)>180,
   Area=[Area;Area];
   k=[k;k(2:end)+k(end)-1];
   ncst=[ncst;[ncst(2:end,1)+360 ncst(2:end,2)]];
   % Ditto.
   nn=ncst(:,1)>MAP_VAR_LIST.longs(2);
   ncst(nn,1)=(ncst(nn,1)-MAP_VAR_LIST.longs(2))/100+MAP_VAR_LIST.longs(2);
  end;
end;

if length(varargin)>0 & strcmp(varargin(1),'patch'),
  optn='patch';
else
  optn='line';
end;


switch optn,
 case 'patch',

  switch MAP_VAR_LIST.rectbox,
    case 'on',
      xl=MAP_VAR_LIST.xlims;
      yl=MAP_VAR_LIST.ylims;
      [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
      oncearound=4;
    case 'off',
      xl=MAP_VAR_LIST.longs;
      yl=MAP_VAR_LIST.lats;
      X=ncst(:,1);
      Y=ncst(:,2);
      [X,Y]=mu_util('clip','on',X,xl(1),X<xl(1),Y);
      [X,Y]=mu_util('clip','on',X,xl(2),X>xl(2),Y);
      [Y,X]=mu_util('clip','on',Y,yl(1),Y<yl(1),X);
      [Y,X]=mu_util('clip','on',Y,yl(2),Y>yl(2),X);
      oncearound=4;
    case 'circle',
      rl=MAP_VAR_LIST.rhomax;
      [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
      oncearound=2*pi;    
  end;

  for i=1:length(k)-1,
    x=X(k(i)+1:k(i+1)-1);
    fk=finite(x);
    if any(fk),
      y=Y(k(i)+1:k(i+1)-1);
      nx=length(x);
      if Area(i)<0, x=flipud(x);y=flipud(y); fk=flipud(fk); end;
%clf
%line(x,y,'color','m');
      st=find(diff(fk)==1)+1;
      ed=find(diff(fk)==-1);
%length(x),
%ed,
%st
       if length(st)<length(ed) | isempty(st), st=[ 1;st]; end;
       if length(ed)<length(st),               ed=[ed;nx]; end;
%ed
%st
      if  ed(1)<st(1),
        if c_edge(x(1),y(1))==9999,
          x=x([(ed(1)+1:end) (1:ed(1))]);
          y=y([(ed(1)+1:end) (1:ed(1))]);
          fk=finite(x);
          st=find(diff(fk)==1)+1;
          ed=[find(diff(fk)==-1);nx];
          if length(st)<length(ed), st=[1;st]; end
        else
          ed=[ed;nx];
          st=[1;st];
        end;
      end;
%ed
%st
      % Get rid of 2-point curves (since often these are out-of-range lines with
      % both endpoints outside the region of interest)
      k2=(ed-st)<3;
      ed(k2)=[]; st(k2)=[];
%%[XX,YY]=m_ll2xy(x(st),y(st),'clip','off');
%line(x,y,'color','r','linest','-');
%line(x(st),y(st),'marker','o','color','r','linest','none');
%line(x(ed),y(ed),'marker','o','color','g','linest','none');
      edge1=c_edge(x(st),y(st));
      edge2=c_edge(x(ed),y(ed));
%-edge1*180/pi
%-edge2*180/pi
      mi=1;
      while length(st)>0,
        if mi==1,
          xx=x(st(1):ed(1));
          yy=y(st(1):ed(1));
        end;
        estart=edge2(1);
        s_edge=edge1;
        s_edge(s_edge<estart)=s_edge(s_edge<estart)+oncearound;
%s_edge,estart
        [md,mi]=min(s_edge-estart);
        switch MAP_VAR_LIST.rectbox,
          case {'on','off'},
            for e=floor(estart):floor(s_edge(mi)),
              if e==floor(s_edge(mi)), xe=x(st([mi mi])); ye=y(st([mi mi])); 
              else  xe=xl; ye=yl; end;
              switch rem(e,4),
                case 0,
                  xx=[xx; xx(end*ones(10,1))];
                  yy=[yy; yy(end)+(ye(2)-yy(end))*[1:10]'/10 ];
                case 1,
                  xx=[xx; xx(end)+(xe(2)-xx(end))*[1:10]'/10 ];
                  yy=[yy; yy(end*ones(10,1))];
                case 2,
                  xx=[xx; xx(end*ones(10,1))];
                  yy=[yy; yy(end)+(ye(1)-yy(end))*[1:10]'/10;];
                case 3,
                  xx=[xx; xx(end)+(xe(1)-xx(end))*[1:10]'/10 ];
                  yy=[yy; yy(end*ones(10,1))];
              end;
            end;
          case 'circle',
            if estart~=9999,
%s_edge(mi),estart
              xx=[xx; rl*cos(-[(estart:.1:s_edge(mi))]')];
              yy=[yy; rl*sin(-[(estart:.1:s_edge(mi))]')];
            end;
        end;
        if mi==1, % joined back to start
          if strcmp(MAP_VAR_LIST.rectbox,'off'), 
            [xx,yy]=m_ll2xy(xx,yy,'clip','off'); 
          end;
%disp(['paused-1 ' int2str(i)]);pause;
          if Area(i)<0,
            patch(xx,yy,varargin{2:end},'facecolor',get(gca,'color'));
          else
            patch(xx,yy,varargin{2:end});
          end;
          ed(1)=[];st(1)=[];edge2(1)=[];edge1(1)=[];
        else
          xx=[xx;x(st(mi):ed(mi))];
          yy=[yy;y(st(mi):ed(mi))];
          ed(1)=ed(mi);st(mi)=[];ed(mi)=[];
          edge2(1)=edge2(mi);edge2(mi)=[];edge1(mi)=[];
        end;
%%disp(['paused-2 ' int2str(i)]);pause;
      end;
    end; 
  end;  

 otherwise,
  [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','on');
 
  % Get rid of 2-point lines (these are probably clipped lines spanning the window)
  fk=finite(X);        
  st=find(diff(fk)==1)+1;
  ed=find(diff(fk)==-1);
  k=find((ed-st)==1);
  X(st(k))=NaN;

  line(X,Y,varargin{:}); 

end;

%-----------------------------------------------------------------------
function edg=c_edge(x,y);
% C_EDGE tests if a point is on the edge or not. If it is, it is
%        assigned a value representing it's position oon the perimeter
%        in the clockwise direction. For x/y or lat/long boxes, these
%        values are
%           0 -> 1 on left edge
%           1 -> 2 on top
%           2 -> 3 on right edge
%           3 -> 4 on bottom
%        For circular boxes, these values are the -ve of the angle
%        from center.

global MAP_VAR_LIST

switch MAP_VAR_LIST.rectbox,
  case 'on',
    xl=MAP_VAR_LIST.xlims;
    yl=MAP_VAR_LIST.ylims;
  case 'off',
    xl=MAP_VAR_LIST.longs;
    yl=MAP_VAR_LIST.lats;
  case 'circle',
    rl2=MAP_VAR_LIST.rhomax^2;
end;

edg=9999+zeros(length(x),1);
tol=1e-10;

switch MAP_VAR_LIST.rectbox,
  case {'on','off'},
    i=abs(x-xl(1))<tol;
    edg(i)=(y(i)-yl(1))/diff(yl);

    i=abs(x-xl(2))<tol;
    edg(i)=2+(yl(2)-y(i))/diff(yl);

    i=abs(y-yl(1))<tol;
    edg(i)=3+(xl(2)-x(i))/diff(xl);

    i=abs(y-yl(2))<tol;
    edg(i)=1+(x(i)-xl(1))/diff(xl);

  case 'circle',
    i=abs(x.^2 + y.^2 - rl2)<tol;
    edg(i)=-atan2(y(i),x(i));   % use -1*angle so that numeric values
                                % increase in CW direction
end;



%%
function [ncst,k,Area]=get_coasts(optn,file);
%
%  GET_COASTS  Loads various GSHHS coastline databases and does some preliminary
%              processing to get things into the form desired by the patch-filling
%              algorithm.


global MAP_PROJECTION MAP_VAR_LIST

llim=rem(MAP_VAR_LIST.longs(1)+360,360)*1e6;
rlim=rem(MAP_VAR_LIST.longs(2)+360,360)*1e6;
tlim=MAP_VAR_LIST.lats(2)*1e6;
blim=MAP_VAR_LIST.lats(1)*1e6;

mrlim=rem(MAP_VAR_LIST.longs(2)+360+180,360)-180;
mllim=rem(MAP_VAR_LIST.longs(1)+360+180,360)-180;
mtlim=MAP_VAR_LIST.lats(2);
mblim=MAP_VAR_LIST.lats(1);

switch optn(1),
  case 'f',   % 'full' (undecimated) database
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
  case 'h',
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
  case 'i',
    ncst=NaN+zeros(492283,2);Area=zeros(41520,1);k=ones(41521,1);
  case 'l',
    ncst=NaN+zeros(101023,2);Area=zeros(10768,1);k=ones(10769,1);
  case 'c',
    ncst=NaN+zeros(14872,2);Area=zeros(1868,1);k=ones(1869,1);
end;
fid=fopen(file,'r','ieee-be');

if fid==-1,
  warning(sprintf(['Coastline file ' file ...
          ' not found \n(Have you installed it? See the User''s Guide for details)' ...
          '\n ---Using default coastline instead']));
  load m_coasts
  return
end;


Area2=Area;

[A,cnt]=fread(fid,9,'int32');

l=0;
while cnt>0,

 C=fread(fid,A(2)*2,'int32');
% pause;
 
 a=rlim>llim;
 b=A(9)<65536;  %%rem(A(4)+360e6,360e6)<rem(A(5)+360e6,360e6); Cross boundary?
 c=llim<rem(A(5)+360e6,360e6);
 d=rlim>rem(A(4)+360e6,360e6);
 
 if a&(b&c&d | ~b&(c|d)) | ~a&(~b | (b&(c|d))),
 
   l=l+1;
 
   x=C(1:2:end)*1e-6;y=C(2:2:end)*1e-6;

   %  make things continuous (join edges that cut across 0-meridian)

   dx=diff(x);
   if A(9)>65536 | any(dx)>356 | any(dx<356),
     x=x-360*cumsum([x(1)>180;(dx>356) - (dx<-356)]);
   end;

   % Antarctic is a special case - extend contour to make nice closed polygon
   % that doesn't surround the pole.   
   if abs(x(1))<1 & abs(y(1)+68.9)<1,
     y=[-89.9;-78.4;y(x<=-180);y(x>-180);   -78.4;-89.9*ones(18,1)];
     x=[  180; 180 ;x(x<=-180)+360;x(x>-180);-180; [-180:20:160]'];
   end;

   % First and last point should be the same.
   
   if x(end)~=x(1) | y(end)~=y(1), x=[x;x(1)];y=[y;y(1)]; end;

   % get correct curve orientation for patch-fill algorithm.
   
   Area2(l)=sum( diff(x).*(y(1:(end-1))+y(2:end))/2 );
   Area(l)=A(8)/10;

   if rem(A(3),2)==0; 
     Area(l)=-abs(Area(l)); 
     if Area2(l)>0, x=x(end:-1:1);y=y(end:-1:1); end;
   else
     if Area2(l)<0, x=x(end:-1:1);y=y(end:-1:1); end; 
   end;

   % Here we try to reduce the number of points.
   
   xflag=0;
   if max(x)>180, % First, save original curve for later if we anticipate
     sx=x;sy=y;   % a 180-problem.
     xflag=1;
   end;
   
   % Look for points outside the lat/long boundaries, and then decimate them
   % by a factor of about 20 (don't get rid of them completely because that
   % can sometimes cause problems when polygon edges cross curved map edges).
   
   tol=.2;   
   nn=y>mtlim+tol;  nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
   nn=nn & rem(1:length(nn),20)'~=0;
   x(nn)=[];y(nn)=[];
   nn=y<mblim-tol;  nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
   nn=nn & rem(1:length(nn),20)'~=0;
   x(nn)=[];y(nn)=[];
   if mrlim>mllim,
     nn=x>mrlim+tol;nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
   nn=nn & rem(1:length(nn),20)'~=0;
     x(nn)=[];y(nn)=[];
     nn=x<mllim-tol;nn=logical(nn-([0;diff(nn)]>0)-([diff(nn);0]<0));nn([1 end])=0;
   nn=nn & rem(1:length(nn),20)'~=0;
     x(nn)=[];y(nn)=[];
   end;
 %% plot(x,y);pause;
   k(l+1)=k(l)+length(x)+1;
   ncst(k(l)+1:k(l+1)-1,:)=[x,y];
   ncst(k(l+1),:)=[NaN NaN];
   
   % This is a little tricky...the filling algorithm expects data to be in the
   % range -180 to 180 deg long. However, there are some land parts that cut across
   % this divide so they appear at +190 but not -170. This causes problems later...
   % so as a kludge I replicate some of the problematic features at 190-360=-170.
   % Small islands are just duplicated, for the Eurasian landmass I just clip
   % of the eastern part.
   
   if xflag,
     l=l+1;Area(l)=Area(l-1); 
     if abs(Area(l))>1e5,
       nn=find(sx>180);nn=[nn;nn(1)];
       k(l+1)=k(l)+length(nn)+1;
       ncst(k(l)+1:k(l+1)-1,:)=[sx(nn)-360,sy(nn)];
     else   % repeat the island at the other edge.
       k(l+1)=k(l)+length(sx)+1;
       ncst(k(l)+1:k(l+1)-1,:)=[sx-360,sy];
     end;
     ncst(k(l+1),:)=[NaN NaN];
   end;
 end;
 
 
 [A,cnt]=fread(fid,9,'int32');

end;

fclose(fid);

%%plot(ncst(:,1),ncst(:,2));pause;clf;
%size(ncst)
%size(Area)
%size(k)


ncst((k(l+1)+1):end,:)=[];
Area((l+1):end)=[];
k((l+2):end)=[];


%size(ncst)
%size(Area)
%size(k)






























