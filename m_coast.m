function m_coast(varargin);
% M_COAST Add a coastline to a given map.
%         M_COAST draw a coastline as either filled patches (slow) or
%         lines (fast) on a given coastline. It uses a coastline database with
%         a resolution of about 1/4 degree. 
%
%         M_COAST( (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_COAST('patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%    
%         See also M_PROJ, M_GRID     

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
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

load m_coasts

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
  elseif MAP_VAR_LIST.longs(2)>180,
   Area=[Area;Area];
   k=[k;k(2:end)+k(end)-1];
   ncst=[ncst;[ncst(2:end,1)+360 ncst(2:end,2)]];
  end;
end;

if nargin>0 & strcmp(varargin(1),'patch'),
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
      if Area(i)<0, x=flipud(x);y=flipud(y); end;
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


