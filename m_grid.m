function m_grid(varargin);
% M_GRID make a grid on a map.
%        M_GRID('parameter','value',...) with any number (or no)
%        optional parameters is used to draw a lat/long grid for a
%        previously initialized map projection.
%
%        The optional parameters allow the user
%        to control the look of the grid. These parameters are listed
%        by MGRID_('get'), with defualt parameters in M_GRID('set');
%
%        see also M_PROJ

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%  19/6/97 - set visibility of titles and so forth to 'on' (they
%            default to 'off' when axes visibility is turned off)
%  2/11/97 - for M5.1, the old way of making the patch at the bottom (i.e.
%            by rearranging the axes children) instead causes matlab to loose
%            track of titles. Try a different fix.
% 11/01/98 - Added way of making longitude lines cut off to prevent crowding near poles (you have
%            to specify a vector for allowabale latitudes for this to work).
% 16/02/98 - Made a little fudge to allow the user to fully specify grid location
%            without getting the edge points. It doesn't quite work if only *one* edge
%            point is desired....but I hope it will be OK that way.
% 19/02/98 - PC-users complain about layers getting out of order! Their fault for using
%            such an awful OS...however (with help from Eric Firing) I think I have
%            a fix.
%  7/04/98 - Another fix to grid locations to not automatically add edge points
%            (as requested by EF)

% Note that much of the work in generating line data 
% is done by calls to the individual projections - 
% most of M_GRID is concerned with the mechanics of plotting


% These structures are initialized by m_proj()

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;


% Otherwise we are drawing a grid!

% Default parameters for grid

xtick=6;
ytick=6;
xlabels=NaN;
ylabels=NaN;
gcolor='k';
glinestyle=':';
glinewidth=get(gca,'linewidth');
gbox='on';
gfontsize=get(gca,'fontsize');
gfontname=get(gca,'fontname');
gxaxisloc=get(gca,'xaxislocation'); 
gyaxisloc=get(gca,'yaxislocation');
gtickdir=get(gca,'tickdir'); 
gticklen=get(gca,'ticklen'); gticklen=gticklen(1); 
gxticklabeldir='middle';
gyticklabeldir='end';

% Parse parameter list for options. I really should do some
% error checking here, but...

k=1;
while k<=length(varargin),
  switch lower(varargin{k}(1:3)),
    case 'box',
      gbox=varargin{k+1};
    case 'xti',
      if length(varargin{k})==5,
        xtick=varargin{k+1};
      else
        xlabels=varargin{k+1};
      end;
    case 'yti',
      if length(varargin{k})==5,
        ytick=varargin{k+1};
      else
        ylabels=varargin{k+1};
      end;
    case 'xla',
      gxticklabeldir=varargin{k+1};
    case 'yla',
      gyticklabeldir=varargin{k+1};
    case 'col',
      gcolor=varargin{k+1};
    case 'lin',
      switch lower(varargin{k}(1:5)),
         case 'linew',
           glinewidth=varargin{k+1};
         case 'lines',
           glinestyle=varargin{k+1};
      end;
    case 'fon',
       switch lower(varargin{k}(1:5)),
         case 'fonts',
           gfontsize=varargin{k+1};
         case 'fontn',
           gfontname=varargin{k+1};
         end;
    case 'xax',
      gxaxisloc=varargin{k+1};
    case 'yax',
      gyaxisloc=varargin{k+1};
    case 'tic',
      switch lower(varargin{k}(1:5)),
        case 'tickl',
           gticklen=varargin{k+1};
        case 'tickd',
           gtickdir=varargin{k+1};
        end;
    case {'get','usa'},
      disp('      ''box'',( ''on'' | ''off'' )');
      disp('      ''xtick'',( num | [value1 value2 ...])');
      disp('      ''ytick'',( num | [value1 value2 ...])');
      disp('      ''xticklabels'',[label1;label2 ...]');
      disp('      ''yticklabels'',[label1;label2 ...]');
      disp('      ''xlabeldir'', ( ''middle'' | ''end'' )');
      disp('      ''ylabeldir'', ( ''end'' | ''middle'' )');
      disp('      ''ticklen'',value');
      disp('      ''tickdir'',( ''in'' | ''out'' )');
      disp('      ''color'',colorspec');
      disp('      ''linewidth'', value');
      disp('      ''linestyle'', ( linespec | ''none'' )');
      disp('      ''fontsize'',value');
      disp('      ''fontname'',name');
      disp('      ''XaxisLocation'',( ''bottom'' | ''middle'' | ''top'' ) ');
      disp('      ''YaxisLocation'',( ''left'' | ''middle'' | ''right'' ) ');
      return;
    case 'set',
      disp(['      box = ' gbox]);
      disp(['      xtick = ' num2str(xtick)]);
      disp(['      ytick = ' num2str(ytick)]);
      disp(['      ticklen = ' num2str(gticklen)]);
      disp(['      tickdir = ' gtickdir]);
      disp(['      xlabeldir = ' gxticklabeldir]);
      disp(['      ylabeldir = ' gyticklabeldir]);
      disp(['      color = ' gcolor]);
      disp(['      linewidth = ' num2str(glinewidth)]);
      disp(['      linestyle = ' glinestyle]);
      disp(['      fontsize = ' num2str(gfontsize)]);
      disp(['      fontname = ' gfontname]);
      disp(['      XaxisLocation = ' gxaxisloc]);
      disp(['      YaxisLocation = ' gyaxisloc]);
      return;
  end;
  k=k+2;
end;     


% Draw the plot box

[X,Y]=feval(MAP_PROJECTION.routine,'box');

if strcmp(gbox,'on');
  line(X(:),Y(:),'linest','-','linewi',glinewidth,'color',gcolor,'tag','m_box','clip','off');
end;

% Axes background - to defeat the inverthardcopy, I need a non-white border (the edgecolor),
% but sneakily I can set it's width to (effectively) 0 so it doesn't actually show!
%
% Now, I used to set this at a large (negative) zdata, but this didn't work for PC users,
% so now I just draw a patch
patch('xdata',X(:),'ydata',Y(:),'facecolor',get(gca,'color'),...
      'edgecolor','k','linest','none','tag','m_color');

% Now I set it at the bottom of the children list so it gets drawn first (i.e. doesn't
% cover anything)
 show=get(0, 'ShowHiddenHandles');
 set(0, 'ShowHiddenHandles', 'on');
 hh=get(gca,'children');
 htags = get(hh,'tag');
 k = strmatch('m_color',htags);
 hht = hh;
 hh(k) = [];
 hh = [hh;hht(k)];
 set(gca,'children',hh);
 set(0, 'ShowHiddenHandles', show);


% X-axis labels and grid

if ~isempty(xtick),

 % Tricky thing - if we are drawing a map with the poles, its nasty when the lines get too close
 % together. So we can sort of fudge this by altering MAP_VAR_LIST.lats to be slightly smaller,
 % and then changing it back again later.
 fudge_north='n';fudge_south='n';
 if ~isempty(ytick) & length(ytick)>1,
  if MAP_VAR_LIST.lats(2)==90, 
    fudge_north='y';
    MAP_VAR_LIST.lats(2)=ytick(end);
  end;
  if MAP_VAR_LIST.lats(1)==-90, 
    fudge_south='y';
    MAP_VAR_LIST.lats(1)=ytick(1);
  end;
 end;

 [X,Y,lg,lgI]=feval(MAP_PROJECTION.routine,'xgrid',xtick,gxaxisloc);
 [labs,scl]=m_labels('lon',lg,xlabels);

 % Draw the grid. Every time we draw something, I first reshape the matrices into a long
 % row so that a) it plots faster, and b) all lines are given the same handle (which cuts
 % down on the number of children hanging onto the axes).

 [n,m]=size(X);
 line(reshape([X;NaN+ones(1,m)],(n+1)*m,1),reshape([Y;NaN+ones(1,m)],(n+1)*m,1),...
      'linest',glinestyle,'color',gcolor,'linewidth',0.1,'tag','m_xgrid');

 % Get the tick data
 [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir);

 % Draw ticks if labels are on top or bottom (not if they are in the middle)

 if strcmp(gxticklabeldir,'middle'),
  if lgI==size(X,1) & strcmp(gxaxisloc,'top'),  % Check to see if the projection supports this option.
   vert='bottom';horiz='center';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+90;
  elseif lgI==1 & strcmp(gxaxisloc,'bottom')
   vert='top';horiz='center';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi-90;
  else
   vert='middle';horiz='center';lgIp1=lgI+1;drawticks=0;
   xx=X(lgI,:); yy=Y(lgI,:);rotang=atan2(Y(lgIp1,:)-Y(lgI,:),X(lgIp1,:)-X(lgI,:))*180/pi-90;
  end;
 else
  if lgI==size(X,1) & strcmp(gxaxisloc,'top'),  % Check to see if the projection supports this option.
   vert='middle';horiz='left';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+180;
  elseif lgI==1 & strcmp(gxaxisloc,'bottom')
   vert='middle';;horiz='right';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi;
  else
   vert='top';;horiz='center';lgIp1=lgI+1;drawticks=0;
   xx=X(lgI,:); yy=Y(lgI,:);rotang=atan2(Y(lgIp1,:)-Y(lgI,:),X(lgIp1,:)-X(lgI,:))*180/pi;
  end;
 end;

 if drawticks,
   [n,m]=size(ltx);
   line(reshape([ltx;NaN+ones(1,m)],(n+1)*m,1),reshape([lty;NaN+ones(1,m)],(n+1)*m,1),...
        'linest','-','color',gcolor,'linewidth',glinewidth,'tag','m_xticks-lower','clip','off');
   line(reshape([utx;NaN+ones(1,m)],(n+1)*m,1),reshape([uty;NaN+ones(1,m)],(n+1)*m,1),...
        'linest','-','color',gcolor,'linewidth',glinewidth,'tag','m_xticks-upper','clip','off');
 end;

 % Add the labels! (whew)

 ik=1:size(X,2);

 for k=ik,
   [rotang(k), horizk, vertk] = upright(rotang(k), horiz, vert);
   text(xx(k),yy(k),labs{k},'horizontal',horizk,'vertical',vertk, ...
        'rot',rotang(k),'fontsize',gfontsize*scl(k),'color',gcolor,...
        'tag','m_xticklabel');
 end;

 if fudge_north=='y',
   MAP_VAR_LIST.lats(2)=90;
 end;
 if fudge_south=='y',
   MAP_VAR_LIST.lats(1)=-90;
 end;

end;

if ~isempty(ytick),
 % Y-axis labels and grid

 [X,Y,lt,ltI]=feval(MAP_PROJECTION.routine,'ygrid',ytick,gyaxisloc);
 [labs,scl]=m_labels('lat',lt,ylabels);

 % Draw the grid
 [n,m]=size(X);
 line(reshape([X;NaN+ones(1,m)],(n+1)*m,1),reshape([Y;NaN+ones(1,m)],(n+1)*m,1),...
      'linest',glinestyle,'color',gcolor,'linewidth',0.1,'tag','m_ygrid');

 % Get the tick data
 [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir);

 % Draw ticks if labels are on left or right (not if they are in the middle)
 if strcmp(gyticklabeldir,'end'),
  if ltI==size(X,1) & strcmp(gyaxisloc,'right'),  % Check to see if the projection supports this option.
   horiz='left';vert='middle';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+180;
  elseif ltI==1 & strcmp(gyaxisloc,'left');
   horiz='right';vert='middle';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi;
  else
   horiz='center';vert='top';ltIp1=ltI+1;drawticks=0;
   xx=X(ltI,:); yy=Y(ltI,:);rotang=atan2(Y(ltIp1,:)-Y(ltI,:),X(ltIp1,:)-X(ltI,:))*180/pi;
  end;
 else
  if ltI==size(X,1) & strcmp(gyaxisloc,'right'),  % Check to see if the projection supports this option.
   horiz='center';vert='top';drawticks=1;
   xx=utx(1,:);yy=uty(1,:);rotang=atan2(diff(uty),diff(utx))*180/pi+270;
  elseif ltI==1 & strcmp(gyaxisloc,'left');
   horiz='center';vert='bottom';drawticks=1;
   xx=ltx(1,:);yy=lty(1,:);rotang=atan2(diff(lty),diff(ltx))*180/pi+90;
  else
   horiz='left';vert='middle';ltIp1=ltI+1;drawticks=0;
   xx=X(ltI,:); yy=Y(ltI,:);rotang=atan2(Y(ltIp1,:)-Y(ltI,:),X(ltIp1,:)-X(ltI,:))*180/pi+90;
  end;
 end;

 if drawticks,
   [n,m]=size(ltx);
   line(reshape([ltx;NaN+ones(1,m)],(n+1)*m,1),reshape([lty;NaN+ones(1,m)],(n+1)*m,1),...
        'linest','-','color',gcolor,'linewidth',glinewidth,'tag','m_yticks-left','clip','off');
   line(reshape([utx;NaN+ones(1,m)],(n+1)*m,1),reshape([uty;NaN+ones(1,m)],(n+1)*m,1),...
        'linest','-','color',gcolor,'linewidth',glinewidth,'tag','m_yticks-right','clip','off');
 end;

 % Finally - the labels!
 ik=1:size(X,2);

 for k=ik,
   [rotang(k), horizk, vertk] = upright(rotang(k), horiz, vert);
   text(xx(k),yy(k),labs{k},'horizontal',horizk,'vertical',vertk,...
        'rot',rotang(k),'fontsize',gfontsize*scl(k),'color',gcolor,'tag','m_yticklabels');
 end;

end;

% Give a 1-1 aspect ratio and get rid of the matlab-provided axes stuff.

set(gca,'visible','off',...
        'dataaspectratio',[1 1 1],...
        'xlim',MAP_VAR_LIST.xlims,...
        'ylim',MAP_VAR_LIST.ylims);

set(get(gca,'title'),'visible','on');
set(get(gca,'xlabel'),'visible','on');
set(get(gca,'ylabel'),'visible','on');

%-------------------------------------------------------------
% upright simply turns tick labels right-side up while leaving
% their positions unchanged.
% Sat  98/02/21 Eric Firing
%
function   [rotang, horiz, vert] = upright(rotang, horiz, vert);
if rotang > 180, rotang = rotang - 360; end
if rotang < -180, rotang = rotang + 360; end
if rotang > 90,
   rotang = rotang - 180;
elseif rotang < -90,
   rotang = 180 + rotang;
else
   return    % no change needed.
end
switch horiz(1)
   case 'l'
      horiz = 'right';
   case 'r'
      horiz = 'left';
end
switch vert(1)
   case 't'
      vert = 'bottom';
   case 'b'
      vert = 'top';
end
  

%--------------------------------------------------------------------------
function [L,fs]=m_labels(dir,vals,uservals);
% M_LONLABEL creates longitude labels
%         Default values are calculated automatically when the grid is 
%         generated. However, the user may wish to specify the labels
%         as either numeric values or as strings (in the usual way
%         for axes).
%
%         If auto-labelling occurs, minutes are labelled in a different
%         (smaller) fontsize than even degrees.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

% If the user has specified [] (i.e. no labels), we return blanks.

if isempty(uservals), 
  L=cellstr(char(' '*ones(length(vals),1)));
  fs=1.0*ones(length(L),1);
  return;
end;

% If the user has specified strings, we merely need to make
% sure that there are enough to cover all ticks.

if any(isstr(uservals)), 
  L=cellstr( uservals((rem([0:length(vals)-1],length(uservals))+1),:) );
  fs=1.0*ones(length(L),1);
  return;
end;

% Otherwise we are going to have to generate labels from numeric
% data.

if length(uservals)==1 & isnan(uservals),  % use default values
  vals=vals(:)'; % make it a row.
else                                       % or ones provided
  lv=length(vals);
  vals=uservals(:)';
  while length(vals)<lv,
    vals=[vals uservals(:)'];
  end;
end;

% longitudes and latitudes have some differences....
if findstr(dir,'lat'), 
  labname=['S';'N';' '];
else
  labname=['W';'E';' '];
  vals=rem(vals+540,360)-180;
end;

i=[vals<0;vals>0;vals==0];  % get the 'names' (i.e. N/S or E/W)
vals=abs(vals);             % Convert to +ve values

L=cell(length(vals),1);
fs=ones(length(vals),1);

% For each label we have different options:
%  1 - even degrees are just labelled as such.
%  2 - ticks that fall on even minutes are just labelled as even minutes
%      in a smaller fontsize.
%  3 - fractional minutes are labelled to 2 decimal places in the
%      smaller fontsize.
for k=1:length(vals),
  if rem(vals(k),1)==0,
    nam=find(i(:,k));
    L{k}=sprintf([' %3.0f^o' labname(nam) ' '],vals(k));
  elseif abs(vals*60-round(vals*60))<.01,
    L{k}=sprintf([' %2.0f'' '],rem(vals(k),1)*60);
    fs(k)=0.75;
  else
    L{k}=sprintf([' %2.2f'' '],rem(vals(k),1)*60);
    fs(k)=0.75;
  end;
end;

% In most cases, the map will have at least one tick with an even degree label,
% but for very small regions (<1 degree in size) this won't happen so we
% want to force one label to show degrees *and* minutes.

if ~any(fs==1),  
 k=round(length(vals)/2);
 nam=find(i(:,k));
 L{k}={sprintf([' %3.0f^o' labname(nam) ' '],fix(vals(k))),...
       sprintf([' %2.2f'' '],rem(vals(k),1)*60)};
 fs(k)=1;
end;


%---------------------------------------------------------
function [ltx,lty,utx,uty]=maketicks(X,Y,gticklen,gtickdir);
% MAKETICKS makes the axis ticks.
%           AXes ticks are based on making short lines at
%           the end of the grid lines X,Y.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

global MAP_VAR_LIST

tlen=gticklen*max( diff(MAP_VAR_LIST.xlims),diff(MAP_VAR_LIST.ylims));

lx=sqrt((X(2,:)-X(1,:)).^2+(Y(2,:)-Y(1,:)).^2);

if strcmp(gtickdir,'out'),
  ltx=[X(1,:)-tlen*(X(2,:)-X(1,:))./lx;X(1,:)];
  lty=[Y(1,:)-tlen*(Y(2,:)-Y(1,:))./lx;Y(1,:)];
else
  ltx=[X(1,:);X(1,:)+tlen*(X(2,:)-X(1,:))./lx];
  lty=[Y(1,:);Y(1,:)+tlen*(Y(2,:)-Y(1,:))./lx];
end;

lx=sqrt((X(end,:)-X(end-1,:)).^2+(Y(end,:)-Y(end-1,:)).^2);

if strcmp(gtickdir,'out'),
  utx=[X(end,:)-tlen*(X(end-1,:)-X(end,:))./lx;X(end,:)];
  uty=[Y(end,:)-tlen*(Y(end-1,:)-Y(end,:))./lx;Y(end,:)];
else
  utx=[X(end,:);X(end,:)+tlen*(X(end-1,:)-X(end,:))./lx];
  uty=[Y(end,:);Y(end,:)+tlen*(Y(end-1,:)-Y(end,:))./lx];
end;





