function m_scale(scale_factor),
% M_SCALE Draws the map at a specified scale
%         After M_GRID has been called, the map is sized to fit within
%         the plot region of the current axes. If you want to force the
%         map to appear at a specific scale (e.g. 1:250000), call
%
%         M_SCALE(scale_factor)
%
%         where the map scale is 1:scale_factor. The map will be drawn
%         with its origin at bottom left within the figure limits (which
%         can be set with calls to ORIENT or by setting the 'paperposition'
%         property of the figure).
%
%         M_SCALE('auto') returns to autoscaling.
%
%         see also M_PROJ, M_GRID.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 13/Nov/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% Need earth radius, in centimeters.
erad=637813700; %cm (from WGS-84)

if isstr(scale_factor),
 gca_pos=get(gca,'userdata');
 if ~isempty(gca_pos),
   set(gca,'position',gca_pos);
 end;
else
  
 fu=get(gcf,'units');
 set(gcf,'units','centimeter');
 f_pos=get(gcf,'position');
 set(gcf,'units',fu);

 map_x=diff(get(gca,'xlim'));
 map_y=diff(get(gca,'xlim'));

 if map_x>10, % we are probably using meters (i.e. UTM coords)
   map_x=map_x/scale_factor*100;
   map_y=map_y/scale_factor*100;
 else
   map_x=map_x*erad/scale_factor;
   map_y=map_y*erad/scale_factor;
 end;
 
 
 if map_x> f_pos(3) | map_y>f_pos(4),
   disp('Warning - map larger than current window at this scale');
 elseif map_x< f_pos(3)/2 & map_y<f_pos(4)/2,
   disp('Warning - map much smaller than current window at this scale');
 end;

 if isempty(get(gca,'userdata')),
  set(gca,'userdata',get(gca,'position'));
 end;

 gca_u=get(gca,'unit');
 set(gca,'unit','centi','position',[1 1 map_x map_y],'unit',gca_u)

end;




