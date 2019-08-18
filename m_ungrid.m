function m_ungrid;
% M_UNGRID Removes a grid;
%          M_UNGRID deletes a map grid, but leaves any plotted
%          data

% Rich Pawlowicz (rich@ocgy.ubc.ca) 4/Apr/97 
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

hh=get(gca,'children');

things=get(hh,'userdata');

for i=1:length(hh),
  if ~isempty(things{i}) & strmatch('m_',things{i}),
   delete(hh(i));
 end;
end;

set(gca,'visible','on');
