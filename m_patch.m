function h=m_patch(long,lat,varargin);
% M_PATCH Create patches on a map
%    M_PATCH(LONG,LAT,C) is a drop-in replacement for PATCH that uses 
%    longitude/latitude coordinates to draw a patch on the current map. 
%    See PATCH for more details about the way in which patch colours and 
%    properties should be specified.
%
%    Currently you cannot specify C to be other than a string or 1x3 RGB
%    vector.
%
%    See also M_LINE, M_LL2XY

% Rich Pawlowicz (rich@ocgy.ubc.ca) 3/Sep/98
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

[m,n]=size(long);

if m==1 & n>1,
  h=mu_coast('vector',[long' lat';long(1) lat(1)],'patch',varargin{:});
elseif m>1 & n==1,
  h=mu_coast('vector',[long lat;long(1) lat(1)],'patch',varargin{:});
else
  h=mu_coast('vector',[reshape([long;long(1,:);NaN+ones(1,n)],(m+2)*n,1),...
                     reshape([lat;lat(1,:);NaN+ones(1,n)],(m+2)*n,1)],'patch',varargin{:});
end;
