function [cs,h]=m_contour(long,lat,varargin);
%  M_CONTOUR Draws contour lines on a map
%    M_CONTOUR(LONG,LAT,DATA,...) draw contours on a map. Behavior
%    is the same as for CONTOUR except that LONG and LAT vectors or
%    matrices must be specified.
%
%    [CS,H]=M_CONTOUR(...) returns a contour matrix C and a vector
%    H of handles to LINE or PATCH objects for use by CLABEL.
%
%    See also CONTOUR

% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if min(size(long))==1 & min(size(lat))==1,
 [long,lat]=meshgrid(long,lat);
end;

[X,Y]=m_ll2xy(long,lat,'clip','on');
[cs,h]=contour(X,Y,varargin{:});
