function [varargout] = pcontinents
% PCONTINENTS 	Plots continental outlines.

  % load topography of the Earth (from topo.dat) 
  load topo
  
  % wrap longitude axis such that Greenwhich meridian is at 0 and 360
  topo      = [topo(:, 181:360), topo(:, 1:180), topo(:, 181:360), topo(:, 1:180)];
  
  % plot zero contour of topography
  hold on
  [cs,h] = contour(-179:540, -89:90, topo, [0 0], 'k-', 'linewidth',1);
  set(h, 'LineWidth',0.5);
  grid off
  hold off
  
  if nargout == 1
    varargout{1} = h;
  elseif nargout == 2
     varargout{1} = cs;
     varargout{2} = h;
  end