function nrbplot_color (color, width, nurbs, subd)
%
% NRBPLOT: Plot a NURBS curve or surface, or the boundary of a NURBS volume.
%
% INPUT:
%
%   nrb		: NURBS curve, surface or volume, see nrbmak.
%
%   npnts	: Number of evaluation points, for a surface or volume, a row
%       vector with the number of points along each direction.
%

nargs = nargin;
if nargs < 4
  error ('Need a NURBS to plot and the number of subdivisions!');
end

% Default values
light='off';
cmap='summer';
colormap (cmap);

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
  knt = nurbs.knots;
  order = nurbs.order;
  p = nrbeval (nurbs, {linspace(knt{1}(order(1)),knt{1}(end-order(1)+1),subd(1)) ...
                       linspace(knt{2}(order(2)),knt{2}(end-order(2)+1),subd(2))});
  if (strcmp (light,'on'))
    % light surface
    surfl (squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)));
    shading interp;
  else
    surf (squeeze (p(1,:,:)), squeeze (p(2,:,:)), squeeze (p(3,:,:)));
    shading faceted;
  end
 elseif (size (nurbs.knots,2) == 3) % plot the boundaries of a NURBS volume
  bnd = nrbextract (nurbs);
  hold_flag = ishold;
  nrbplot (bnd(1), subd(2:3), varargin{:});
  hold on
  nrbplot (bnd(2), subd(2:3), varargin{:});
  nrbplot (bnd(3), subd([1 3]), varargin{:});
  nrbplot (bnd(4), subd([1 3]), varargin{:});
  nrbplot (bnd(5), subd(1:2), varargin{:});
  nrbplot (bnd(6), subd(1:2), varargin{:});

  if (~hold_flag)
    hold off
  end

 else
  error ('nrbplot: some argument is not correct')
 end
else
  % plot a NURBS curve
  order = nurbs.order;
  p = nrbeval (nurbs, linspace (nurbs.knots(order), nurbs.knots(end-order+1), subd));

  if (any (nurbs.coefs(3,:)))
    % 3D curve
    plot3 (p(1,:), p(2,:), p(3,:));
    grid on;
  else
    % 2D curve
    plot (p(1,:), p(2,:), 'color', color, 'linewidth', width);
  end
end
axis equal;

end
