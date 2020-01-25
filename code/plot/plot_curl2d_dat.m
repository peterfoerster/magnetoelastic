% INPUT:
%
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     npts:        cell array with coordinates of points along each parametric direction

function plot_curl2d_dat (u, space, geometry, npts, filename, iptc)

ndim = numel (space.knots);

if (nargin < 4 || isempty (npts))
  npts = 51 * ones (1, ndim);
else
  if (iscell (npts))
    vtk_pts = npts;
    npts = cellfun (@numel, vtk_pts);
  elseif (numel (npts) == 1)
    npts = npts * ones (1, ndim);
  end
end

if (ndim == 3)
  ncuts = 2 * ones (1, ndim);
end

if (~exist ('vtk_pts', 'var'))
  for idim = 1:ndim
    vtk_pts{idim} = linspace (space.knots{idim}(1), space.knots{idim}(end), npts(idim));
  end
end

[eu, F] = sp_eval (u, space, geometry, vtk_pts, 'gradient');
eu = squeeze(sqrt( eu(1,:,:).^2 + eu(2,:,:).^2 ));
rdim = size (F, 1);

if (ndim == 1)
  if (rdim == 1)
    plot (F, eu)
  elseif (rdim == 2)
    plot3 (F(1,:), F(2,:), eu)
  else
    error ('The plot on 3D curves has not been implemented yet')
  end
elseif (ndim == 2)
  if (rdim == 2)
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    surf (X, Y, eu)
    % write_dat3D([filename '_' num2str(iptc) '.dat'], X, Y, eu);
  elseif (rdim == 3)
    [X, Y, Z]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)));
    surf (X, Y, Z, eu)
  end
elseif (ndim == 3)
  hold_flag = ishold;

  for idim = 1:ndim
    plot_pts = vtk_pts;
    plot_pts{idim} = linspace(space.knots{idim}(1), space.knots{idim}(end), ncuts(idim)+2);
    [eu, F] = sp_eval (u, space, geometry, plot_pts);
    indices = {1:npts(1), 1:npts(2), 1:npts(3)};

    if (ncuts(idim) > 0)
      cuts = 2:ncuts(idim)+1;
    else
      cuts = 1:ncuts(idim)+2;
    end
    for ii = cuts
      indices{idim} = ii;
      surf (squeeze(F(1,indices{:})), squeeze(F(2,indices{:})), squeeze(F(3,indices{:})), squeeze(eu(indices{:})));
      hold on
    end
  end

  if (~hold_flag)
    hold off;
  end
end

end
