% INPUT:
%
%     A:           vector of dof weights
%     space:       object defining the discrete space (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     npts:        cell array with coordinates of points along each parametric direction

function plot_curl2d_dat (A, space, geometry, npts, filename, iptc)

ndim = numel (space.knots);

if (iscell (npts))
    vtk_pts = npts;
    npts = cellfun (@numel, vtk_pts);
elseif (numel (npts) == 1)
    npts = npts * ones (1, ndim);
end

if (~exist ('vtk_pts', 'var'))
  for idim = 1:ndim
    vtk_pts{idim} = linspace (space.knots{idim}(1), space.knots{idim}(end), npts(idim));
  end
end

[eA, F] = sp_eval (A, space, geometry, vtk_pts, 'gradient');
eA = squeeze(sqrt( eA(1,:,:).^2 + eA(2,:,:).^2 ));
rdim = size (F, 1);

if (ndim == 2)
  if (rdim == 2)
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    surf (X, Y, eA)
    % write_dat3D([filename '_' num2str(iptc) '.dat'], X, Y, eA);
  end
end
end
