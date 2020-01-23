% INPUT:
%
%     u:           vector of dof weights
%     space:       object representing the space of discrete functions (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     npts:        number of points along each parametric direction where to evaluate
%     pts:         cell array with the coordinates along each parametric direction of the points where to evaluate
%     filename:    name of the output file.
%     fieldnames:  how to name the saved variables in the vtk file
%
% OUTPUT:
%
%    none

function sp_to_vtk_curl2d (u, space, geometry, npts, filename, fieldname)

  % eu = [u_x1; u_x2]
  [eu, F] = sp_eval (u, space, geometry, npts, 'gradient');
  % curl(u), if only a z-component is present
  b = [eu(2,:,:); -eu(1,:,:)];

  msh_to_vtk (F, b, filename, fieldname);

end
