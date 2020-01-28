% INPUT:
%
%     A:           vector of dof weights
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

function sp_to_vtk_curl2d (A, space, geometry, npts, filename, fieldname)

  [eA, F] = sp_eval (A, space, geometry, npts, 'gradient');
  % curl(A), if only a z-component is present
  b = [eA(2,:,:); -eA(1,:,:)];

  msh_to_vtk (F, b, filename, fieldname);

end
