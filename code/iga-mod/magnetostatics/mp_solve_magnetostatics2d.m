% MP_SOLVE_MAGNETOSTATICS2D: Solve the 2d magnetostatic problem in a multipatch geometry.
%
% Function to solve the weak problem
%
%     S (1/mu22 * A_x1 * At_x1 + 1/mu11 * A_x2 * At_x2) dx = S (f) dx    in Omega
%                                                    A x n = g           on Gamma_N
%                                                        u = h           on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_mag:          name of the file containing the geometry
%    - nmnn_sides_mag:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides_mag: sides with Dirichlet boundary condition
%    - mu:               magnetic permeability as cell array of function handles (diagonal components of material tensor)
%    - f_mag:            source term (current density), one dimensional
%    - g_mag:            function for Neumann condition
%    - h_mag:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: array of geometry structures (see geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space:    multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch)
%  A:        the computed degrees of freedom

function [geometry, msh, space, A] = ...
              mp_solve_magnetostatics2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_mag);
npatch = numel (geometry);

msh = cell (1, npatch);
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

msh   = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp

% Compute and assemble the matrices
mat = op_ms2d_mp (space, space, msh, mu);
rhs = op_f_v_mp_mod (space, msh, f_mag);

% Apply Neumann boundary conditions
Nbnd = cumsum ([0, boundaries.nsides]);
for iref = nmnn_sides_mag
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
  gref = @(varargin) g_mag(varargin{:},iref);
  rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
  rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn;
end

% Apply Dirichlet boundary conditions
A = zeros (space.ndof, 1);
[A_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_mod (space, msh, h_mag, drchlt_sides_mag);
A(drchlt_dofs) = A_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs)*A_drchlt;

% Solve the linear system
A(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
