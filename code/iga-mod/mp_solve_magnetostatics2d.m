% MP_SOLVE_MAGNETOSTATICS_2D: Solve the 2d magnetostatic problem in a multipatch geometry.
%
% Function to solve the weak problem
%
%     S (nu11*A_x2*At_x2 + nu22*A_x1*At_x1) dx = S (f) dx    in Omega
%                                            ? = g           on Gamma_N
%                                            u = h           on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - nu:           magnetic reluctivity as cell array of fuction handles (diagonal components of material tensor)
%    - f:            source term, one dimensional
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
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
%  u:        the computed degrees of freedom
%

function [geometry, msh, space, u] = ...
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
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
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

msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp

% Compute and assemble the matrices
stiff_mat = op_gradu_gradv_mp_mstatic (space, space, msh, nu_mag);
rhs = op_f_v_mp_mod (space, msh, f);

% Apply Neumann boundary conditions
Nbnd = cumsum ([0, boundaries.nsides]);
for iref = nmnn_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
  gref = @(varargin) g(varargin{:},iref);
  rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
  rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn;
end

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
