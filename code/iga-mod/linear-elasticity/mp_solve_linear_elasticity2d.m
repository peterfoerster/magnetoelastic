% MP_SOLVE_LINEAR_ELASTICITY2D: Solve a 2d linear elasticity problem in a multipatch domain.
%
% Function to solve the problem
%
%     S ( 1/(E1-nu12^2*E2)(E1^2*u1_x1*v1_x1 +
%         nu12*E2*E1*(u1_x1*v2_x2 + u2_x2*v1_x1) + E1*E2*u2_x2*v2_x2 +
%         G12*(E1-nu12^2*E2)*(u1_x2 + u2_x1)*(v1_x2 + v2_x1)) ) dx = S (f) dx   in Omega
%                                                                ? = g          on Gamma_N
%                                                                u = h          on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_mec:          name of the file containing the geometry
%    - nmnn_sides_mec:   sides with Neumann boundary condition
%    - drchlt_sides_mec: sides with Dirichlet boundary condition
%    - E:                cell array of function handles for Young's modulus
%    - nu:               Poisson's ratio
%    - G:                pre-magnetization
%    - f_mec:            source term
%    - h_mec:            function for Dirichlet boundary condition
%    - g_mec:            function for Neumann condition (if nmnn_sides is not empty)
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

function [geometry, msh, space, u] = ...
              mp_solve_linear_elasticity2d (problem_data, method_data)

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
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_mec);
npatch = numel (geometry);

for iptc = 1:npatch
  degelev  = max (degree - (geometry(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geometry(iptc).nurbs, degelev);
  [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

  nurbs    = geo_load (nrbkntins (nurbs, nknots));
  fields   = fieldnames (nurbs);
  for ifld = 1:numel (fields)
      geometry(iptc).(fields{ifld}) = nurbs.(fields{ifld});
  end

% Construct msh structure
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (geometry(iptc).nurbs.knots, rule);
  msh{iptc} = msh_cartesian (geometry(iptc).nurbs.knots, qn, qw, geometry(iptc));

% Construct space structure
  sp_scalar = sp_nurbs (geometry(iptc).nurbs, msh{iptc});
  scalar_spaces = cell (msh{iptc}.rdim, 1);
  for idim = 1:msh{iptc}.rdim
    scalar_spaces{idim} = sp_scalar;
  end
  sp{iptc} = sp_vector (scalar_spaces, msh{iptc});
end

msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp scalar_spaces

% Compute and assemble the matrices (filter NaN and fill with zeros?)
mat = op_le2d_mp (space, space, msh, E, nu, G);
rhs = op_f_v_mp (space, msh, f_mec);

% Apply Neumann boundary conditions
Nbnd = cumsum ([0, boundaries.nsides]);
for iref = nmnn_sides_mec
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
  gref = @(varargin) g_mec(varargin{:},iref);
  rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
  rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn;
end

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_mod (space, msh, h_mec, drchlt_sides_mec);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
end
