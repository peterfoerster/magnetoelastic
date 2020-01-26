% SOLVE_LINEAR_ELASTICITY2D: Solve a 2d linear elasticity problem on a NURBS domain.
%
% The function solves the linear elasticity problem
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
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - E:            cell array of function handles for Young's modulus
%    - nu:           Poisson's ratio
%    - G:            pre-magnetization
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
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
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete basis functions (see sp_vector)
%  u:        the computed degrees of freedom

function [geometry, msh, sp, u] = ...
              solve_linear_elasticity2d (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry = mp_geo_load (geo_name);
% only take specific patch into account
geometry = geometry(iptc);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
space_scalar = sp_nurbs (nurbs, msh);
scalar_spaces = repmat ({space_scalar}, 1, msh.rdim);
sp = sp_vector (scalar_spaces, msh);
clear space_scalar scalar_spaces

% Assemble the matrices
mat    = op_linear_elasticity2d_tp (sp, sp, msh, E, nu_mech, G, iptc);
rhs    = op_f_v_tp (sp, msh, f_mech);

% Apply Neumann boundary conditions
for iside = nmnn_sides_mech
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
  gside = @(varargin) g_mech(varargin{:},iside);
  dofs = sp.boundary(iside).dofs;
  rhs(dofs) = rhs(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), gside);
end

% Apply Dirichlet boundary conditions
u = zeros (sp.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h_mech, drchlt_sides_mech);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - mat (int_dofs, drchlt_dofs) * u_drchlt;

% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
