% SOLVE_WEAK_COUPLING1D: Solve a weak magnetoelastic problem on a one-dimensional domain.
%
%  Function to solve the weak problem
%
%     <rho * w'',v> + a(w,v) = l(v)     in Omega = F(0,1)
%
%     with w1(0,t) = w2(0,t) = 0   and   w2'(0,t) = 0,
%
%  where <rho * w'',v> = b^2 S (rho v * w'') dx1,
%                a(w,v) = 1/2*E*A * S (w1' * v1') dx1 + 1/2*E*I * S (w2'' * v2'') dx1,
%                <l,v>  = sum(D1k * Bk) * (A * S (v1') dx1 - b^3/2 * S (v2'') dx1).
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - drchlt_sides: sides with Dirichlet boundary conditions
%    - nmmn_sides:   sides with Neumann boundary conditions
%    - b:            width
%    - rho:          material density
%    - E:            Young's modulus
%    - I:            axial moment of inertia
%    - A:            cross-sectional area
%    - B:            magnetic flux density
%    - g:            function for Neumann boundary conditions
%    - h:            function for Dirichlet boundary conditions
%
%  method_data: a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_2d)
%  space:    space object that defines the discrete basis functions (see sp_vector_2d)
%  w:        the computed degrees of freedom

function [geometry, msh, sp, w] = solve_weak_coupling1d (problem_data, method_data)

   data_names = fieldnames (problem_data);
   for iopt=1:numel(data_names)
      eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
   end
   data_names = fieldnames (method_data);
   for iopt=1:numel(data_names)
      eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
   end

   geometry = geo_load (geo_name);
   [knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
   rule     = msh_gauss_nodes (nquad);
   [qn, qw] = msh_set_quad_nodes (zeta, rule);
   msh      = msh_cartesian (zeta, qn, qw, geometry);
   % one two-dimensional scalar space
   sp = sp_bspline (knots, degree, msh);

   % assemble matrix
   mat = op_wmec1d_tp (sp, sp, msh, b, rho, E, A, I);
keyboard
return
   % assemble rhs
   rhs = op_f_v_tp (sp, msh, f);

% Apply homogeneous 1st Dirichlet boundary conditions
%  and 1st Neumann boundary conditions
drchlt_dofs = [];
for iside = 1:2*msh.ndim
  if (drchlt1_ends(iside))
    drchlt_dofs = [drchlt_dofs, sp.boundary(iside).dofs];
  elseif (nmnn1_ends(iside))
    rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + P;
  end
end

% Apply 2nd Dirichlet boundary conditions by using the Lagrange multipliers method
%  and 2nd Neumann boundary conditions
n_d2 = sum (drchlt2_ends);
C = sparse (numel(drchlt2_ends), sp.ndof);
for iside = 1:2*msh.ndim
  if (drchlt2_ends(iside) || nmnn2_ends(iside))
    msh_aux = msh_boundary_side_from_interior (msh, iside);
    sp_side = sp_precompute (sp.constructor (msh_aux), msh_aux, 'gradient', true);

    if (drchlt2_ends(iside))
      C(iside, sp_side.connectivity) = reshape (sp_side.shape_function_gradients(:,:,:), 1, sp_side.nsh);
    elseif (nmnn2_ends(iside))
      rhs(sp_side.connectivity) = rhs(sp_side.connectivity) - ...
          g(iside) * reshape (sp_side.shape_function_gradients(:,:,:), sp_side.nsh, 1);
    end
  end
end
stiff_mat = [stiff_mat,         C(drchlt2_ends, :).'; ...
             C(drchlt2_ends,:), sparse(n_d2, n_d2)];
rhs(sp.ndof+(1:n_d2)) = 0;
u = zeros (sp.ndof + n_d2, 1);
int_dofs = setdiff (1:(sp.ndof+n_d2), drchlt_dofs);

% Solve the static problem
K = stiff_mat(int_dofs, int_dofs);
F = rhs(int_dofs);
u(int_dofs) = K\F;
u = u(1:sp.ndof);

end
