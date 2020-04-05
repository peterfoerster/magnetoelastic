% SOLVE_WEAK_COUPLING1D: Solve a weakly coupled magnetoelastic problem on a one-dimensional domain.
%
%  Function to solve the weak problem
%
%     <rho * w'',v> + a(w,v) = l(v)     in Omega = F(0,1)
%
%     with w1(0,t) = w2(0,t) = 0   and   w2'(0,t) = 0,
%
%  where <rho * w'',v> = b^2 S rho( v1 * w1'' + v2 * w2'' ) dx1,
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
%  msh1/2:   mesh objects that define the quadrature rules (see msh_2d)
%  space1/2: space objects that define the discrete basis functions (see sp_scalar)
%  w:        the computed degrees of freedom

function [geometry, msh1, space1, msh2, space2, w] = solve_weak_coupling1d (problem_data, method_data)

   data_names = fieldnames (problem_data);
   for iopt=1:numel(data_names)
      eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
   end
   data_names = fieldnames (method_data);
   for iopt=1:numel(data_names)
      eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
   end

   geometry = geo_load (geo_name);

   % two one-dimensional meshes
   [knots1, zeta] = kntrefine (geometry.nurbs.knots, nsub(1)-1, degree(1), regularity(1));
   rule     = msh_gauss_nodes (nquad(1));
   [qn, qw] = msh_set_quad_nodes (zeta, rule);
   msh1     = msh_cartesian (zeta, qn, qw, geometry, 'boundary', true, 'der2', true);

   [knots2, zeta] = kntrefine (geometry.nurbs.knots, nsub(2)-1, degree(2), regularity(2));
   rule     = msh_gauss_nodes (nquad(2));
   [qn, qw] = msh_set_quad_nodes (zeta, rule);
   msh2     = msh_cartesian (zeta, qn, qw, geometry, 'boundary', true, 'der2', true);

   % two one-dimensional scalar spaces
   space1 = sp_bspline (knots1, degree(1), msh1);
   space2 = sp_bspline (knots2, degree(2), msh2);

   % assemble matrices
   mat1 = op_wmec1d_w1_tp (space1, space1, msh1, b, rho, E, A, I);
   mat2 = op_wmec1d_w2_tp (space2, space2, msh2, b, rho, E, A, I);

   % assemble rhs
   rhs1 = op_l_v1_tp (space1, msh1, D, B, A);
   rhs2 = op_l_v2_tp (space2, msh2, D, B, b);

   % apply Dirichlet boundary conditions
   drchlt1_dofs = [];
   for iside=drchlt1_sides
      drchlt1_dofs = [drchlt1_dofs, space1.boundary(iside).dofs];
   end

   % drchlt2_dofs = [];
   % for iside=drchlt2_sides
   %    drchlt2_dofs = [drchlt2_dofs, space2.boundary(iside).dofs];
   % end
   % int_dofs2 = setdiff(1:space2.ndof, drchlt2_dofs);

   % apply secondary Neumann condition
   % for iside=nmnn2_sides
   %    msh2_side   = msh_boundary_side_from_interior (msh2, iside);
   %    space2_side = sp_precompute (space2.constructor(msh2_side), msh2_side, 'gradient', true);
   %    rhs2(space2_side.connectivity) = rhs2(space2_side.connectivity) - g(iside) * reshape (space2_side.shape_function_gradients(:,:,:), space2_side.nsh, 1);
   % end

   % apply secondary Dirichlet condition
   nL = numel(drchlt2_sides);
   L  = sparse (nL, space2.ndof);
   for iside=drchlt2_sides
      msh2_side   = msh_boundary_side_from_interior (msh2, iside);
      space2_side = sp_precompute (space2.constructor (msh2_side), msh2_side, 'gradient', true);
      L(iside, space2_side.connectivity) = reshape(space2_side.shape_function_gradients(:,:,:), 1, space2_side.nsh);
   end
   mat2 = [mat2, L(drchlt2_sides,:).';
           L(drchlt2_sides,:), sparse(nL,nL)];
   rhs2(space2.ndof+(1:nL)) = 0;
   int_dofs2 = 1:(space2.ndof+nL);

   mat = [mat1, sparse(size(mat1,1), size(mat2,2));
          sparse(size(mat2,1), size(mat1,2)), mat2];
   rhs  = [rhs1; rhs2];

   % solve the system
   w = zeros(size(mat,2),1);
   int_dofs1 = setdiff (1:space1.ndof, drchlt1_dofs);
   int_dofs  = [int_dofs1, space1.ndof+int_dofs2];

   mat_dofs = mat(int_dofs, int_dofs);
   rhs_dofs = rhs(int_dofs);
   w(int_dofs) = mat_dofs\rhs_dofs;
   w = w(1:(space1.ndof+space2.ndof));
end
