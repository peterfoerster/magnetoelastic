% MP_SOLVE_COUPLING2D_FD: Solve the 2d coupled magnetoelastic (quasistatic, frequency domain) problem in a multipatch domain.
%
% Function to solve the problem (check matrix formulation)
%
%       matrix:      DoFs:    rhs:
%     (A     -C)    (u)      (0)
%     (-C^T   B)    (A)      (j)
%
% where A = a(u,v), B = b(A,At), C = c(v,A) and j is the discrete current density.
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_mag/mec:          name of the file containing the geometry
%    - nmnn_sides_mag/mec:   sides with Neumann boundary condition
%    - drchlt_sides_mag/mec: sides with Dirichlet boundary condition
%    - mu:                   magnetic permeability
%    - kappa:                complex conductivity
%    - f_mag:                magnetic source term (current density phasor)
%    - g_mag:                function for Neumann condition
%    - h_mag:                function for Dirichlet boundary condition
%    - E:                    Young's modulus
%    - nu:                   Poisson's ratio
%    - G:                    pre-magnetization
%    - tau:                  coefficient for Neumann condition
%    - f_mech:               source term
%    - g_mech:               function for Neumann condition
%    - h_mech:               function for Dirichlet boundary condition
%    - f:                    coupling coefficients
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
%  geometry:      array of geometry structures (see geo_load)
%  msh_mec/mag:   multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space_mec/mag: multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch)
%  u/A:           the computed degrees of freedom

function [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d_fd (problem_data, method_data)
data_names=fieldnames(problem_data);
for iopt=1:numel(data_names)
   eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names=fieldnames(method_data);
for iopt=1:numel(data_names)
   eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

[geometry_mag, boundaries_mag, interfaces_mag, ~, boundary_interfaces_mag] = mp_geo_load (geo_mag);
[geometry_mec, boundaries_mec, interfaces_mec, ~, boundary_interfaces_mec] = mp_geo_load (geo_mec);
npatch = numel(geometry_mec);

% construct spaces
msh     = cell(1,npatch);
sp_mag  = cell(1,npatch);
sp_mec  = cell(1,npatch);
for iptc=1:npatch
   degelev = max(degree - (geometry_mec(iptc).nurbs.order-1), 0);
   nurbs   = nrbdegelev (geometry_mec(iptc).nurbs, degelev);
   [~, ~, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
   nurbs  = geo_load (nrbkntins (nurbs, nknots));
   fields = fieldnames(nurbs);
   for ifld=1:numel(fields)
      geometry_mec(iptc).(fields{ifld}) = nurbs.(fields{ifld});
   end

   rule      = msh_gauss_nodes (nquad);
   [qn, qw]  = msh_set_quad_nodes (geometry_mec(iptc).nurbs.knots, rule);
   msh{iptc} = msh_cartesian (geometry_mec(iptc).nurbs.knots, qn, qw, geometry_mec(iptc));

   % magnetic scalar space (grad-preserving)
   sp_mag{iptc} = sp_bspline (geometry_mec(iptc).nurbs.knots, geometry_mec(iptc).nurbs.order-1, msh{iptc});

   % mechanic vector space (grad-preserving)
   sp_scalar = sp_nurbs (geometry_mec(iptc).nurbs, msh{iptc});
   scalar_spaces = cell(msh{iptc}.rdim,1);
   for idim=1:msh{iptc}.rdim
      scalar_spaces{idim} = sp_scalar;
   end
   sp_mec{iptc} = sp_vector (scalar_spaces, msh{iptc});
end

msh_mag   = msh_multipatch (msh, boundaries_mag);
space_mag = sp_multipatch (sp_mag, msh_mag, interfaces_mag, boundary_interfaces_mag);
msh_mec   = msh_multipatch (msh, boundaries_mec);
space_mec = sp_multipatch (sp_mec, msh_mec, interfaces_mec, boundary_interfaces_mec);
clear msh sp_mag sp_mec sp scalar_spaces

% magnetic problem
B_mat = op_mqs2d_mp (space_mag, space_mag, msh_mag, mu, kappa);
rhs_mag = op_f_v_mp_mod (space_mag, msh_mag, f_mag);

% boundary conditions
Nbnd = cumsum([0, boundaries_mag.nsides]);
for iref=nmnn_sides_mag
   iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
   gref = @(varargin) g_mag(varargin{:},iref);
   rhs_mag_nmnn = op_f_v_mp (space_mag.boundary, msh_mag.boundary, gref, iref_patch_list);
   rhs_mag(space_mag.boundary.dofs) = rhs_mag(space_mag.boundary.dofs) + rhs_mag_nmnn;
end
A = zeros(space_mag.ndof,1);
[A_drchlt, drchlt_dofs_mag] = sp_drchlt_l2_proj (space_mag, msh_mag, h_mag, drchlt_sides_mag);
A(drchlt_dofs_mag) = A_drchlt;
int_dofs_mag = setdiff(1:space_mag.ndof, drchlt_dofs_mag);
rhs_mag(int_dofs_mag) = rhs_mag(int_dofs_mag) - B_mat(int_dofs_mag,drchlt_dofs_mag) * A_drchlt;

% mechanic problem
A_mat = op_le2d_mp (space_mec, space_mec, msh_mec, E, nu ,G);
rhs_mec = op_f_v_mp (space_mec, msh_mec, f_mec);

% boundary conditions
Nbnd = cumsum([0, boundaries_mec.nsides]);
for iref=nmnn_sides_mec
   iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
   gref = @(varargin) g_mec(varargin{:},iref);
   rhs_mec_nmnn = op_f_v_mp (space_mec.boundary, msh_mec.boundary, gref, iref_patch_list);
   rhs_mec(space_mec.boundary.dofs) = rhs_mec(space_mec.boundary.dofs) + rhs_mec_nmnn;
end
u = zeros(space_mec.ndof,1);
[u_drchlt, drchlt_dofs_mec] = sp_drchlt_l2_proj (space_mec, msh_mec, h_mec, drchlt_sides_mec);
u(drchlt_dofs_mec) = u_drchlt;
int_dofs_mec = setdiff(1:space_mec.ndof, drchlt_dofs_mec);
rhs_mec(int_dofs_mec) = rhs_mec(int_dofs_mec) - A_mat(int_dofs_mec,drchlt_dofs_mec) * u_drchlt;

% coupling
C_mat = op_mec2d_mp (space_mec, space_mag, msh_mag, f);

A_mat = A_mat(int_dofs_mec,int_dofs_mec);
B_mat = B_mat(int_dofs_mag,int_dofs_mag);
C_mat = C_mat(int_dofs_mec,int_dofs_mag);

rhs_mag = rhs_mag(int_dofs_mag);
rhs_mec = rhs_mec(int_dofs_mec);

mat = [A_mat -C_mat; -C_mat.' B_mat];
rhs = [rhs_mec; rhs_mag];

DoFs = mat \ rhs;
u(int_dofs_mec) = DoFs(1:length(int_dofs_mec));
A(int_dofs_mag) = DoFs(length(int_dofs_mec)+1:length(int_dofs_mec)+length(int_dofs_mag));
end
