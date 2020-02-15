% MP_SOLVE_COUPLING2D: Solve the 2d coupled magnetoelastic (static) problem in a multipatch domain.
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
%    - f_mag:                magnetic source term (current density)
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

function [geometry_mec, msh_mec, space_mec, msh_mag, space_mag] = mp_rebuild_coupling2d (problem_data, method_data)
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
end
