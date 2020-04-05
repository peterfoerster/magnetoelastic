% INPUT:
%
%  space: object representing the multipatch space of trial functions (see sp_multipatch)
%  msh:   object containing the domain partition and the quadrature rule (see msh_multipatch)
%  h:     function handle to compute the Dirichlet condition
%  refs:  boundary references on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions

function [u, dofs] = sp_drchlt_l2_proj_mod (space, msh, h, refs, varargin)

   M = spalloc (space.boundary.ndof, space.boundary.ndof, 3*space.boundary.ndof);
   rhs = zeros(space.boundary.ndof,1);

   boundaries = msh.boundaries;
   Nbnd = cumsum ([0, boundaries.nsides]);
   bnd_dofs = [];
   for iref=refs
      iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
      href = @(varargin) h(varargin{:}, iref);
      f_one = @(varargin) ones (size(varargin{1}));

      M = M + op_u_v_mp (space.boundary, space.boundary, msh.boundary, f_one, iref_patch_list);
      rhs = rhs + op_f_v_mp (space.boundary, msh.boundary, href, iref_patch_list);

      boundary_gnum = space.boundary.gnum;

      % apply union to individual bnd_dofs
      bnd_refs = [];
      for irpl=iref_patch_list
         bnd_refs = union (bnd_refs, boundary_gnum{irpl});
      end

      bnd_dofs = union (bnd_dofs, boundary_dofs);
   end

   u = M(bnd_dofs,bnd_dofs) \ rhs(bnd_dofs);
   dofs = space.boundary.dofs(bnd_dofs);

   if (~isempty (space.boundary.boundary_orientation))
      u = u .* space.boundary.boundary_orientation(bnd_dofs).';
   end
end
