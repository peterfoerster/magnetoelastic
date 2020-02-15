% INPUT:
%
%   space:      object representing the function space (see sp_multipatch)
%   msh:        object defining the domain partition and the quadrature rule (see msh_multipatch)
%   f:          function handle to compute the source function
%   patch_list: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   rhs: assembled right-hand side

function rhs = op_f_v_mp_mod (space, msh, f, patch_list)
   if (nargin < 4)
      patch_list = 1:msh.npatch;
   end

   if (space.npatch ~= msh.npatch)
      error ('op_f_v_mp_mod: the number of patches does not coincide');
   end

   rhs = zeros(space.ndof,1);
   for iptc=patch_list
      rhs_loc = op_f_v_tp_mod (space.sp_patch{iptc}, msh.msh_patch{iptc}, f, iptc);

      if (~isempty (space.dofs_ornt))
         rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
      end

      rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
   end
end
