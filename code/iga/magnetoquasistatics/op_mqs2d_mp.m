% INPUT:
%
%   spA:     object representing the space of trial functions (see sp_multipatch)
%   spAt:    object representing the space of test functions (see sp_multipatch)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   mu:      cell array of function handles to compute the magnetic permeability.
%   kappa:   complex conductivity.
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix

function mat = op_mqs2d_mp (spA, spAt, msh, mu, kappa, patch_list)
   if (nargin < 6)
      patch_list = 1:msh.npatch;
   end

   if ((spA.npatch ~= spAt.npatch) || (spA.npatch ~= msh.npatch))
      error('op_mqs2d_mp: the number of patches does not coincide')
   end

   ncounter = 0;
   for iptc=patch_list
      [rs, cs, vs] = op_mqs2d_tp (spA.sp_patch{iptc}, spAt.sp_patch{iptc}, msh.msh_patch{iptc}, mu, kappa, iptc);

      rows(ncounter+(1:numel(rs))) = spAt.gnum{iptc}(rs);
      cols(ncounter+(1:numel(rs))) = spA.gnum{iptc}(cs);

      if (~isempty(spAt.dofs_ornt))
         vs = spAt.dofs_ornt{iptc}(rs)' .* vs;
      end
      if (~isempty(spA.dofs_ornt))
         vs = vs .* spA.dofs_ornt{iptc}(cs)';
      end

      vals(ncounter+(1:numel(rs))) = vs;
      ncounter = ncounter + numel(rs);
   end

   mat = sparse (rows, cols, vals, spAt.ndof, spA.ndof);
   clear rows cols vals rs cs vs
end
