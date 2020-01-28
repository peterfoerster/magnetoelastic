% INPUT:
%
%   spu:      object representing the space of trial functions (see sp_multipatch)
%   spv:      object representing the space of test functions (see sp_multipatch)
%   msh:      object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   E, nu, G: function handles to compute the coefficients
%   patches:  list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix

function mat = op_le2d_mp (spu, spv, msh, E, nu, G, patch_list)

  if (nargin < 7)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spv.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_su_ev_mp: the number of patches does not coincide')
  end

  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_le2d_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, E, nu, G, iptc);
    rows(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  mat = sparse (rows, cols, vals, spv.ndof, spu.ndof);
  clear rows cols vals rs cs vs

end
