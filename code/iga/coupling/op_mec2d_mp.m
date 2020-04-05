% INPUT:
%
%   spacev:  object representing the space of mechanical test functions (see sp_vector)
%   spaceA:  object representing the space of magnetic trial functions (see sp_scalar)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   f:       function handles to compute the coefficients
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix

function mat = op_mec2d_mp (spv, spA, msh, f, patch_list)

  if (nargin < 7)
    patch_list = 1:msh.npatch;
  end

  if ((spv.npatch ~= spA.npatch) || (spv.npatch ~= msh.npatch))
    error ('op_mec2d_mp: the number of patches does not coincide')
  end

  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_mec2d_tp (spv.sp_patch{iptc}, spA.sp_patch{iptc}, msh.msh_patch{iptc}, f, iptc);
    rows(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spA.gnum{iptc}(cs);
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  mat = sparse (rows, cols, vals, spv.ndof, spA.ndof);
  clear rows cols vals rs cs vs

end
