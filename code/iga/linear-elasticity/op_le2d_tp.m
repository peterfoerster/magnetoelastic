% INPUT:
%
%   spu:      object representing the space of trial functions (see sp_vector)
%   spv:      object representing the space of test functions (see sp_vector)
%   msh:      object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   E, nu, G: function handles to compute the coefficients
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_le2d_tp (spu, spv, msh, E, nu, G, iptc)

  mat = spalloc (spv.ndof, spu.ndof, 5*spu.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);

    spu_col = sp_evaluate_col (spu, msh_col, 'value', false, 'gradient', true);
    spv_col = sp_evaluate_col (spv, msh_col, 'value', false, 'gradient', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    for iE=1:numel(E)
      E_iel{iE} = E{iE} (x{:}, iptc);
    end
    for inu=1:numel(nu)
      nu_iel{inu} = nu{inu} (x{:}, iptc);
    end
    for iG=1:numel(G)
      G_iel{iG} = G{iG} (x{:}, iptc);
    end

    mat = mat + op_le2d (spu_col, spv_col, msh_col, E_iel, nu_iel, G_iel);
  end

  if (nargout == 1)
    varargout{1} = mat;
  elseif (nargout == 3)
    [rows, cols, vals] = find (mat);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
