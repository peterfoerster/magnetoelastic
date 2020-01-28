% INPUT:
%
%   spacev:   object representing the space of mechanical test functions (see sp_vector)
%   spaceA:   object representing the space of magnetic trial functions (see sp_vector)
%   msh:      object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   fc: function handles to compute the coefficients
%
% OUTPUT:
%
%   Ct:     assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_linear_elasticity2d_tp (spacev, spaceA, msh, fc, iptc)
  % (1,u)->v, (2,v)->A
  Ct = spalloc (spaceA.ndof, spacev.ndof, 5*spacev.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);

    spv_col = sp_evaluate_col (spacev, msh_col, 'value', false, 'gradient', true);
    spA_col = sp_evaluate_col (spaceA, msh_col, 'value', false, 'gradient', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    for ifc=1:numel(fc)
      fc_iel{ifc} = fc{ifc} (x{:}, iptc);
    end

    Ct = Ct + op_linear_elasticity2d (spv_col, spA_col, msh_col, fc_iel);
  end

  if (nargout == 1)
    varargout{1} = Ct;
  elseif (nargout == 3)
    [rows, cols, vals] = find (Ct);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
