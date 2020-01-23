% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_scalar)
%   spv:     object representing the space of test functions (see sp_scalar)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%   coeff:   cell array of function handles for coefficients
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_gradu_gradv_tp_mstatic (space1, space2, msh, coeff, iptc)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);
    if (nargin == 5)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      for ic=1:numel(coeff)
         coeffs{ic} = coeff{ic} (x{:}, iptc);
      end
    else
      coeffs = ones (msh_col.nqn, msh_col.nel);
    end

    A = A + op_gradu_gradv_mstatic (sp1_col, sp2_col, msh_col, coeffs);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
