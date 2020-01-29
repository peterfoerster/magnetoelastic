% INPUT:
%
%   spA:  object representing the space of trial functions (see sp_scalar)
%   spAt: object representing the space of test functions (see sp_scalar)
%   msh:  object defining the domain partition and the quadrature rule (see msh_cartesian)
%   mu:   cell array of function handles for magnetic permeability
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_ms2d_tp (spA, spAt, msh, mu, iptc)

  mat = spalloc (spAt.ndof, spA.ndof, 3*spA.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col  = msh_evaluate_col (msh, iel);
    spA_col  = sp_evaluate_col (spA, msh_col, 'value', false, 'gradient', true);
    spAt_col = sp_evaluate_col (spAt, msh_col, 'value', false, 'gradient', true);

   for idim = 1:msh.rdim
     x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
   end
   for imu=1:numel(mu)
      mu_iel{imu} = mu{imu} (x{:}, iptc);
   end

    mat = mat + op_ms2d (spA_col, spAt_col, msh_col, mu_iel);
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
