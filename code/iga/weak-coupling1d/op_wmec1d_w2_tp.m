% INPUT:
%
%   spw: object representing the space of trial functions (see sp_scalar)
%   spv: object representing the space of test functions (see sp_scalar)
%   msh: object defining the domain partition and the quadrature rule (see msh_cartesian)
%   b:   width
%   rho: density
%   E:   Young's modulus
%   A:   cross-sectional area
%   I:   axial moment of inertia
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_wmec1d_w2_tp (spw, spv, msh, b, rho, E, A, I)

   mat = spalloc (spv.ndof, spw.ndof, 3*spw.ndof);

   for iel=1:msh.nel_dir(1)
      msh_col = msh_evaluate_col (msh, iel);
      spw_col = sp_evaluate_col (spw, msh_col, 'value', false, 'hessian', true);
      spv_col = sp_evaluate_col (spv, msh_col, 'value', true, 'hessian', true);

      mat = mat + op_wmec1d_w2 (spw_col, spv_col, msh_col, b, rho, E, A, I);
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
