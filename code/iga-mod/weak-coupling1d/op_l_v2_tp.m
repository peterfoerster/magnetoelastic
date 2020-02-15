% INPUT:
%
%   space: object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   D:     coupling coefficient?
%   B:     magnetic flux density
%   b:     width
%
% OUTPUT:
%
%   rhs: assembled right-hand side

function rhs = op_l_2v_tp (space, msh, D, B, b)

   rhs = zeros(space.ndof,1);

   for iel=1:msh.nel_dir(1)
      msh_col = msh_evaluate_col (msh, iel);
      sp_col  = sp_evaluate_col (space, msh_col, 'value', false, 'hessian', true);

      rhs = rhs + op_l_v2 (sp_col, msh_col, D, B, b);
   end
end
