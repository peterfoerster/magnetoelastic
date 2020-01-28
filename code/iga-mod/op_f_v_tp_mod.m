% INPUT:
%
%   space: object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   f:     function handle to compute the source function
%
% OUTPUT:
%
%   rhs: assembled right-hand side

function rhs = op_f_v_tp_mod (space, msh, f, iptc)

  rhs = zeros (space.ndof, 1);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    rhs = rhs + op_f_v (sp_col, msh_col, f (x{:}, iptc));
  end

end
