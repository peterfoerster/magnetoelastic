% INPUT:
%
%   spA:  structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spAt: structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:  structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   mu:   cell array of permeabilities
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_gradu_gradv_ms (spA, spAt, msh, mu)
  gradA  = reshape (spA.shape_function_gradients, spA.ncomp, [], ...
		   msh.nqn, spA.nsh_max, msh.nel);
  gradAt = reshape (spAt.shape_function_gradients, spAt.ncomp, [], ...
		   msh.nqn, spAt.nsh_max, msh.nel);

  ndir = size (gradA, 2);

  rows   = zeros (msh.nel * spA.nsh_max * spAt.nsh_max, 1);
  cols   = zeros (msh.nel * spA.nsh_max * spAt.nsh_max, 1);
  values = zeros (msh.nel * spA.nsh_max * spAt.nsh_max, 1);

  % weights in physical space
  jacdet_weights = msh.jacdet .* msh.quad_weights;

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      % gradA_iel (ndim x nqn x 1 x nsh)
      gradA_iel  = reshape (gradA(:,:,:,1:spA.nsh(iel),iel), spA.ncomp*ndir, msh.nqn, 1, spA.nsh(iel));
      gradAt_iel = reshape (gradAt(:,:,:,1:spAt.nsh(iel),iel), spAt.ncomp*ndir, msh.nqn, spAt.nsh(iel), 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
      mu11_iel   = reshape (mu{1}(:,iel), [1,msh.nqn,1,1]);
      mu22_iel   = reshape (mu{2}(:,iel), [1,msh.nqn,1,1]);

      tmp1 = bsxfun (@times, jacdet_iel./mu22_iel, gradA_iel(1,:,:,:));
      tmp1 = bsxfun (@times, tmp1, gradAt_iel(1,:,:,:));

      tmp2 = bsxfun (@times, jacdet_iel./mu11_iel, gradA_iel(2,:,:,:));
      tmp2 = bsxfun (@times, tmp2, gradAt_iel(2,:,:,:));

      values(ncounter+(1:spA.nsh(iel)*spAt.nsh(iel))) = reshape (sum (tmp1 + tmp2, 2), spAt.nsh(iel), spA.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spAt.connectivity(:,iel), spA.connectivity(:,iel));
      rows(ncounter+(1:spA.nsh(iel)*spAt.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spA.nsh(iel)*spAt.nsh(iel))) = cols_loc;
      ncounter = ncounter + spA.nsh(iel)*spAt.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv_ms: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spAt.ndof, spA.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_gradu_gradv_ms: wrong number of output arguments')
  end

end
