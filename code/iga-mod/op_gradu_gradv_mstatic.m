% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: cell array of coefficients
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_gradu_gradv_mstatic (spu, spv, msh, coeff)
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
		   msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  % weights in physical space
  jacdet_weights = msh.jacdet .* msh.quad_weights;
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      % gradu_iel (ndim,nqn,1,nsh)
      gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel));
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
      nu11 = reshape (coeff{1}(:,iel), [1,msh.nqn,1,1]);
      nu22 = reshape (coeff{2}(:,iel), [1,msh.nqn,1,1]);

      % nu11*u_x2*v_x2 + nu22*u_x1*v_x1
      jacdet_gradu11 = bsxfun (@times, jacdet_iel.*nu22, gradu_iel(1,:,:,:));
      jacdet_gradu22 = bsxfun (@times, jacdet_iel.*nu11, gradu_iel(2,:,:,:));

      tmp11 = bsxfun (@times, jacdet_gradu11, gradv_iel(1,:,:,:));
      tmp22 = bsxfun (@times, jacdet_gradu22, gradv_iel(2,:,:,:));

      % sum over the quadrature nodes and reshape to match the shape function blocks in A
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp11+tmp22, 2), spv.nsh(iel), spu.nsh(iel));

      % compute the corresponding indices of A (shape function blocks)
      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_gradu_gradv: wrong number of output arguments')
  end

end
