% INPUT:
%
%   spu:      structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%   spv:      structure representing the space of test functions (see sp_vector/sp_evaluate_col)
%   msh:      structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   E, nu, G: coefficients evaluated at the quadrature points
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_linear_elasticity2d (spu, spv, msh, E, nu, G)
  % gradu (ncomp x rdim x msh_col.nqn x nsh_max x msh_col.nel)
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  % weights in physical space and coefficients
  jacdet_weights = msh.jacdet .* msh.quad_weights;
  E1   = E{1};
  E2   = E{2};
  nu12 = nu{1};
  G12  = G{1};
  jacdet_weights = 1 ./ (E1 - (nu12.^2) .* E2) .* jacdet_weights;

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp, ndir, msh.nqn, 1, spu.nsh(iel));
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp, ndir, msh.nqn, spv.nsh(iel), 1);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,1,msh.nqn,1,1]);
      E1_iel     = reshape (E1(:,iel), [1,1,msh.nqn,1,1]);
      E2_iel     = reshape (E2(:,iel), [1,1,msh.nqn,1,1]);
      nu12_iel   = reshape (nu12(:,iel), [1,1,msh.nqn,1,1]);
      G12_iel    = reshape (G12(:,iel), [1,1,msh.nqn,1,1]);

      tmp1 = bsxfun(@times, (E1_iel.^2), gradu_iel(1,1,:,:,:));
      tmp1 = bsxfun(@times, tmp1, gradv_iel(1,1,:,:,:));

      tmp21 = bsxfun(@times, gradu_iel(1,1,:,:,:), gradv_iel(2,2,:,:,:));
      tmp22 = bsxfun(@times, gradu_iel(2,2,:,:,:), gradv_iel(1,1,:,:,:));
      tmp2  = bsxfun(@times, nu12_iel .* E2_iel .* E1_iel, tmp21 + tmp22);

      tmp3 = bsxfun(@times, E1_iel .* E2_iel, gradu_iel(2,2,:,:,:));
      tmp3 = bsxfun(@times, tmp3, gradv_iel(2,2,:,:,:));

      tmp41 = gradu_iel(1,2,:,:,:) + gradu_iel(2,1,:,:,:);
      tmp42 = gradv_iel(1,2,:,:,:) + gradv_iel(2,1,:,:,:);
      tmp4  = bsxfun(@times, tmp41, tmp42);
      tmp4  = bsxfun(@times, G12_iel .* (E1_iel - (nu12_iel.^2) .* E2_iel), tmp4);

      tmp = bsxfun(@times, jacdet_iel, tmp1 + tmp2 + tmp3 + tmp4);

      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp, 3), spv.nsh(iel), spu.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_linear_elasticity2d: singular map in element number %d', iel)
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
    error ('op_linear_elasticity2d: wrong number of output arguments')
  end

end
