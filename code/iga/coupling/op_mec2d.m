% INPUT:
%
%   spv: structure representing the space of mechanical test functions (see sp_vector/sp_evaluate_col)
%   spA: structure representing the space of magnetic trial functions (see sp_scalar/sp_evaluate_col)
%   msh: structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   f:   coefficients evaluated at the quadrature points
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_mec2d (spv, spA, msh, f)
  % gradv (ncomp x rdim x msh_col.nqn x nsh_max x msh_col.nel)
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
  gradA = reshape (spA.shape_function_gradients, spA.ncomp, [], msh.nqn, spA.nsh_max, msh.nel);

  ndir = size (gradv, 2);

  rows = zeros (msh.nel * spv.nsh_max * spA.nsh_max, 1);
  cols = zeros (msh.nel * spv.nsh_max * spA.nsh_max, 1);
  values = zeros (msh.nel * spv.nsh_max * spA.nsh_max, 1);

  % weights and coefficients in physical space
  jacdet_weights = msh.jacdet .* msh.quad_weights;
  f11  = f{1};
  f21  = f{2};
  f62  = f{3};

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp, ndir, msh.nqn, spv.nsh(iel), 1);
      gradA_iel = reshape (gradA(:,:,:,1:spA.nsh(iel),iel), spA.ncomp, ndir, msh.nqn, 1, spA.nsh(iel));

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,1,msh.nqn,1,1]);
      f11_iel    = reshape (f11(:,iel), [1,1,msh.nqn,1,1]);
      f21_iel    = reshape (f21(:,iel), [1,1,msh.nqn,1,1]);
      f62_iel    = reshape (f62(:,iel), [1,1,msh.nqn,1,1]);

      tmp11 = bsxfun(@times, f11_iel, gradv_iel(1,1,:,:,:));
      tmp12 = bsxfun(@times, f21_iel, gradv_iel(2,2,:,:,:));
      tmp1  = bsxfun(@times, gradA_iel(1,2,:,:,:), tmp11 + tmp12);

      tmp2 = gradv_iel(1,2,:,:,:) + gradv_iel(2,1,:,:,:);
      tmp2 = bsxfun(@times, f62_iel, tmp2);
      tmp2 = bsxfun(@times, gradA_iel(1,1,:,:,:), tmp2);

      tmp = 1/2 * bsxfun(@times, jacdet_iel, tmp1 - tmp2);

      values(ncounter+(1:spv.nsh(iel)*spA.nsh(iel))) = reshape (sum (tmp, 3), spv.nsh(iel), spA.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spA.connectivity(:,iel));
      rows(ncounter+(1:spv.nsh(iel)*spA.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spv.nsh(iel)*spA.nsh(iel))) = cols_loc;
      ncounter = ncounter + spv.nsh(iel)*spA.nsh(iel);

      % catch NaN values originating from non-involved patches
      idx = isnan(values);
      values(idx) = 0;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_mec2d: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spA.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_mec2d: wrong number of output arguments')
  end

end
