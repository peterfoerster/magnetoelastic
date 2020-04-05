% INPUT:
%
%   spw: structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv: structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh: structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
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

function varargout = op_wmec1d_w1 (spw, spv, msh, b, rho, E, A, I)
   % (msh_col.nqn x nsh_max x msh_col.nel)
   valv = reshape(spv.shape_functions, 1, msh.nqn, spv.nsh_max, msh.nel);
   % (ncomp x msh_col.nqn x nsh_max x msh_col.nel)
   derw = reshape(spw.shape_function_gradients, 1, msh.nqn, spw.nsh_max, msh.nel);
   derv = reshape(spv.shape_function_gradients, 1, msh.nqn, spv.nsh_max, msh.nel);
   % (rdim x rdim x msh_col.nqn x nsh_max x msh_col.nel)
   der2w = reshape(spw.shape_function_hessians(1,1,:,:), 1, msh.nqn, spw.nsh_max, msh.nel);

   rows   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   cols   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   values = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);

   jacdet_weights = msh.jacdet .* msh.quad_weights;

   ncounter = 0;
   for iel=1:msh.nel
      if (all(msh.jacdet(:,iel)))
         valv_iel = reshape(valv(:,:,1:spv.nsh(iel),iel), 1, msh.nqn, spv.nsh(iel), 1);
         derw_iel = reshape(derw(:,:,1:spw.nsh(iel),iel), 1, msh.nqn, 1, spw.nsh(iel));
         derv_iel = reshape(derv(:,:,1:spv.nsh(iel),iel), 1, msh.nqn, spv.nsh(iel), 1);
         der2w_iel = reshape(der2w(:,:,1:spw.nsh(iel),iel), 1, msh.nqn, 1, spw.nsh(iel));

         jacdet_iel = reshape(jacdet_weights(:,iel), [1,msh.nqn,1,1]);

         tmp1 = bsxfun (@times, valv_iel, der2w_iel);
         tmp1 = bsxfun (@times, b^2*jacdet_iel*rho, tmp1);
         tmp2 = bsxfun (@times, derw_iel, derv_iel);
         tmp2 = bsxfun (@times, 1/2*E*A*jacdet_iel, tmp2);

         values(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = reshape(sum(tmp1 + tmp2, 2), spv.nsh(iel), spw.nsh(iel));

         [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spw.connectivity(:,iel));
         rows(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = rows_loc;
         cols(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = cols_loc;
         ncounter = ncounter + spw.nsh(iel)*spv.nsh(iel);

      else
         warning ('geopdes:jacdet_zero_at_quad_node', 'op_wmec1d_w1: singular map in element number %d', iel)
      end
   end

   if (nargout == 1 || nargout == 0)
      varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), values(1:ncounter), spv.ndof, spw.ndof);
   elseif (nargout == 3)
      varargout{1} = rows(1:ncounter);
      varargout{2} = cols(1:ncounter);
      varargout{3} = values(1:ncounter);
   else
      error ('op_wmec1d_w1: wrong number of output arguments')
   end
end
