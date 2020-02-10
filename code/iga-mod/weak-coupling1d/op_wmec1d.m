% INPUT:
%
%   spw:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
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

function varargout = op_gradgradu_gradgradv (spw, spv, msh, b, rho, E, A, I)

   derw = reshape (spw.shape_function_gradients, spw.ncomp, [], msh.nqn, spw.nsh_max, msh.nel);
   derv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

   der2w = reshape (spw.shape_function_hessians, spw.ncomp, [], msh.nqn, spw.nsh_max, msh.nel);
   der2v = reshape (spv.shape_function_hessians, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

   ndir   = size(der2w,2);
   rows   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   cols   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   values = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);

   jacdet_weights = msh.jacdet .* msh.quad_weights;
keyboard
   ncounter = 0;
   for iel=1:msh.nel
      if (all(msh.jacdet(:,iel)))
         der2u_iel = reshape (der2u(:,:,:,1:spu.nsh(iel),iel), spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel));
         der2v_iel = reshape (der2v(:,:,:,1:spv.nsh(iel),iel), spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1);

         jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

         jacdet_der2u = bsxfun (@times, jacdet_iel, der2u_iel);
         tmp1 = sum (bsxfun (@times, jacdet_der2u, der2v_iel), 1);
         values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel));

         [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
         rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
         cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
         ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

      else
         warning ('geopdes:jacdet_zero_at_quad_node', 'op_wmec1d: singular map in element number %d', iel)
      end
   end

   if (nargout == 1 || nargout == 0)
      varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), values(1:ncounter), spv.ndof, spw.ndof);
   elseif (nargout == 3)
      varargout{1} = rows(1:ncounter);
      varargout{2} = cols(1:ncounter);
      varargout{3} = values(1:ncounter);
   else
      error ('op_wmec1d: wrong number of output arguments')
   end
end
