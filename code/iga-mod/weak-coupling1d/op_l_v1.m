% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   D:     coupling coefficient?
%   B:     magnetic flux density
%   A:     cross-sectional area
%
% OUTPUT:
%
%   rhs: assembled right-hand side

function rhs = op_l_v1 (spv, msh, D, B, A)

   derv = reshape(spv.shape_function_gradients, 1, msh.nqn, spv.nsh_max, msh.nel);
   rhs  = zeros(spv.ndof,1);

   jacdet_weights = msh.jacdet .* msh.quad_weights;

   for iel=1:msh.nel
      if (all(msh.jacdet(:,iel)))
         derv_iel   = reshape(derv(:,:,1:spv.nsh(iel),iel), 1, msh.nqn, spv.nsh(iel));
         jacdet_iel = reshape(jacdet_weights(:,iel), [1,msh.nqn,1]);

         tmp = bsxfun (@times, D*B*A*jacdet_iel, derv_iel);

         rhs_loc = sum(tmp, 2);
         rhs(spv.connectivity(1:spv.nsh(iel), iel)) = rhs(spv.connectivity(1:spv.nsh(iel), iel)) + rhs_loc(:);
      else
         warning ('geopdes:jacdet_zero_at_quad_node', 'op_l_v1: singular map in element number %d', iel)
      end
   end
end
