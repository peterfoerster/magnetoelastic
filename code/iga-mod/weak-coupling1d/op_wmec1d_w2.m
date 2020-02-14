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

function varargout = op_wmec1d (spw, spv, msh, b, rho, E, A, I)
   % (rdim x rdim x msh_col.nqn x nsh_max x msh_col.nel)
   der2w = reshape(spw.shape_function_hessians(1,1,:,:), 1, msh.nqn, spw.nsh_max, msh.nel);
   der2v = reshape(spv.shape_function_hessians(1,1,:,:), 1, msh.nqn, spv.nsh_max, msh.nel);

   rows   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   cols   = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);
   values = zeros(msh.nel*spw.nsh_max*spv.nsh_max,1);

   jacdet_weights = msh.jacdet .* msh.quad_weights;

   ncounter = 0;
   for iel=1:msh.nel
      if (all(msh.jacdet(:,iel)))
         der2w_iel = reshape(der2w(:,:,1:spw.nsh(iel),iel), 1, msh.nqn, 1, spw.nsh(iel));
         der2v_iel = reshape(der2v(:,:,1:spv.nsh(iel),iel), 1, msh.nqn, spv.nsh(iel), 1);

         jacdet_iel = reshape(jacdet_weights(:,iel), [1,msh.nqn,1,1]);

         tmp = bsxfun (@times, der2w_iel, der2v_iel);
         tmp = bsxfun (@times, 1/2*E*I*jacdet_iel, tmp);

         values(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = reshape(sum(tmp, 2), spv.nsh(iel), spw.nsh(iel));

         [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spw.connectivity(:,iel));
         rows(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = rows_loc;
         cols(ncounter+(1:spw.nsh(iel)*spv.nsh(iel))) = cols_loc;
         ncounter = ncounter + spw.nsh(iel)*spv.nsh(iel);

      else
         warning ('geopdes:jacdet_zero_at_quad_node', 'op_wmec1d_w2: singular map in element number %d', iel)
      end
   end

   if (nargout == 1 || nargout == 0)
      varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), values(1:ncounter), spv.ndof, spw.ndof);
   elseif (nargout == 3)
      varargout{1} = rows(1:ncounter);
      varargout{2} = cols(1:ncounter);
      varargout{3} = values(1:ncounter);
   else
      error ('op_wmec1d_w2: wrong number of output arguments')
   end
end
