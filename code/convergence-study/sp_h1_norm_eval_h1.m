% evaluate h1

function [l2_ref, l2, h1s_ref, h1s, errh1_elem] = sp_h1_norm_eval_h1 (sp, msh, u_ref, space_ref, u, space, geometry)
% get quadrature points
msh_pts_x = msh.qn{1};
msh_pts_y = msh.qn{2};
msh_pts_y = reshape( msh_pts_y, 1, numel(msh_pts_y) );

% evaluate solution
grad_valu_ref = sp_eval(u_ref, space_ref, geometry, {msh_pts_x, msh_pts_y}, 'gradient');
grad_valu = sp_eval(u, space, geometry, {msh_pts_x, msh_pts_y}, 'gradient');

% reshape the solution in two steps to exactly match the original
grad_valu_ref = reshape(grad_valu_ref, msh.rdim, msh.nqn, msh.nel);
grad_valu = reshape(grad_valu, msh.rdim, msh.nqn, msh.nel);

grad_valu_ref = reshape (grad_valu_ref, sp.ncomp, msh.rdim, msh.nqn, msh.nel);
grad_valu = reshape (grad_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);

% compute determinant and quadrature weights
w = msh.quad_weights .* msh.jacdet;
% evaluate l2
[l2_ref, l2, errl2_elem] = sp_l2_norm_eval_l2 (sp, msh, u_ref, space_ref, u, space, geometry);

errh1s_elem = sum (reshape (sum (sum ((grad_valu - grad_valu_ref).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
errh1_elem  = sqrt (errl2_elem.^2 + errh1s_elem);

h1s_elem = sum (reshape (sum (sum ((grad_valu).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
h1s = sqrt (sum (h1s_elem));
h1s_elem_ref = sum (reshape (sum (sum ((grad_valu_ref).^2, 1), 2), [msh.nqn, msh.nel]) .* w);
h1s_ref = sqrt (sum (h1s_elem_ref));
end
