function [errh1, errl2] = sp_errh1 (geometry, msh_ref, space_ref, u_ref, space, u)
   errh1 = 0; errl2 = 0;
% keyboard
   for iel=1:msh_ref.nel_dir(1)
      msh_ref_col = msh_evaluate_col (msh_ref, iel);
     [errh1_el, errl2_el] = sp_errh1_el (geometry, msh_ref_col, space_ref, u_ref, space, u);

     errh1 = errh1 + errh1_el;
     errl2 = errl2 + errl2_el;
   end
end
