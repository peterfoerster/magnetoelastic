function [errh1, errl2] = sp_errh1 (geometry, msh_ref, space_ref, u_ref, space, u)
   errh1 = 0; errl2 = 0;

   for iel=1:msh_ref.nel_dir(1)
      msh_ref_col = msh_evaluate_col (msh_ref, iel);
      sp_ref_col  = sp_evaluate_col (space_ref, msh_col, 'value', true, 'gradient', true);
      sp_col      = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', true);

     [errh1_el, errl2_el] = sp_errh1_el (geometry, msh_ref_col, sp_ref_col, u_ref, sp_col, u);

     errh1 = errh1 + errh1_el;
     errl2 = errl2 + errl2_el;
   end
end
