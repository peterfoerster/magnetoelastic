%  INPUT:
%
%  A:          vector of dof weights
%  space:      object representing the space of discrete functions (see sp_multipatch)
%  geometry:   geometry structure (see geo_load)
%  npts:       number of points along each parametric direction where to evaluate
%  filename:   name of the output file.
%  fieldnames: how to name the saved variables in the vtk file
%  omega:      angular frequency
%  T:          vector of timesteps

function sp_to_vtk_mp_curl2d_fd (A, space, geometry, npts, filename, fieldname, omega, T)
   fid = fopen([filename '.pvd'], 'w');
   fprintf(fid, '<?xml version="1.0"?>\n');
   fprintf(fid, '<VTKFile type="Collection" version="0.1">\n');
   fprintf(fid, '<Collection>\n');

   i_cmp = complex(0, 1);
   for it=1:length(T)
      A_td = real(A*exp(i_cmp*omega*T(it)));

      for iptc=1:space.npatch
         filename_ptc = [filename '_' num2str(it) '_' num2str(iptc) '.vts'];
         fprintf(fid, '<DataSet timestep="%i" part="%i" file="%s"/>\n', it, iptc, filename_ptc);
         if(isempty(space.dofs_ornt))
            sp_to_vtk_curl2d (A_td(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc), npts, filename_ptc, fieldname);
         else
            sp_to_vtk (A_td(space.gnum{iptc}) .* space.dofs_ornt{iptc}', space.sp_patch{iptc}, geometry(iptc), npts, filename_ptc, fieldname);
         end
      end
   end

   fprintf(fid, '</Collection>\n');
   fprintf(fid, '</VTKFile>');
   fclose(fid);
end
