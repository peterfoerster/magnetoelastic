% INPUT:
%
%  A:        vector of dof weights
%  space:    object defining the discrete space (see sp_multipatch)
%  geometry: an array of geometry structures (see mp_geo_load)
%  npts:     cell array with coordinates of points along each parametric direction

function plot_curl2d_mp (A, space, geometry, npts)
   hold_flag = ishold;
   for iptc=1:space.npatch
      if (isempty(space.dofs_ornt))
         plot_curl2d_dat (A(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc), npts, iptc);
      end
      hold on
   end

   if (~hold_flag)
      hold off
   end
end
