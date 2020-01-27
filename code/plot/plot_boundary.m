function [] = plot_boundary (geometry, boundaries)
   width = 4;
   nsub  = 8;
   for ib=1:length(boundaries)
      switch(ib)
         % homogeneous Dirichlet
         case{1}
            for ip=1:length(boundaries(ib).patches)
               nrbplot_color('b', width, geometry(boundaries(ib).patches(ip)).boundary(boundaries(ib).faces(ip)).nurbs, nsub);
            end
         % inhomogeneous Neumann
         case{2}
            for ip=1:length(boundaries(ib).patches)
               nrbplot_color('r', width, geometry(boundaries(ib).patches(ip)).boundary(boundaries(ib).faces(ip)).nurbs, nsub);
            end
         % homogeneous Neumann
         case{3}
            for ip=1:length(boundaries(ib).patches)
               nrbplot_color('g', width, geometry(boundaries(ib).patches(ip)).boundary(boundaries(ib).faces(ip)).nurbs, nsub);
            end
      end
   end%for
end
