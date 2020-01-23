function [] = plot_geometry (geometry, boundaries)
   for iptc=1:length(geometry)
      hold on;
      nrbkntplot(geometry(iptc).nurbs);
      x = geo_nurbs(geometry(iptc).nurbs, geometry(iptc).dnurbs, geometry(iptc).dnurbs2, {0.5,0.5}, 0, geometry(iptc).rdim);
      text(x(1), x(2), num2str(iptc));
      hold off;
   end

   hold on;
   plot_boundary(geometry, boundaries);
   hold off;
   view(2);
end
