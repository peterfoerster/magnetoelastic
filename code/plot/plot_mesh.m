function [] = plot_mesh (geometry, nsub)
   for iptc=1:length(geometry)
      nurbs = geometry(iptc).nurbs;
      [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, nurbs.order-2);
      nurbs    = geo_load (nrbkntins (nurbs, nknots));
      fields   = fieldnames (nurbs);
      for ifld = 1:numel (fields)
         geometry(iptc).(fields{ifld}) = nurbs.(fields{ifld});
      end
   end

   for iptc=1:length(geometry)
      hold on;
      nrbkntplot(geometry(iptc).nurbs);
      hold off;
   end
end
