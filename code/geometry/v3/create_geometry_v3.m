pkg load nurbs; clf;
filename = 'magnetoelastic_v3';

[plate] = create_plate_v2();
[outer_boundary] = create_outerboundary_v2 (plate);
[domain] = discretize_domain_v2 (plate, outer_boundary);

ptcs = create_ptcs_v3 (plate, outer_boundary, domain);

hold on;
for icrv=1:length(ptcs)
   nrbkntplot(ptcs(icrv));
   x = nrbeval(ptcs(icrv), [0.5 0.5]);
   text(x(1), x(2), num2str(icrv));
end
hold off;

write_geometryfile_v3 (ptcs, filename);
