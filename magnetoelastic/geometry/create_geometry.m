pkg load nurbs; clf;
filename = 'magnetoelastic.txt';

[plate] = create_plate();

hold on;
for icrv=1:length(plate)
   nrbplot(plate(icrv), 10);
   x = nrbeval(plate(icrv), 0.5);
   text(x(1), x(2), num2str(icrv));
end
hold off;

[outer_boundary] = create_outerboundary(plate);

hold on;
for icrv=1:length(outer_boundary)
   nrbplot(outer_boundary(icrv), 10);
   x = nrbeval(outer_boundary(icrv), 0.5);
   text(x(1), x(2), num2str(icrv));
end
hold off;

[domain] = discretize_domain (plate, outer_boundary);

hold on;
for icrv=1:length(domain)
   nrbplot(domain(icrv), 10);
   x = nrbeval(domain(icrv), 0.5);
   text(x(1), x(2), num2str(icrv));
end
hold off;

ptcs = create_ptcs (plate, outer_boundary, domain);

hold on;
for icrv=1:length(ptcs)
   nrbkntplot(ptcs(icrv));
   x = nrbeval(ptcs(icrv), [0.5 0.5]);
   text(x(1), x(2), num2str(icrv));
end
hold off;

write_geometryfile (ptcs, filename);
