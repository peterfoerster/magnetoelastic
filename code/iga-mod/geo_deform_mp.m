% INPUT:
%
%     u:         vector of dof weights for the displacement
%     space:     structure representing the space of discrete functions
%     geometry:  geometry structure. It must contain a nurbs substructure
%                 as given in the NURBS toolbox.
%
% OUTPUT:
%
%     new_geometry: geometry structure for the deformed domain.

function new_geom = geo_deform_mp (u, space, geometry);

  for iptc=1:numel(geometry)

  if (~isfield (geometry, 'nurbs'))
    error ('geo_deform: only NURBS-based geometries are allowed')
  end

  if (space.ncomp == 2)
    nurbs = geometry(iptc).nurbs;
    ndof_comp = prod (space.sp_patch{iptc}.ndof_dir, 2);

    u1 = reshape (u(1:ndof_comp(1)), [1, space.sp_patch{iptc}.ndof_dir(1,:)]);
    u2 = reshape (u(ndof_comp(1)+[1:ndof_comp(2)]), [1, space.sp_patch{iptc}.ndof_dir(2,:)]);
    weights = nurbs.coefs(4,:,:);

    nurbs.coefs(1,:,:) = nurbs.coefs(1,:,:) + u1 .* weights;
    nurbs.coefs(2,:,:) = nurbs.coefs(2,:,:) + u2 .* weights;

    new_geom(iptc).nurbs = nurbs;

    new_geom(iptc).map      = @(PTS) geo_2d_nurbs (new_geom(iptc).nurbs, PTS, 0);
    new_geom(iptc).map_der  = @(PTS) geo_2d_nurbs (new_geom(iptc).nurbs, PTS, 1);
    new_geom(iptc).map_der2 = @(PTS) geo_2d_nurbs (new_geom(iptc).nurbs, PTS, 2);

  elseif (space.ncomp == 3)
    nurbs = geometry.nurbs;
    ndof_comp = prod (space.ndof_dir, 2);

    u1 = reshape (u(1:ndof_comp(1)), [1, space.ndof_dir(1,:)]);
    u2 = reshape (u(ndof_comp(1)+[1:ndof_comp(2)]), [1, space.ndof_dir(2,:)]);
    u3 = reshape (u(ndof_comp(1)+ndof_comp(2)+[1:ndof_comp(3)]), [1, space.ndof_dir(3,:)]);
    weights = nurbs.coefs(4,:,:,:);

    nurbs.coefs(1,:,:,:) = nurbs.coefs(1,:,:,:) + u1 .* weights;
    nurbs.coefs(2,:,:,:) = nurbs.coefs(2,:,:,:) + u2 .* weights;
    nurbs.coefs(3,:,:,:) = nurbs.coefs(3,:,:,:) + u3 .* weights;

    new_geom.nurbs = nurbs;

    new_geom.map     = @(PTS) geo_3d_nurbs (new_geom.nurbs, PTS, 0);
    new_geom.map_der = @(PTS) geo_3d_nurbs (new_geom.nurbs, PTS, 1);

  end
end%for

end
