% INPUT:
%
%     u:        vector of dof weights for the displacement
%     space:    structure representing the space of discrete functions
%     geometry: geometry structure. It must contain a nurbs substructure as given in the NURBS toolbox.
%
% OUTPUT:
%
%     geometry_def: geometry structure for the deformed domain.

function geometry_def = geo_deform_mp (u, space, geometry);
   for iptc=1:numel(geometry)
      geometry_def(iptc) = geo_deform (u(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc));
   end
end
