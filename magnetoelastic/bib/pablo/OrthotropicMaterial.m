classdef OrthotropicMaterial
    properties
        E1;
        E2;
        E3;
        mu12;
        mu13;
        mu23;
        nu13;
        nu12;
        nu23;
        %         compute_orientation_basis; % To be set by a user as @(F, DF) ...;
        angle_xy = 0.0;
    end
    methods
        
        function D = compute_D_base(obj, n_pts)
            % Stiffness matrix in [x,y,z] coordinate system.
            D = zeros(6, 6);
            nu31 = obj.nu13 * obj.E3 / obj.E1;
            nu21 = obj.nu12 * obj.E2 / obj.E1;
            nu32 = obj.nu23 * obj.E3 / obj.E2;
            D(1:3, 1:3) = inv(...
                [1 / obj.E1,     -nu21 / obj.E2, -nu31 / obj.E3;
                -obj.nu12 / obj.E1,         1 / obj.E2, -nu32 / obj.E3;
                -obj.nu13 / obj.E1, -obj.nu23 / obj.E2,     1 / obj.E3]);
            
            D(4:6, 4:6) = diag([obj.mu12, obj.mu13, obj.mu23]);
            D = repmat(D, [1, 1, n_pts]);
            return;
        end
        
        function D = compute_D(obj, DF)
            D = full_to_voigt(obj.compute_C(DF));
            return;
        end
        
        function C = compute_C(obj, DF)
            
            %             if isa(obj.compute_orientation_basis, 'function_handle')
            %                 e = obj.compute_orientation_basis(F, DF);
            %             else
            %                 e = obj.compute_in_plane_orthonormal_basis(DF);
            %             end
            
            e = compute_change_basis_matrices(DF);
            e_rot = zeros(size(e));
            e_rot(:, 1, :) = cos(obj.angle_xy) * e(:, 1, :) + sin(obj.angle_xy) * e(:, 2, :);
            e_rot(:, 2, :) = -sin(obj.angle_xy) * e(:, 1, :) + cos(obj.angle_xy) * e(:, 2, :);
            e_rot(:, 3, :) = e(:, 3, :);
            
            s = size(e_rot);
            
            ss = 1 : numel(s);
            ss(2) = 1;
            ss(1) = 2;
            e_rot = permute(e_rot, ss);
            
            n_pts = prod(s(3:end));
            D = obj.compute_D_base(n_pts);
            C = rotate_tensor(voigt_to_full(D), e_rot);
            C = reshape(C, [3, 3, 3, 3, s(3:end)]);
            
            return;
            
        end
        
        
        % This rotation is buggy ...
        %         function R = compute_material_orientation(obj, e)
        %
        %             e_rot = zeros(size(e));
        %             e_rot(:, 1, :) = cos(obj.angle_xy) * e(:, 1, :) + sin(obj.angle_xy) * e(:, 2, :);
        %             e_rot(:, 2, :) = -sin(obj.angle_xy) * e(:, 1, :) + cos(obj.angle_xy) * e(:, 2, :);
        %             e_rot(:, 3, :) = e(:, 3, :);
        %
        %             error('The code below may be buggy. It does not work. Do not use it.');
        %
        %             lx = e_rot(1, 1, :);
        %             ly = e_rot(2, 1, :);
        %             lz = e_rot(3, 1, :);
        %             rx = e_rot(1, 2, :);
        %             ry = e_rot(2, 2, :);
        %             rz = e_rot(3, 2, :);
        %             tx = e_rot(1, 3, :);
        %             ty = e_rot(2, 3, :);
        %             tz = e_rot(3, 3, :);
        %
        %             s = size(e_rot);
        %             R = zeros([6, 6, s(3:end)]);
        %
        %             R(1, 1, :) = lx .* lx;
        %             R(1, 2, :) = rx .* rx;
        %             R(1, 3, :) = tx .* tx;
        %             R(1, 4, :) = 2.0 * lx .* ly;
        %             R(1, 5, :) = 2.0 * rx .* ry;
        %             R(1, 6, :) = 2.0 * tx .* ty;
        %             R(2, 1, :) = ly .* ly;
        %             R(2, 2, :) = ry .* ry;
        %             R(2, 3, :) = ty .* ty;
        %             R(2, 4, :) = 2.0 * ly .* lz;
        %             R(2, 5, :) = 2.0 * ry .* rz;
        %             R(2, 6, :) = 2.0 * ty .* tz;
        %             R(3, 1, :) = lz .* lz;
        %             R(3, 2, :) = rz .* rz;
        %             R(3, 3, :) = tz .* tz;
        %             R(3, 4, :) = 2.0 * lx .* lz;
        %             R(3, 5, :) = 2.0 * rx .* rz;
        %             R(3, 6, :) = 2.0 * tx .* tz;
        %             R(4, 1, :) = lx .* rx;
        %             R(4, 2, :) = lx .* tx;
        %             R(4, 3, :) = rx .* tx;
        %             R(4, 4, :) = lx .* ry + ly .* rx;
        %             R(4, 5, :) = tx .* ly + ty .* lx;
        %             R(4, 6, :) = rx .* ty + ry .* tx;
        %             R(5, 1, :) = ly .* ry;
        %             R(5, 2, :) = ly .* ty;
        %             R(5, 3, :) = ry .* ty;
        %             R(5, 4, :) = lz .* rx + lx .* rz;
        %             R(5, 5, :) = tz .* lx + tx .* lz;
        %             R(5, 6, :) = rz .* tx + rx .* tz;
        %             R(6, 1, :) = lz .* rz;
        %             R(6, 2, :) = lz .* tz;
        %             R(6, 3, :) = rz .* tz;
        %             R(6, 4, :) = ly .* rz + lz .* ry;
        %             R(6, 5, :) = ty .* lz + tz .* ly;
        %             R(6, 6, :) = ry .* tz + rz .* ty;
        %
        %         end
        
        %         function e = compute_in_plane_orthonormal_basis(obj, DF)
        %
        %             s = size(DF);
        %             if numel(s) < 3
        %                 nel = 1;
        %                 nqn = 1;
        %             elseif numel(s) < 4
        %                 nel = 1;
        %                 nqn = s(3);
        %             else
        %                 nel = s(4);
        %                 nqn = s(3);
        %             end
        %
        %
        %             g1 = reshape(DF(:, 1, :, :), [3, nqn, nel]);
        %             g2 = reshape(DF(:, 2, :, :), [3, nqn, nel]);
        %             norm_g1 = repmat(sqrt(sum(g1 .* g1, 1)), [3, 1, 1]);
        %             e1 = reshape(g1 ./ norm_g1, [1, 3, nqn, nel]);
        %             g1_x_g2 = zeros(3, nqn, nel);
        %             g1_x_g2(1, :, :) = g1(2, :, :) .* g2(3, :, :) - g1(3, :, :) .* g2(2, :, :);
        %             g1_x_g2(2, :, :) = g1(3, :, :) .* g2(1, :, :) - g1(1, :, :) .* g2(3, :, :);
        %             g1_x_g2(3, :, :) = g1(1, :, :) .* g2(2, :, :) - g1(2, :, :) .* g2(1, :, :);
        %             norm_g1_x_g2 = repmat(sqrt(sum(g1_x_g2 .* g1_x_g2, 1)), [3, 1, 1]);
        %             e3 = reshape(g1_x_g2 ./ norm_g1_x_g2, [1, 3, nqn, nel]);
        %             e2 = zeros(size(e1));
        %             e2(1, 1, :, :) = e3(1, 2, :, :) .* e1(1, 3, :, :) - e3(1, 3, :, :) .* e1(1, 2, :, :);
        %             e2(1, 2, :, :) = e3(1, 3, :, :) .* e1(1, 1, :, :) - e3(1, 1, :, :) .* e1(1, 3, :, :);
        %             e2(1, 3, :, :) = e3(1, 1, :, :) .* e1(1, 2, :, :) - e3(1, 2, :, :) .* e1(1, 1, :, :);
        %
        %
        %             e = zeros(3, 3, nqn, nel);
        %
        %             e(1, :, :, :) = e1;
        %             e(2, :, :, :) = e2;
        %             e(3, :, :, :) = e3;
        %
        %             e = permute(e, [2, 1, 3, 4]);
        %
        %         end
    end
end

function Crot = rotate_tensor(C, R)
sR = size(R);
sC = size(C);
Crot = zeros(size(C));
if prod(sR(3:end)) ~= prod(sC(5:end))
    error('Invalid dimensions');
end

n = prod(sR(3:end));
for q = 1 : n
    R_ = R(:, :, q);
    C_ = C(:, :, :, :, q);
    for i = 1 : 3
        for j = 1 : 3
            for k = 1 : 3
                for l = 1 : 3
                    v = 0.0;
                    for a = 1 : 3
                        for b = 1 : 3
                            for c = 1 : 3
                                for d = 1 : 3
                                    v = v + C_(a, b, c, d) * R_(a, i) * R_(b, j) * R_(c, k) * R_(d, l);
                                end
                            end
                        end
                    end
                    Crot(i, j, k, l, q) = v;
                end
            end
        end
    end
end


return;

end

% function Drot = rotate_voigt(D, R)
%
% Drot = full_to_voigt(rotate_tensor(voigt_to_full(D), R));
%
% return;
%
% end


function i = to_Voigt_index(k, l)
if k == l
    i = k;
else
    i = k + l + 1;
end

end


function C = voigt_to_full(D)
s = size(D);
if numel(s) > 2
    C = zeros([3, 3, 3, 3, s(3:end)]);
else
    C = zeros([3, 3, 3, 3]);
end
for i = 1 : 3
    for j = 1 : 3
        for k = 1 : 3
            for l = 1 : 3
                C(i, j, k, l, :) = D(to_Voigt_index(i, j), to_Voigt_index(k, l), :);
            end
        end
    end
end
return;
end

function D = full_to_voigt(C)
s = size(C);
if numel(s) > 4
    D = zeros([6, 6, s(5:end)]);
else
    D = zeros([6, 6]);
end

D(1, 1, :) = C(1, 1, 1, 1, :);
D(1, 2, :) = C(1, 1, 2, 2, :);
D(1, 3, :) = C(1, 1, 3, 3, :);
D(1, 4, :) = C(1, 1, 1, 2, :);
D(1, 5, :) = C(1, 1, 1, 3, :);
D(1, 6, :) = C(1, 1, 2, 3, :);
D(2, 1, :) = C(2, 2, 1, 1, :);
D(2, 2, :) = C(2, 2, 2, 2, :);
D(2, 3, :) = C(2, 2, 3, 3, :);
D(2, 4, :) = C(2, 2, 1, 2, :);
D(2, 5, :) = C(2, 2, 1, 3, :);
D(2, 6, :) = C(2, 2, 2, 3, :);
D(3, 1, :) = C(3, 3, 1, 1, :);
D(3, 2, :) = C(3, 3, 2, 2, :);
D(3, 3, :) = C(3, 3, 3, 3, :);
D(3, 4, :) = C(3, 3, 1, 2, :);
D(3, 5, :) = C(3, 3, 1, 3, :);
D(3, 6, :) = C(3, 3, 2, 3, :);
D(4, 1, :) = C(1, 2, 1, 1, :);
D(4, 2, :) = C(1, 2, 2, 2, :);
D(4, 3, :) = C(1, 2, 3, 3, :);
D(4, 4, :) = C(1, 2, 1, 2, :);
D(4, 5, :) = C(1, 2, 1, 3, :);
D(4, 6, :) = C(1, 2, 2, 3, :);
D(5, 1, :) = C(1, 3, 1, 1, :);
D(5, 2, :) = C(1, 3, 2, 2, :);
D(5, 3, :) = C(1, 3, 3, 3, :);
D(5, 4, :) = C(1, 3, 1, 2, :);
D(5, 5, :) = C(1, 3, 1, 3, :);
D(5, 6, :) = C(1, 3, 2, 3, :);
D(6, 1, :) = C(2, 3, 1, 1, :);
D(6, 2, :) = C(2, 3, 2, 2, :);
D(6, 3, :) = C(2, 3, 3, 3, :);
D(6, 4, :) = C(2, 3, 1, 2, :);
D(6, 5, :) = C(2, 3, 1, 3, :);
D(6, 6, :) = C(2, 3, 2, 3, :);

end

% Creado:     Pablo Antolin
% Modificado: 2011 Dec 12 - 12:50
% -*- coding: utf-8 -*-

function C = compute_change_basis_matrices(DF)


s = size(DF);
if numel(s) < 3
    nel = 1;
    nqn = 1;
    if (size(DF) ~= [3, 3])
        error('Invalid dimensions');
    end
elseif numel(s) < 4
    nel = 1;
    nqn = s(3);
    if (size(DF) ~= [3, 3, nqn])
        error('Invalid dimensions');
    end
else
    nel = s(4);
    nqn = s(3);
    if (size(DF) ~= [3, 3, nqn, nel])
        error('Invalid dimensions');
    end
end

npts = nqn * nel;


g1 = reshape(DF(:, 1, :), [3, npts]);
g2 = reshape(DF(:, 2, :), [3, npts]);
norm_g1 = repmat(sqrt(sum(g1 .* g1, 1)), [3, 1]);
e1 = reshape(g1 ./ norm_g1, [1, 3, npts]);
g1_x_g2 = zeros(3, nqn, nel);
g1_x_g2(1, :) = g1(2, :) .* g2(3, :) - g1(3, :) .* g2(2, :);
g1_x_g2(2, :) = g1(3, :) .* g2(1, :) - g1(1, :) .* g2(3, :);
g1_x_g2(3, :) = g1(1, :) .* g2(2, :) - g1(2, :) .* g2(1, :);
norm_g1_x_g2 = repmat(sqrt(sum(g1_x_g2 .* g1_x_g2, 1)), [3, 1]);
e3 = reshape(g1_x_g2 ./ norm_g1_x_g2, [1, 3, npts]);
e2 = zeros(size(e1));
e2(1, 1, :) = e3(1, 2, :) .* e1(1, 3, :) - e3(1, 3, :) .* e1(1, 2, :);
e2(1, 2, :) = e3(1, 3, :) .* e1(1, 1, :) - e3(1, 1, :) .* e1(1, 3, :);
e2(1, 3, :) = e3(1, 1, :) .* e1(1, 2, :) - e3(1, 2, :) .* e1(1, 1, :);


C = zeros(3, 3, npts);

C(1, :, :) = e1;
C(2, :, :) = e2;
C(3, :, :) = e3;

C = permute(C, [2, 1, 3]);

if numel(s) < 3
  C = reshape(C, [3, 3]);
elseif numel(s) < 4
  C = reshape(C, [3, 3, nqn]);
else
  C = reshape(C, [3, 3, nqn, nel]);
end

end

