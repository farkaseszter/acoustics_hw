clear;

% Import the mesh and its parameters
sphere = gmsh_import_mesh('../meshes/sphere.msh');
params = gmsh_import_params('../meshes/sphere.geo');

% Create field mesh
field = mesh_create_slab([0 -2 -3.5; 0 2 -3.5; 0 2 3.5; 0 -2 3.5], [200 200]);
%field = mesh_create_slab([params.r*1.1 0; params.r*5 0; params.r*5 2*pi; params.r*1.1 2*pi], [100 350]);
%field.nodes = [ zeros(size(field.nodes, 1), 1),field.nodes(:,1).*cos(field.nodes(:,2)),field.nodes(:,1).*sin(field.nodes(:,2)) ];
field = mesh_section_circular(field, params.r,'none');

% Create plane wave source
plane_w = mesh_create_slab([3 4 -7; 3 -4 -7; -3 -4 -7; -3 4 -7], [50 50]);


% Plot the surface and field meshes
%figure;
%plot_mesh(speaker); hold on;
%plot_mesh(field);
%plot_mesh(plane_w);

%% Determine constants - article observes ka E 0.1, 2, 4 for a = 1 m
c = 343;            % Speed of sound [m/s]
rho0 = 1.21;        % Average density [kg/m^3]
k = 4;
f = k*c/(2*pi());
%% Define BC-s

% Reflection BC on sphere - Neumann BC, dp/dn = 0
% Overcomplicated code to test nonzero values
[cnt, nrm] = mesh_element_props(sphere);
sel = cnt(:,3)==cnt(:,3);
n_elements = size(sphere.elements, 1);    % Number of elements
a = zeros(n_elements, 3);
a(sel, 3) = 0;                            % Acceleration in z [m/s^2]
an_2 = dot(a, nrm, 2);                    % Acceleration in z [m/s^2]
% dp/dn on the surface
qs_2 = - rho0 * an_2;

% Plane wave BC - Dirichlet BC, p0 is supposed to be 1
[cntw, nrmw] = mesh_element_props(plane_w);
sel_w = cntw(:,3)==cntw(:,3);
n_elements_w = size(plane_w.elements, 1);     % Number of elements
kvec = [0; 0; 1]*k;
% p0 scaled up so waves in field reach 1 Pa
p_1 = 14.5*exp(1i * cntw * kvec);

%% Set up the BEM matrices for the BC-s
fprintf('BEM surface assembly ... '); tic;
[Gs11, Hs11] = bem_matrices(k, plane_w);
[Gs12, Hs12] = bem_matrices(k, sphere,cntw);
[Gs21, Hs21] = bem_matrices(k, plane_w,cnt);
[Gs22, Hs22] = bem_matrices(k, sphere);
fprintf('Ready in %.2f s\n', toc);


%% Solve the system
% Solve the surface equations
fprintf('Solving boundary equation ... '); tic;
Z = [ Hs12, -Gs11; (-0.5*eye(size(Hs22, 1))+Hs22), -Gs21];
rhs = [(0.5*eye(size(Hs11, 1))-Hs11) Gs12; -Hs21 Gs22]*[p_1;qs_2];

% Solve system
sol = Z \ rhs;
q1 = sol(n_elements+1:end);
p2 = sol(1:n_elements);
fprintf('Ready in %.2f s\n', toc);

% Solve the equations for the radiated field
fprintf('BEM field assembly ... '); tic;
[Gf_sp, Hf_sp] = bem_matrices(k, sphere, field.nodes);
[Gf_pw, Hf_pw] = bem_matrices(k, plane_w, field.nodes);
fprintf('Ready in %.2f s\n', toc);

%% Calculate pressure in field
fprintf('Computing field ... '); tic;
pf = Hf_sp * p2 - Gf_sp * qs_2 + Hf_pw * p_1 - Gf_pw*q1;
pf_sphere = Hf_sp * p2 - Gf_sp * qs_2;
pf_pw = Hf_pw * p_1 - Gf_pw*q1;
fprintf('Ready in %.2f s\n', toc);

%% Rotation before plot
rot = [1 0 0; 0 0 1; 0 -1 0];
sphere.nodes = sphere.nodes*rot;
field.nodes = field.nodes*rot;

%% Plot the solution
figure;
%plot_mesh(speaker, real(ps));
hold on
h = plot_mesh(field, real(pf));
%h = plot_mesh(field, real(pf_sphere));
%h = plot_mesh(field, real(pf_pw));
set(h, 'EdgeColor', 'none');
caxis([min(real(pf)) max(real(pf))]);
axis off;
ylim([-4 4]);
zlim([-5 5]);
view([-45 10]);
cb = colorbar('SouthOutside');
ylabel(cb, 'Pressure field Re(p) [Pa]', 'FontSize', 12);
xlabel('x [m]');
zlabel('z [m]');

%% Calculate field at r = 1.5 m
query_points = mesh_section_circle(field, 1.5*params.r,'any');
query_points = mesh_section_circular(query_points, 1.5*params.r,'nall');
figure;
g = plot_mesh(query_points,real(pf));
set(g, 'EdgeColor', 'none');
cb = colorbar('SouthOutside');
axis off;