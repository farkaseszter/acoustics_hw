clear;

% Import the mesh and its parameters
sphere = gmsh_import_mesh('../meshes/sphere.msh');
params = gmsh_import_params('../meshes/sphere.geo');

% Create field mesh
field = mesh_create_slab([-3.5 -3 0; -3.5 3 0; 3.5 3 0; 3.5 -3 0], [100 100]);

field = mesh_section_circular(field, params.r,'none');

% Plot the surface and field meshes
fig = figure;
t = tiledlayout(1,2);
nexttile(1);
plot_mesh(sphere); hold on;
plot_mesh(field);
view(0, 90)
xlabel('$x [m]$','interpreter','latex')
ylabel('$y [m]$','interpreter','latex')
nexttile(2);
plot_elem_normals(sphere);
xlabel('$x [m]$','interpreter','latex')
ylabel('$y [m]$','interpreter','latex')
t.TileSpacing = 'compact';
t.Padding = 'compact';
fig.Position = [100 100 900 400];

%% Determine constants - article observes ka E 0.1, 2, 4 for a = 1 m
c = 343;            % Speed of sound [m/s]
rho0 = 1.21;        % Average density [kg/m^3]
k = 4;
f = k*c/(2*pi());
%% Define BC-s
% Reflection BC on sphere - Neumann BC, dp/dn = 0
[cnt, nrm] = mesh_element_props(sphere);
direction = [1 0 0];
[ps_inc, qs_inc] = incident_field('plane', direction, cnt, nrm, k);
% qs_tot = qs_scat + qs_inc = 0
qs_scat = -qs_inc;

%% Set up the BEM matrices for the BC-s
fprintf('BEM surface assembly ... '); tic;
[Gs, Hs] = bem_matrices(k, sphere);
fprintf('Ready in %.2f s\n', toc);


%% Solve the system
% Solve the surface equations
fprintf('Solving boundary equation ... '); tic;
ps_scat = (Hs - 0.5*eye(size(Hs, 1))) \ (Gs * qs_scat);
fprintf('Ready in %.2f s\n', toc);

%% 
% Solve the equations for the radiated field
fprintf('BEM field assembly ... '); tic;
[Gf, Hf] = bem_matrices(k, sphere, field.nodes);
fprintf('Ready in %.2f s\n', toc);

%% 
[cntf, nrmf] = mesh_element_props(field);
[pf_inc] = incident_field('plane', direction, field.nodes, nrmf, k );

%% Calculate pressure in field
fprintf('Computing field ... '); tic;
pf_scat = Hf * ps_scat - Gf * qs_scat;
pf= pf_inc + pf_scat;
fprintf('Ready in %.2f s\n', toc);

%% Analytical
[ps_scat_a, pf_scat_a, err] = ana_sphere_scat(k, sphere, field.nodes, 1);
pf_a = pf_inc  + pf_scat_a;

%% Plot the solution
fig_res = figure;
t_res = tiledlayout(2,2);
nexttile(1);
hold on
plot_mesh(sphere, real(ps_scat));
h = plot_mesh(field, real(pf));
set(h, 'EdgeColor', 'none');
caxis([min(real(pf)) max(real(pf))]);
axis off;
view([0, 90]);
cb = colorbar('SouthOutside');
xlabel('$x [m]$','interpreter','latex')
ylabel('$Pressure field Re(p) [Pa]$','interpreter','latex')
title('Total pressure field - numerical solution', 'interpreter', 'latex');

nexttile(2);
hold on
plot_mesh(sphere, real(ps_scat_a));
h = plot_mesh(field, real(pf_a));
set(h, 'EdgeColor', 'none');
caxis([min(real(pf)) max(real(pf))]);
axis off;
view([0, 90]);
cb = colorbar('SouthOutside');
xlabel('$x [m]$','interpreter','latex')
ylabel('$Pressure field Re(p) [Pa]$','interpreter','latex')
title('Total pressure field - analytical solution', 'interpreter', 'latex');

nexttile(3);
hold on
plot_mesh(sphere, real(ps_scat));
h = plot_mesh(field, real(pf_scat));
set(h, 'EdgeColor', 'none');
caxis([min(real(pf)) max(real(pf))]);
axis off;
view([0, 90]);
cb = colorbar('SouthOutside');
xlabel('$x [m]$','interpreter','latex')
ylabel('$Pressure field Re(p) [Pa]$','interpreter','latex')
title('Scattered pressure field - numerical solution', 'interpreter', 'latex');

nexttile(4);
hold on
plot_mesh(sphere, real(ps_scat_a));
h = plot_mesh(field, real(pf_scat_a));
set(h, 'EdgeColor', 'none');
caxis([min(real(pf)) max(real(pf))]);
axis off;
view([0, 90]);
cb = colorbar('SouthOutside');
xlabel('$x [m]$','interpreter','latex')
ylabel('$Pressure field Re(p) [Pa]$','interpreter','latex')
t_res.TileSpacing = 'compact';
t_res.Padding = 'compact';
title('Scattered pressure field - analytical solution', 'interpreter', 'latex');

fig_res.Position = [100 100 800 800];

%% Calculate field at r = 1.5 m
nodes = field.nodes;

radius = 1.5;
tolerance = 0.032; % Allow slight deviations from exact radius

% Calculate distance from origin (assuming circle centered at [0,0,0])
distances = sqrt(nodes(:,1).^2 + nodes(:,2).^2);
circle_node_indices = find(abs(distances - radius) < tolerance);
% Get theta angles (azimuthal positions) of these nodes
theta = atan2(nodes(circle_node_indices, 2), nodes(circle_node_indices, 1));
% Sort nodes by theta for a continuous polar plot
[theta, sort_idx] = sort(theta);
press = abs(pf)-abs(pf_inc);
press_a = abs(pf_a)-abs(pf_inc);
relative_pressures = (press(circle_node_indices(sort_idx)));
relative_pressures_a = (press_a(circle_node_indices(sort_idx)));

figure;
polarplot(theta, relative_pressures, 'k-', 'LineWidth', 1.5);
hold on;
polarplot(theta, relative_pressures_a, 'ro', 'LineWidth', 1.5,'MarkerSize', 2);

% Customize plot
thetaticklabels(0:30:330);
rlim([-1, max(relative_pressures)]);
grid on;

% Add reference line
polarplot(theta, zeros(size(theta)), 'b--', 'LineWidth', 1); 
legend('BEM', 'Analytical','Incident Wave', 'Location', 'northeastoutside');
set(gcf, 'Color', 'w');


