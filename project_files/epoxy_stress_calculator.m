% cf = carbon-fiber
% al = aluminum
% ep = epoxy
% Assume SI units for all
% To run successfully, install 'Symbolic Math Toolbox' Add-On in Matlab


% Define material properties
E_al = 70e9; % Young's modulus of aluminum
v_al = 0.33; % Poisson's ratio of aluminum (approx.)

E_cf = 9.76e10; % Young's modulus of carbon-fiber tubing
sigma_comp_cf = 9.65e8; % Compressive yield strength of carbon-fiber tubing
v_cf = 0.06; % Poisson's ratio of aluminum (approx.)

%{
    From Loctite E-120HP properties:
    - Tensile Strength = 41 MPa
    - Elongation = 10%
    Since Young's modulus is not provided, we therefore approximate Young's 
    modulus based on these properties (actual 'E' will be at least this steep)
%}
E_ep = 41e6/.1; % (approx. see above)
v_ep = 0.30; % Poisson's ratio of epoxy (approx.)
G_ep = E_ep/(2*(1+v_ep)); % Shear modulus of epoxy (assuming it is isotropic) (approx.)
%{
    Conservative estimate based on documented shear strength for Loctite E-120HP on aluminum (anodised);
    no documented shear strength (tensile and compression tests with our own carbon-fiber tubing and
    aluminum insert will further define this) 
%}
tau_strength_ep = 14e6;


% Define applied loads and dimensions
N_t = 4058; % input('Enter applied tensile load on aluminum insert: ');
N_c = 3957; % input('Enter applied compression load on aluminum insert: ');
id_cf = .022225; % input('Enter carbon-fiber tubing ID: ');
od_cf = .025273; % input('Enter carbon-fiber tubing OD: ');
l_cf = .1651; % input('Enter carbon-fiber tubing length: ');
h = .0002; % input('Enter bond gap: ');
l_b = 0.0254; % input('Enter bond length: ');


% Calculations
c = l_b/2;
R = id_cf / 2; % Radius at bonding (with cf)
A_cf = pi() * (od_cf^2 - id_cf^2) / 4; % Cross-sectional area of carbon-fiber tubing
A_al = pi() * (id_cf/2 - h)^2; % Cross-sectional area of aluminum insert
I_cf = pi*(od_cf^4 - id_cf^4)/64; % Area moment of inertia for circular tube
% beta = (E_cf*A_cf) / (E_cf*A_cf + E_al*A_al);
k = 2*pi()*R*G_ep/h; % Coefficient in second-order differential equations
% alpha = sqrt(k*(E_cf*A_cf+E_al*A_al)/(E_cf*A_cf*E_al*A_al));
k_r = sqrt(I_cf/A_cf); % Radius of gyration for carbon-fiber tubing
C = 2; % End condition for buckling, specific to our suspension configuration
P_critical = A_cf*C*pi^2*E_cf/(l_cf/k_r)^2; % Critical load before tube becomes unstable (buckles)
sigma_compressive = N_c / A_cf; % Normal stress due to compressive loading

syms y(x)

%{
    Set-up 2nd order diff. eq. of f(x), where f(x) defines the ratio of the tensile load shared by the aluminum
    and cf each (see my Word doc in control arms research folder). To solve this, we must symbolically represent
    it as a 1st order diff. eq., as follows
%}
[V] = odeToVectorField(diff(y,2) - k*(E_cf*A_cf+E_al*A_al)*y/(E_cf*A_cf*E_al*A_al) == -k/(E_al*A_al)); 

M = matlabFunction(V, 'vars', {'x', 'Y'});

%{
Here I set initial conditions where at y(-c) = 1 and y(c) = 0, where y(x) is proportional to the
reaction force felt by the carbon fiber tubing since at y(-c), the span of the aluminum insert
span has ended, after which the load is then fully felt by the carbon-fiber tubing
%}
bc = @(Ya, Yb) [Ya(1) - 1; Yb(1)];

guess = @(x) [(c-x)/(2*c); -1/(2*c)];

solinit = bvpinit(linspace(-c, c, 1000), guess); % Generates the displacement space
sol = bvp4c(M, bc, solinit); % 4th order accuracy

dfdx = @(x) deval(sol, x, 2);

% For the following equations, refer to the research doc I wrote for the derivation
tau = @(x) - N_t/(2*pi()*R) * dfdx(x);

strain_r = @(x) (N_t*R/h)*((v_al*(1-deval(sol, x, 1))/(E_al*A_al))-v_cf*deval(sol, x, 1)/(E_cf*A_cf));
strain_theta = @(x) 0.5*N_t*(v_al*(deval(sol, x, 1)-1)/(E_al*A_al)-v_cf*deval(sol, x, 1)/(E_cf*A_cf));
strain_x = @(x) 0.5*N_t*(deval(sol, x, 1)/(E_cf*A_cf)+(1-deval(sol, x, 1))/(E_al*A_al));

sigma_x = @(x) (1-v_ep)*E_ep*strain_x(x)/((1+v_ep)*(1-2*v_ep)) + v_ep*E_ep*(strain_r(x)+strain_theta(x))/((1+v_ep)*(1-2*v_ep));
sigma_r = @(x) (1-v_ep)*E_ep*strain_r(x)/((1+v_ep)*(1-2*v_ep)) + v_ep*E_ep*(strain_x(x)+strain_theta(x))/((1+v_ep)*(1-2*v_ep));
sigma_theta = @(x) (1-v_ep)*E_ep*strain_theta(x)/((1+v_ep)*(1-2*v_ep)) + v_ep*E_ep*(strain_x(x)+strain_r(x))/((1+v_ep)*(1-2*v_ep));

figure
tiledlayout(2,2)
nexttile
fplot(@(x) deval(sol,x/1000,1), [-c c]*1000)
xlabel('x (mm)')
ylabel('f(x)')
grid on

nexttile 
fplot(@(x) dfdx(x/1000), [-c c]*1000);
xlabel('x (mm)')
ylabel('df/dx (1/mm)')

nexttile
fplot(@(x) tau(x/1000)/1e6, [-c c]*1000)
xlabel('x (mm)')
ylabel('Epoxy Shear Stress (MPa)')

saveas(gcf, "f(x)_and_shear_stress.png");

figure
tiledlayout(3,1)
nexttile
fplot(@(x) sigma_x(x/1000), [-c c]*1000)
xlabel('x (mm)')
ylabel('sigma_x (Pa)')

nexttile
fplot(@(x) sigma_r(x/1000), [-c c]*1000)
xlabel('x (mm)')
ylabel('sigma_r (Pa)')

nexttile
fplot(@(x) sigma_theta(x/1000), [-c c]*1000)
xlabel('x (mm)')
ylabel('sigma_{theta} (Pa)')
saveas(gcf, "normal_stresses.png");

tau_max = max(abs(tau([-c c])));
fprintf('\n\n\tEpoxy Shearing FoS: %.3f \n', tau_strength_ep/tau_max);
fprintf('\tCarbon-Fiber Tubing Compression FoS: %.3f \n', sigma_comp_cf/sigma_compressive);
fprintf('\tEuler Buckling on Carbon-Fiber Tubing FoS: %.3f \n', P_critical/N_c);

