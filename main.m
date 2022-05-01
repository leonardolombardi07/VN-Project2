%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning Workspace & Importing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; % Closes all figures
clear; % Clears variables from workspace
clc; % Clears all text from command window

% Importing functions from folders with given relative path
addpath(genpath("./analytical"));
addpath(genpath("./myklestad"));
addpath(genpath("./kumai"));
addpath(genpath("./ship"));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants to be reused throughout this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = 200e9; % Young Modulus [N/m^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1) Validation of the Mykelstad method by considering a uniform beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Beam constants
L_beam = 30; % Length of beam [m]
p_beam = 7850; % Density of beam [kg/m^3]
B_beam = 1.2; % Width of beam [m]
H_beam = 0.4; % Thickness of beam [m]
A_beam = B_beam * H_beam; % Sectional Area of beam [m^2]
m_beam = p_beam * A_beam * L_beam; % Mass of beam [kg]
I_beam = B_beam * (H_beam^3) / 12; % Area moment of inertia of beam section [m^4]
EI_beam = E * I_beam; % [N.m^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 a) Solve the uniform beam problem of (a) with the Mykelstad
% method. Carry out a convergence analysis. Calculate the first three
% bending natural frequencies and Eigen functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1.a - I) Solving the uniform beam problem of (a) with the Mykelstad method
% 2.1.a - II) Convergence analysis
% 2.1.a - III) Calculating the first three bending natural frequencies and Eigen functions

% Vector with increasing numbers of stations, for convergence analysis
convergent_num_of_stations = [10];

% for i = 1:length(convergent_num_of_stations)
% num_of_stations = convergent_num_of_stations(i);
num_of_stations = 10;
num_of_fields = num_of_stations + 1;

% Field length discretization
dx = (L_beam / (num_of_fields)); % Length of Field [m]
dxs = dx * ones(num_of_fields, 1); % Field Length Vector

% EIs discretization
EIis = EI_beam * ones(num_of_fields, 1); % Lumped Stifness Vector

% Mass discretization
mi = m_beam / num_of_stations; % Mass per station [kg/st]
mis = mi * ones(num_of_stations, 1); % Lumped Masses Vector

% [v, wn] = myklestad_clamped_clamped(mis, EIis, dxs);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 b) Comparing the above results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical results for comparison
x = 0:dx:L_beam;
[v_analytical_beam, wn_analytical_beam] = analytical(x, EI_beam, p_beam, A_beam, L_beam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2) Calculating the vertical flexural vibrations of a ship using
% Mykelstad method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ship constants
Lpp = 157.5; % Length between perpendiculars [m]
LOA = 40; % Length over all [m]
B = 22.86; % Beam [m]
D = 11; % Depth [m]
T = 8.55; % Draft [m]
displacement = 18036.9; % [t]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 a) Calculating the first three natural frequencies and modes for the
% non-uniform ship vertical flexural vibrations with Mykelstad method.
% Carry out a convergence analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.2.a - I) Calculating the first three natural frequencies and modes for the
% non-uniform ship vertical flexural vibrations with Mykelstad method
% 2.2.a - II) Convergence analysis

% Vector with increasing numbers of stations, for convergence analysis
convergent_num_of_stations = [10];

% for i = 1:length(convergent_num_of_stations)
% num_of_stations = convergent_num_of_stations(i);
num_of_stations = 10;
num_of_fields = num_of_stations - 1;

% Field length discretization. Considering a uniform ship here
dx = (Lpp / (num_of_fields)); % Length of Field [m]
dxs = dx * ones(num_of_fields, 1); % Field Length Vector

% Vector with the position of every station
x = 0:dx:Lpp;

% EIs discretization
% EIis = get_lumped_stiffness_vector(x);

% Mass discretization
% mi = get_lumped_mass_vector(x);
% added_mi = get_added_mass_vector(x);
% mis = mi + added_mi;

% [v_ship, wn_ship] = myklestad_free_free(mis, EIis, dxs);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 b) Calculating the first three natural frequencies by using a simple
% formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wn = kumai(B, D, T, Lpp, displacement)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 c) Compare the estimations of natural frequency and natural modes
% from the FEM method (1.3), the results from Mykelstad method (2.2a) and
% the empirical calculation (2.2b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
