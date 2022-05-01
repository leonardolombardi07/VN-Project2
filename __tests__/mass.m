addpath("../ship");

Lpp = 157.5; % Length between perpendiculars [m]

num_of_stations = 10
num_of_fields = num_of_stations - 1;

% Field length discretization. Considering a uniform ship here
dx = (Lpp / (num_of_fields)); % Length of Field [m]
dxs = dx * ones(num_of_fields, 1); % Field Length Vector

% Vector with the position of every station
x_ship = 0:dx:Lpp;

% EIs discretization
EIis = get_lumped_stiffness_vector(x_ship);

% Mass discretization
mis = get_lumped_mass_vector(x_ship); % we consider the addded mass here
