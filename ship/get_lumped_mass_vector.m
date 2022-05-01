function [mi] = get_lumped_mass_vector(x)
    % x is a vector with the position of each station
    % mi is the vector with corresponding lumped masses for each station

    % Position of stations
    x_of_stations = [0; 7.88; 15.75; 23.63; 31.50; 39.38; 47.25; 55.13; 63.00; 70.88; ...
                    78.75; 86.83; 94.50; 102.38; 110.25; 118.13; 126.00; 133.88; 141.75; ...
                    149.63; 157.50];

    % Mass distribution [kg/m]
    W_per_station = [17; 38; 59; 73; 88; 103; 117; 132; 146; 146; 146; 146; 146; ...
                    135; 123; 111; 99; 88; 76; 51; 27] .* 1000;

    % Regression to obtain the polynomial responsible for describing weigth per length against x
    W_per_length_fit = fit(x_of_stations, W_per_station, 'linearinterp'); % Weight per length [kg/m]

    mi = zeros(length(x), 1);

    for i = 1:(length(x) - 1)
        xx = x(i:i + 1);
        yy = W_per_length_fit(xx);
        mi(i) = trapz(xx, yy);
    end

end
