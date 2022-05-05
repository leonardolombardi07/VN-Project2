function [Y] = plot_convergence_analysis(convergent_num_of_stations, L)

    for i = 1:length(convergent_num_of_stations)
        num_of_stations = convergent_num_of_stations(i);
        num_of_fields = num_of_stations - 1;

        % Field length discretization. Considering a uniform ship here
        dx = (L / (num_of_fields)); % Length of Field [m]
        dxs = dx * ones(num_of_fields, 1); % Field Length Vector

        % Vector with the position of every station
        x = 0:dx:L;

        % EIs discretization
        EIis = get_lumped_stiffness_vector(x);

        % Mass discretization
        mis = get_lumped_mass_vector(x); % we consider the addded mass here

        [v_ship, wn_ship] = myklestad_free_free(mis, EIis, dxs);

    end

    Y1 = extract_mode_from_station_vector(v_ship{1});
    Y2 = extract_mode_from_station_vector(v_ship{2});
    Y3 = extract_mode_from_station_vector(v_ship{3});

    figure;
    subplot(3, 1, 1);
    xlabel("x [m]")
    ylabel("Y1")
    title("Mode Vector for first ship natural frequency")
    plot(x, Y1);

    subplot(3, 1, 2);
    plot(x, Y2);
    xlabel("x [m]")
    ylabel("Y1")
    title("Mode Vector for second ship natural frequency")

    subplot(3, 1, 3);
    plot(x, Y3);
    xlabel("x [m]")
    ylabel("Y1")
    title("Mode Vector for third ship natural frequency")
    return;

end
