function [wn] = kumai(B, D, T, Lpp, displacement)
    % Calculate the natural frequencies by using the simple Kumai's formula
    % B = Beam
    % D = Depth
    % T = Draft
    % Lpp = Length between perpendiculars
    % displacement = displacement of ship

    % Area moment of inertia of ship's midsection [m^4]
    I = B * (D^3) / 12;

    % Displacement including virtual mass [t]
    displacement_with_vmass = (1.2 + B / (3 * T)) * displacement;

    % First natural frequency [rad/s]
    wn(1) = (2 * pi / 60) * (3.07 * 10^6) * sqrt(I / (displacement_with_vmass * Lpp^3));

    num_of_natural_modes = 3;
    alpha = 0.845; % considering a tanker ship

    for n = 2:num_of_natural_modes
        wn(n) = wn(1) * (n - 1)^alpha;
    end

end
