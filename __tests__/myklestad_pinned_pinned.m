function [v, wn] = myklestad_pinned_pinned(mi, EIi, dx)
    % Myklestad Method
    %--------------------------------------------------------------------------
    % Iterates through lumped masses over stations and their nearby fields,
    % relating displacement and force at the right end of a station to those at
    % the left end of the same station to find out station vectors
    % Returns a list of station vectors (one for every station) and a vector
    % containing the natural frequencies of the system
    %
    % Input
    % ----------
    %       [mi] :      Lumped Masses Vector                              [n,1]
    %       [EIi] :     Lumped Stifness Vector                         [n+1,1]
    %       [dx]:       Field Length Vector                            [n+1,1]
    %
    % Output
    % ----------
    %       [v]:        Station Vectors                                   [n,1]
    %                    - each station vector has the form [Y, psi, M, Q],
    %                    where "Y" stands for translational displacement, "psi"
    %                    for angular displacement, "M" bending moment and "Q"
    %                    shearing force
    %
    %       [wn]:       Natural Frequencies Vector                        [n,1]

    num_of_stations = length(mi);
    syms w; % Symbolic variable representing natural frequency

    TF = {}; % Field Transfer Matrices
    TS = {}; % Station Transfer Matrices
    Ts = {}; % Transfer Matrices - cell array containing a transfer matrix at each station

    for i = 2:(num_of_stations + 1)
        % Flexibility influence coefficients
        a_YM = (dx(i)^2) / (2 * EIi(i)); % displacement at i+1 due to a unit moment applied at i+1
        a_YQ = (dx(i)^3) / (3 * EIi(i)); % displacement at i+1 due to a unit force applied at i+1
        a_psiM = dx(i) / EIi(i); % slope at i+1 due to a unit moment applied at i+1
        a_psiQ = (dx(i)^2) / (2 * EIi(i)); % slope at i+1 due to a unit force applied at i+1

        TF{i} = [1 dx(i) a_YM -a_YQ / 2;
            0 1 a_psiM -a_psiQ;
            0 0 1 -dx(i);
            0 0 0 1];

        TS{i - 1} = [1 0 0 0;
                0 1 0 0;
                0 0 1 0;
                -w^2 * mi(i - 1) 0 0 1];

        Ts{i - 1} = TF{i} * TS{i - 1};
    end

    a_YM = (dx(1)^2) / (2 * EIi(1));
    a_YQ = (dx(1)^3) / (3 * EIi(1));
    a_psiM = dx(1) / EIi(1);
    a_psiQ = (dx(1)^2) / (2 * EIi(1));
    TF{1} = [1 dx(1) a_YM -a_YQ / 2;
        0 1 a_psiM -a_psiQ;
        0 0 1 -dx(1);
        0 0 0 1]; % Field Matrix at first field

    T = TF{1};

    for i = 1:length(Ts)
        T = Ts{i} * T;
    end

    symbolic_wn = solve(det([T(1, 2) T(1, 4);
                        T(3, 2) T(3, 4)]) == 0, w);

    all_wn = double(subs(symbolic_wn)); % This can contain negative numbers
    wn = sort(all_wn(all_wn > 0)); % Keep only positive numbers, sorted from lowest to highest

    T = double(subs(T, w, wn(1)));

    v = {}; % TODO: prealocate the cells with arrays to improve performance
    Y0 = 0; psi0 = 1; M0 = 0; Q0 = -T(1, 2) * psi0 / T(1, 4);
    v{1} = [Y0; psi0; M0; Q0]; % v0 station vector at start/left edge (wall)
    v{2} = TF{1} * v{1}; % v1
    v{3} = double(subs(Ts{1}, w, wn(1))) * v{2}; % v2
    v{4} = double(subs(Ts{2}, w, wn(1))) * v{3}; % v3
    v{5} = double(subs(Ts{3}, w, wn(1))) * v{4}; % v4
    v{6} = double(subs(Ts{4}, w, wn(1))) * v{5}; % v5
    v{7} = double(subs(Ts{5}, w, wn(1))) * v{6}; % v6
    v{8} = double(subs(Ts{6}, w, wn(1))) * v{7}; % v7
    v{9} = double(subs(Ts{7}, w, wn(1))) * v{8}; % v8
    v{10} = double(subs(Ts{8}, w, wn(1))) * v{9}; % v9
    v{11} = double(subs(Ts{9}, w, wn(1))) * v{10}; % v10
    v{12} = T * v{1}; % station vector at end/right edge (wall)
end
