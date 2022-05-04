function [v, wn] = mykelstad_clamped2(mi, EI, dx)
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
    %       [mi] :      Lumped Masses Vector                           [n,1]
    %       [EI] :      Lumped Stifness Vector                         [n+1,1]
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initializing Variables
    num_of_stations = length(mi);
    syms w;
    TF = {}; % Field Transfer Matrices
    T_S = {}; % Station Transfer Matrices
    T_s_cellArray = {}; % Transfer Matrices - cell array containing a transfer matrix at each station

    % |--FIELD0--STATION1--FIELD1--...STATIONn--FIELDn--|

    % Field 0 will be created dependending on boundary conditions
    % so the loop start from i = 2
    for i = 1:(num_of_stations)
        % Flexibility influence coefficients
        a_YM = (dx(i)^2) / (2 * EI(i)); % displacement at i+1 due to a unit moment applied at i+1
        a_YQ = (dx(i)^3) / (3 * EI(i)); % displacement at i+1 due to a unit force applied at i+1
        a_psiM = dx(i) / EI(i); % slope at i+1 due to a unit moment applied at i+1
        a_psiQ = (dx(i)^2) / (2 * EI(i)); % slope at i+1 due to a unit force applied at i+1

        % Transfer Matrix at field i
        TF{i} = [1 dx(i) a_YM -a_YQ / 2;
            0 1 a_psiM -a_psiQ;
            0 0 1 -dx(i);
            0 0 0 1];

        % Transfer Matrix at station i-1
        T_S{i} = [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            -w^2 * mi(i) 0 0 1];

        % Transfer matrix relating station vector on left side of station i+1
        % to station vector on the left side of station i
        T_s_cellArray{i} = TF{i} * T_S{i};
    end

    % CALCULATING OVERALL TRANSFER MATRIX (T)

    % Initial Field Transfer Matrix
    a_YM = (dx(1)^2) / (2 * EI(1));
    a_YQ = (dx(1)^3) / (3 * EI(1));
    a_psiM = dx(1) / EI(1);
    a_psiQ = (dx(1)^2) / (2 * EI(1));
    TF{1} = [1 dx(1) a_YM -a_YQ / 2;
        0 1 a_psiM -a_psiQ;
        0 0 1 -dx(1);
        0 0 0 1]; % Field Matrix at first field

    % Overall Transfer Matrix (T = T_s_cellArray{n}*T_s_cellArray{n-1}*...T_s_cellArray{2}*T_s_cellArray{1}*TF{1})
    T = TF{1};

    for i = 1:length(T_s_cellArray)
        T = T_s_cellArray{i} * T;
    end

    %% CALCULATING NATURAL FREQUENCIES (wn)

    % This is found by carefull analysis of vRn = T*vL0
    % (where vRn is the station vector at the end edge
    % and vL0 is the station vector at the start edge)
    % Knowing that, for a clamped edge, Y=0, psi=0, M≠0, and Q≠0
    % we find that the determinant below must be equal to 0
    symbolic_wn = solve(det([T(1, 3) T(1, 4);
                        T(2, 3) T(2, 4)]) == 0, w);

    % Converting the natural frequencies from symbolic to numeric values
    all_wn = double(subs(symbolic_wn)); % This can contain negative numbers
    wn = all_wn(all_wn > 0); % Keep only positive numbers

    % Storing only the three first natural frequencies
    wn_sorted = sort(wn);
    wn = wn_sorted(1:3);

    % CALCULATING STATE VECTORS (v) FOR 3 FIRST MODES
    max_num_of_modes = 3;
    v = {};

    for i = 1:max_num_of_modes
        T = double(subs(T, w, wn(i)));

        % We'll have one station vector for every station + 2 station vectors:
        % one at the start/left edge, and one at the end/right edge
        % If we have: |--FIELD0--STATION1--FIELD1--STATION2--FIELD2--|
        % We get:     v0-------v1-----------------v2----------------v3

        % Preallocating the variable that will contain the station vectors
        % with 4x1 vectors

        v{i} = {};

        for j = 1:(num_of_stations + 1)
            v{i}{j} = zeros(4, 1);
        end

        % TODO: determine M0 and Q0 in the right way. How to arbitrate M0, instead of 0?
        Y0 = 0; psi0 = 0; M0 = -1e8; Q0 = -T(1, 3) * M0 / T(1, 4);
        v{i}{1} = [Y0; psi0; M0; Q0]; % station vector at start/left edge (wall)
        v{i}{2} = TF{1} * v{i}{1}; % station vector at first station

        for j = 3:(length(v{i}))
            % Creating a transfer matrix at station with symbolic natural frequency
            % converted to numeric natural frequency

            Ts_iminus1 = double(subs(T_s_cellArray{j - 1}, w, wn(i)));
            v{i}{j} = Ts_iminus1 * v{i}{j - 1};
        end

        %v{i}{length(v)} = T * v{i}{1}; % station vector at end/right edge (wall)
    end

end
