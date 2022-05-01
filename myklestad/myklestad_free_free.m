function [v, wn] = myklestad_free_free(mi, EIi, dx)
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
    %       [EIi] :     Lumped Stifness Vector                         [n-1,1]
    %       [dx]:       Field Length Vector                            [n-1,1]
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
    %
    % Notes
    % ----------
    % For a free-free case, we have something like
    % STATION1--FIELD1--STATION2--FIELD2--STATION3
    % where "|" stands for edges. In this case, we have n stations and n-1 fields

    % Declaring variables
    num_of_stations = length(mi);
    syms w; % Symbolic variable representing natural frequency

    % TODO: prealocate the cells with 4x4 matrices to improve perfomance
    % We need to declare this variables as cell arrays because
    % this is the only way to store matrices in a MATLAB data structure
    % See: https://www.mathworks.com/matlabcentral/answers/496101-matrix-of-matrices-matrix-with-matrices-inside-it
    TF = {}; % Field Transfer Matrices
    TS = {}; % Station Transfer Matrices
    Ts = {}; % Transfer Matrices - cell array containing a transfer matrix at each station

    for i = 1:(num_of_stations - 1)
        % Flexibility influence coefficients
        a_YM = (dx(i)^2) / (2 * EIi(i)); % displacement at i+1 due to a unit moment applied at i+1
        a_YQ = (dx(i)^3) / (3 * EIi(i)); % displacement at i+1 due to a unit force applied at i+1
        a_psiM = dx(i) / EIi(i); % slope at i+1 due to a unit moment applied at i+1
        a_psiQ = (dx(i)^2) / (2 * EIi(i)); % slope at i+1 due to a unit force applied at i+1

        % Transfer Matrix at nearby field
        TF{i} = [1 dx(i) a_YM -a_YQ / 2;
            0 1 a_psiM -a_psiQ;
            0 0 1 -dx(i);
            0 0 0 1];

        % Transfer Matrix at station
        TS{i} = [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            -w^2 * mi(i) 0 0 1];

        % Transfer matrix relating station vector on left side of station i+1
        % to station vector on the left side of station i
        Ts{i} = TF{i} * TS{i};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING OVERALL TRANSFER MATRIX (T)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Overall Transfer Matrix (T = TS{n}*Ts{n-1}*...Ts{2}*Ts{1})
    T = Ts{1};

    for i = 2:(num_of_stations - 1)
        T = Ts{i} * T;
    end

    % TODO: abstract method to calculate a station matrix
    % so we don't repeat this here and inside the loop
    % Final Station Transfer Matrix (TS{n})
    TS{num_of_stations} = [1 0 0 0;
                        0 1 0 0;
                        0 0 1 0;
                        -w^2 * mi(num_of_stations) 0 0 1];

    T = TS{num_of_stations} * T;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING NATURAL FREQUENCIES (wn)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is found by carefull analysis of vRn = T*vL0
    % (where vRn is the station vector at the end edge
    % and vL0 is the station vector at the start edge)
    % Knowing that, for a free edge, Y≠0, psi≠0, M=0, and Q=0
    % we find that the determinant below must be equal to 0
    symbolic_wn = solve(det([T(3, 1) T(3, 2);
                        T(4, 1) T(4, 2)]) == 0, w);

    % Converting the natural frequencies from symbolic to numeric values
    all_wn = double(subs(symbolic_wn)); % This can contain negative numbers
    wn = sort(all_wn(all_wn > 0)); % Keep only positive numbers and sorted from lowest to highest

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING STATE VECTORS (v) FOR MULTIPLE MODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TODO: this can be a parameter of the function
    % or even selecting exactly what modes the user wants
    max_num_of_modes = 3;
    v = {}; % TODO: prealocate the cells with arrays to improve performance

    for i = 1:max_num_of_modes
        % Replacing the symbolic natural frequency with the calculated numeric
        % natural frequency for the overall matrix
        T = double(subs(T, w, wn(i)));

        % We'll have one station vector for every station:
        % If we have: STATION1--FIELD1--STATION2--FIELD2--STATION3
        % We get:     v1----------------v2----------------v3

        % Preallocating the variable that will contain the station vectors
        % with 4x1 vectors
        v{i} = {};

        for j = 1:(num_of_stations)
            v{i}{j} = zeros(4, 1);
        end

        % TODO: determine Y0 and psi0 in the right way. How to arbitrate psi0, instead of 1?
        psi0 = 1; Y0 = -T(3, 2) * psi0 / T(3, 1); M0 = 0; Q0 = 0;
        v{i}{1} = [Y0; psi0; M0; Q0]; % station vector at start/left edge (wall)

        for j = 2:(length(v{i}))
            % Creating a transfer matrix at station with symbolic natural frequency
            % converted to numeric natural frequency
            Ts_iminus1 = double(subs(Ts{j - 1}, w, wn(i)));
            v{i}{j} = Ts_iminus1 * v{i}{j - 1};
        end

    end

end
