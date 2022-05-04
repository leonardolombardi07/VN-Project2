function [v, wn] = myklestad_clamped_clamped(mi, EIi, dx)
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
    %       [EIi] :     Lumped Stifness Vector                          [n+1,1]
    %       [dx]:       Field Length Vector                             [n+1,1]
    %
    % Output
    % ----------
    %       [v]:        Station Vectors                        cell array [n,3]
    %                    - three vectors (for the first three natural
    %                    frequencies) containing station vectors for that given
    %                    frequency
    %                    - each station vector has the form [Y, psi, M, Q],
    %                    where "Y" stands for translational displacement, "psi"
    %                    for angular displacement, "M" bending moment and "Q"
    %                    shearing force
    %
    %       [wn]:       Natural Frequencies Vector                        [n,1]
    %
    % Notes
    % ----------
    % For a clamped-clamped case, we have something like
    % |--FIELD0--STATION1--FIELD1--STATION2--FIELD2--|,
    % where "|" stands for edges. In this case, we have n stations and n+1 fields

    % Declaring variables
    num_of_stations = length(mi);
    syms w; % Symbolic variable representing natural frequency

    % TODO: prealocate the cells with 4x4 matrices to improve perfomance
    % We need to declare this variables as cell arrays because
    % this is the only way to store matrices in a MATLAB data structure
    % See: https://www.mathworks.com/matlabcentral/answers/496101-matrix-of-matrices-matrix-with-matrices-inside-it
    TF = {}; % Field Transfer Matrices
    TS = {}; % Station Transfer Matrices
    T = {}; % Transfer Matrices

    % TODO: make this comment clearer or find a better way to access
    % data from fields
    % We start counting at i=2 because we skip the first field, which is not used
    % to calculate Ts{i}. Remember, for a clamped-clamped case, we have:
    % |--FIELD0--STATION1--FIELD1--STATION2--FIELD2--|
    % For STATION1, we calculate the transfer matrix as Ts{1} = TF{2}*Ts{1}
    % For STATION2, we calculate the transfer matrix as Ts{2} = TF{3}*Ts{2}
    % So there's no need to access FIELD0 in the loop
    for i = 2:(num_of_stations + 1)
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
        TS{i - 1} = [1 0 0 0;
                0 1 0 0;
                0 0 1 0;
                -w^2 * mi(i - 1) 0 0 1];

        % Transfer matrix relating station vector on left side of station i+1
        % to station vector on the left side of station i
        T{i - 1} = TF{i} * TS{i - 1};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING OVERALL TRANSFER MATRIX (T_overall)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TODO: abstract method to calculate a field matrix
    % so we don't repeat this here and inside the loop
    % Initial Field Transfer Matrix
    a_YM = (dx(1)^2) / (2 * EIi(1));
    a_YQ = (dx(1)^3) / (3 * EIi(1));
    a_psiM = dx(1) / EIi(1);
    a_psiQ = (dx(1)^2) / (2 * EIi(1));
    TF{1} = [1 dx(1) a_YM -a_YQ / 2;
        0 1 a_psiM -a_psiQ;
        0 0 1 -dx(1);
        0 0 0 1]; % Field Matrix at first field

    % Overall Transfer Matrix (T_overall = T{n}*T{n-1}*...T{2}*T{1}*TF{1})
    % Unfortunately, in math books the overall transfer matrix is also
    % represented by the letter "T". Here, we are a bit more explicit naming it
    % "T_overall"
    T_overall = TF{1};

    for i = 1:length(T)
        T_overall = T{i} * T_overall;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING NATURAL FREQUENCIES (wn)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is found by careful analysis of vRn = T_overall*vL0
    % (where vRn is the station vector at the end edge
    % and vL0 is the station vector at the start edge)
    % Knowing that, for a clamped edge, Y=0, psi=0, M≠0, and Q≠0
    % we find that the determinant below must be equal to 0
    symbolic_wn = solve(det([T_overall(1, 3) T_overall(1, 4);
                        T_overall(2, 3) T_overall(2, 4)]) == 0, w);

    % Converting the natural frequencies from symbolic to numeric values
    all_wn = double(subs(symbolic_wn)); % This can contain negative numbers
    wn = sort(all_wn(all_wn > 0)); % Keep only positive numbers, sorted from lowest to highest

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
        T_overall_for_wn = double(subs(T_overall, w, wn(i)));

        % We'll have one station vector for every station + 2 station vectors:
        % one at the start/left edge, and one at the end/right edge
        % If we have: |--FIELD0--STATION1--FIELD1--STATION2--FIELD2--|
        % We get:     v0-------v1-----------------v2----------------v3

        % Preallocating the variable that will contain the station vectors
        % with 4x1 vectors

        v{i} = {};

        for j = 1:(num_of_stations + 2)
            v{i}{j} = zeros(4, 1);
        end

        % TODO: determine M0 and Q0 in the right way. How to arbitrate M0, instead of 1?
        % Q0 is also found by analysis of vRn = T_overall*vL0
        Y0 = 0; psi0 = 0; M0 = 1;
        Q0 = -T_overall_for_wn(1, 3) * M0 / T_overall_for_wn(1, 4);
        v{i}{1} = [Y0; psi0; M0; Q0]; % station vector at start/left edge (wall)
        v{i}{2} = TF{1} * v{i}{1}; % station vector at first station

        for j = 3:(length(v{i}) - 1)
            % Replacing the symbolic natural frequency with the calculated numeric
            % natural frequency for the transfer matrix at given station
            T_numeric = double(subs(T{j - 2}, w, wn(i)));
            v{i}{j} = T_numeric * v{i}{j - 1};
        end

        v{i}{length(v{i})} = T_overall_for_wn * v{i}{1}; % station vector at end/right edge (wall)
    end

end
