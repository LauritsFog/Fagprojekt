function [U, Y] = simulatePID(tk, xk, yk, dk, Nk, p, ctrlPar, ctrlStatek, ctrlAlgorithm, simModel, simMethod, observationModel, N, tzero, haltingiter, rampingfunction)
    
    U = zeros(2,N);
    Y = zeros(1,N);
    Y(1) = yk;
    
    Ts = ctrlPar(1);
    
    tpause = haltingiter;
    
    for k = 1:N
        % Times
        tkp1 = tk+Ts;
        
        % Compute manipulated inputs
        [uk, ctrlStatekp1] = ctrlAlgorithm(tk, yk, dk, ctrlPar, ctrlStatek, tzero, tpause, haltingiter, rampingfunction);
        
        if tpause ~= 0
            tpause = tpause-1;
        end
        
        % Time interval
        tspank = linspace(tk, tkp1, Nk+1);

        % Solve initial value problem
        [Tk, Xk] = simMethod(simModel, tspank, xk, [], uk, dk, p);

        % States at the next time step
        tkp1 = Tk(end   );
        xkp1 = Xk(end, :)';

        % Observed variables at the next time step
        ykp1 = observationModel(tkp1, xkp1, p);

        % Store solution
        U(:, k) = uk;
        Y(:, k+1) = ykp1;
        ctrlStatek = ctrlStatekp1;

        % Update initial condition
        tk = tkp1;
        xk = xkp1;
        yk = ykp1;
        
        % Setting disturbance to 0 after intitial effect
        dk = 0;
    end
end
    
        