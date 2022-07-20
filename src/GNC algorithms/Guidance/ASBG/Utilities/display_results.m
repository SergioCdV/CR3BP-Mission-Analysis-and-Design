%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Display results %%
% Function to display the results of the optimization

% Inputs: - scalar exitflag, the optimisation exitflag
%         - structure output, containing information about the optimisation process 
%         - scalar r0, the fundamental length unit
%         - scalar t0, the fundamental time unit
%         - scalar tfapp, the initial estimated time of flight
%         - scalar tf, the final computed time of flight 
%         - scalar dV, the final optimal cost

function display_results(exitflag, output, r0, t0, tfapp, tf, dV)
    % Constants
    days2sec = t0/86400;

    % Print the results of the optimisation
    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);

    % Time of flight results
    fprintf("Initial estimation of flight time: %0.2f days\n", tfapp*days2sec);
    fprintf("Final calculation of flight time: %0.2f days\n", tf*days2sec);

    % Cost results
    fprintf("Final cost %0.2f m/s\n\n", dV*r0/t0);
end