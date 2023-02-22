%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Display results %%
% Function to display the results of the optimization

% Inputs: - scalar exitflag, the optimisation exitflag
%         - string cost, indicating the cost function to be minimized
%         - structure output, containing information about the optimisation process

function display_results(exitflag, cost, output)
    % Print the results of the optimisation
    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);

    % Cost results
    fprintf("Final cost: %0.2f\n", cost);
end