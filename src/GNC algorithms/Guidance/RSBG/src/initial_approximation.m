%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - sampling_distribution, string specifying the sampling grid
%           distribution
%         - structure St, inidicating the target spacecraft position
%           evolution
%         - scalar tfapp, the initial approximation of the time of flight
%         - vector tau, the collocation points to be used 
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collocation basis
%         - string dynamics, indicating the dynamics vectorfield
%           formulation to be used

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - array Capp, the initial estimation of the spacecraft state vector
%          - scalar Napp, the estimated number of revolutions needed
%          - scalar tfapp, the initial initial time of flight

function [Papp, Capp, Napp, tfapp] = initial_approximation(mu, St, tfapp, tau, initial, final, basis, dynamics)
    % Preliminary number of revolutions
    dtheta = final(2)-initial(2);
    if (dtheta < 0)
        dtheta = dtheta + 2*pi; 
    end
    
    Napp = ceil( (dtheta+tfapp*0.5*(initial(4)+final(4)) ) / (2*pi) );
    if (Napp > 0)
        % New initial TOF
        tfapp = tfapp*Napp;
    end 

    % Generate the polynomial basis
    n_init = repmat(3, [1 3]);
    Bapp = state_basis(n_init, tau, basis);

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(length(initial)/2, max(n_init)+1);  
    Papp = boundary_conditions(tfapp, n_init, initial, final, Napp, Papp, Bapp, basis);

    % State vector approximations
    Capp = evaluate_state(Papp, Bapp, n_init);

    % Time-regularized solution 
    switch (dynamics)
        case 'Sundman'
            % Arc-length regularization or generalized Sundman transformation
            r = sundman_radius(mu, tfapp, St, tau, Capp);
            tfapp = tfapp*trapz(tau, r.^(-1));  
            Papp = boundary_conditions(tfapp, n_init, initial, final, Napp, Papp, Bapp, basis);
        
            % State vector approximations
            Capp = evaluate_state(Papp, Bapp, n_init);
        otherwise
    end
end