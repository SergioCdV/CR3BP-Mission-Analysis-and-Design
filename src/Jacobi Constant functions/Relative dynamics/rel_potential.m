%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/06/21
% File: rel_potential.m 
% Issue: 0 
% Validated: 

%% Relative potential %%
% For a given gravitational parameter mu and position vector r, this function computes the 
% potential function associated with that input position vector

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - vector s, a 12x1 array containing the relative synodic position 
%           velocity vectors of the chaser and the target

% Outputs: - scalar U, the potential function associated with the input position vector 

% New versions:

function [U] = rel_potential(mu, s)
    % State variables 
    rt = s(1:3);                        % Position vector of the target mass
    rho = s(7:9);                       % Position vector of the relative particle
    
    % Constants of the problem 
    mu_r(1) = 1-mu;                     % Gravitational parameters of the most massive primary
    mu_r(2) = mu;                       % Gravitational parameters of the least massive primary
    
    % Location of the unsteady primaries 
    R(1:3,1) = [-mu; 0; 0];             % Location of the first primary
    R(1:3,2) = [1-mu; 0; 0];            % Location of the second primary
    Rr(1:3,1) = R(:,1)-rt;              % Location of the unsteady first primary
    Rr(1:3,2) = R(:,2)-rt;              % Location of the unsteady second primary
        
    % Relative potential function 
    U = mu_r.*(1./sqrt(dot(rho-Rr, rho-Rr,1))-rho.'*Rr./sqrt(dot(Rr, Rr,1)));
    U = -sum(U);
end