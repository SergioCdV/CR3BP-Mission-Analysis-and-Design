%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: CoPotentialFunction.m 
% Issue: 0 
% Validated: 

%% Co-orbital potential function %%
% For a given gravitational parameter mu and position vector r, this function computes the 
% potential function associated with that input position vector

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - array s, a 12xN array containing the relative synodic position 
%           velocity vectors of the chaser and the target

% Outputs: - vector U, the potential function associated with the input position vector 

% New versions:

function [U] = CoPotentialFunction(mu, s)
    % State variables 
    rt = s(1:3,:);                      % Position vector of the target mass
    rho = s(7:9,:);                     % Position vector of the relative particle
    
    % Constants of the problem 
    mu_r(1) = 1 - mu;                   % Gravitational parameters of the first primary
    mu_r(2) = mu;                       % Gravitational parameters of the secondary primary
    
    % Location of the unsteady primaries 
    R(1:3,1) = [-mu; 0; 0];             % Location of the first primary
    R(1:3,2) = [1 - mu; 0; 0];          % Location of the second primary
    Rr(1:3,:) = R(:,1) - rt;            % Location of the relative first primary
    Rr(4:6,:) = R(:,2) - rt;            % Location of the relative second primary

    Rr_norm(1,:) = sqrt( dot(Rr(1:3,:), Rr(1:3,:), 1) );
    Rr_norm(2,:) = sqrt( dot(Rr(4:6,:), Rr(4:6,:), 1) );
        
    % Relative potential function 
    U = zeros(2, size(s,2));
    for i = 1:length(mu_r)
        idx = 1 + 3 * (i - 1): 3 * i; 
        rel_pos = rho - Rr(idx,:);
        U(i,:) = mu_r(i) * ( 1 ./ sqrt( dot(rel_pos, rel_pos, 1) ) - dot(rho, Rr(idx,:), 1) ./ Rr_norm(i,:) );
    end

    U = -sum(U, 1);
end