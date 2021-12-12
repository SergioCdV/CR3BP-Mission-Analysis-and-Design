%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/20
% File: dimensionalizer.m 
% Issue: 0 
% Validated: 

%% Dimensionalizer %%
% For a given gravitational parameter mu, distance between the primaries d and synodic period T, 
% this function dimensionalizes/normalizes position or velocity vectors and
% epochs.

% Inputs: - scalar V, the characteristic velocity of the system.
%         - scalar d, the characteristic distance of the system.
%         - scalar T, the synodic period of the system.
%         - string X, an epoch, position or velocity vector.
%         - boolean direction, 0 to normalize and 1 to dimensionalize.

% Outputs: - output sol, the normalized variable.

% New versions: 

function [sol] = dimensionalizer(d, T, V, X, magnitude, direction)
    %Main computation
    if (direction == 0)
        %Normalize input
        switch (magnitude)
            case 'Epoch'
                sol = X*(2*pi/T);
            case 'Position'
                sol = X/d;
            case 'Velocity'
                sol = X/V;
            otherwise
                error('No valid magnitude was selected');
        end
            
    elseif (direction == 1)
        %Dimensionalize input
        switch (magnitude)
            case 'Epoch'
                sol = X*(T/(2*pi));
            case 'Position'
                sol = X*d;
            case 'Velocity'
                sol = X*V;
            otherwise
                error('No valid magnitude was selected');
        end
        
    %Error branch
    else
        error('No valid direction was selected');
    end
end