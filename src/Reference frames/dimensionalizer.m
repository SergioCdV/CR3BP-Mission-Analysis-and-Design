%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/20
% File: dimensionalizer.m 
% Issue: 0 
% Validated: 

%% Dimensionalizer %%
% For a given gravitational parameter mu, distance between the primaries d and synodic period T, 
% this function dimensionalizes/normalizes position or velocity vectors and
% epochs. These variables are imported from a .txt file in the system_constants.m function.

% Inputs: - scalar V, the most massive primary orbital velocity.
%         - scalar d, the distance between primaries.
%         - scalar T, the synodic period of the system.
%         - input X, an epoch, position or velocity vector, specified by
%           string s.
%         - boolean direction, 0 to normalize and 1 to dimensionalize.

% Outputs: - output sol.

% New versions: 

function [sol] = dimensionalizer(d, T, V, X, s, direction)
    %Main computation
    if (direction == 0)
        %Normalize input
        switch (s)
            case 'Epoch'
                sol = X*(2*pi/T);
            case 'Position'
                sol = X/d;
            case 'Velocity'
                sol = X/V;
            otherwise
                sol = [];
        end
            
    elseif (direction == 1)
        %Dimensionalize input
        switch (s)
            case 'Epoch'
                sol = X*(T/(2*pi));
            case 'Position'
                sol = X*d;
            case 'Velocity'
                sol = X*V;
            otherwise
                sol = [];
        end
        
    %Error branch
    else
        sol = [];
    end
end