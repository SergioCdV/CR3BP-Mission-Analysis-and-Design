%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/01/23
% File: ISP_control.m 
% Issue: 0 
% Validated: 25/01/23

%% Impuslive Sequence Pruner %%
% This script contains the function to reduce the control law by means of Lawden's pruner.

% Inputs: - 

% Output: - array dV, the pruned impulses sequence 
%         - scalar cost, the associated L2 cost of the sequence

% New versions: 

function [dV, cost] = ISP_control(Phi, B, dV, J)
    % Constants 
    m = size(B,1);                              % Dimension of the state space
    num_sequence = size(dV,2)-(size(B,1)+1);    % Number of sequence reductions

    % Initial setup
    if (num_sequence >= 0)
        cost = J;                     % Current cost 
        sequence = dV(:,1:m+1);       % Initial sequence

        for i = 0:num_sequence
         % Reduce the sequence 
         [dv, cost] = sequence_reduction(Phi(1:(m+1+i),:), B, sequence, cost);
    
            % New sequence 
            if (i < num_sequence)
                sequence = [dv dV(:,m+1+i)];
            end
        end
    
        % Final sequence 
        dV = dv;
    else
        cost = J; 
    end
end

%% Auxiliary functions 
function [dV, cost] = sequence_reduction(Phi, B, dV, J)
    % Constants 
    m = size(B,1);                              % Dimension of the state space

    % Weights of the sequence 
    V = sqrt(dot(dV,dV,1));

    % Final matrix
    M = reshape(Phi(end,:), [m m]);

    u = zeros(m, size(dV,2));
    for i = 1:size(dV,2)
        R = M*reshape(Phi(i,:), [m m])^(-1)*B;
        u(:,i) = R*dV(:,i);
    end

    u(:, V ~= 0) = u(:, V ~= 0)./V(V ~= 0);

    index = V == 0;
    dumb = u(:,~index);

    if (~isempty(dumb))
        b = -dumb(:,end);
        dumb = [(dumb(:,1:end-1)\b).' 1];
    
        alpha = zeros(1,size(u,2));
        alpha(~index) = dumb;
        Alpha = sum(alpha);
    
        if (Alpha < 0)
            alpha = -alpha;
        end
    
        if (any(isnan(alpha)))
            a = 1; 
        end
    
        beta = alpha(V ~= 0)./V(V ~= 0);
        beta_r = max(beta);
        mu = V(V ~= 0)./beta_r.*(beta_r-beta);
    
        % New impulse sequence 
        dV(:,V ~= 0) = (mu./V(V ~= 0)).*dV(:,V ~= 0); 
    
        % New sequence cost 
        cost = J-Alpha/beta_r;
    else
        cost = V; 
    end
end

