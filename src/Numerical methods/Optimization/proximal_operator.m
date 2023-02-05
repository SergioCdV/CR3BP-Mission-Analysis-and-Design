%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/20
% File: continuation.m 
% Issue: 0 
% Validated: 

%% Proximal operator %%
% This function contains several proximal operators for different norms. 

% Inputs: - vector v, the point for which the proximal operator shall be
%           computed
%         - scalar lambda, the weight of the norm in the minimization
%           problem
%         - string cost_norm, the norm for which the proximal operators shall be
%           computed

% Outputs: - vector x, the arg min of the optimization process 

function [x] = proximal_operator(v, lambda, cost_norm)
    % Branch the different norm solutions 
    switch (cost_norm)
        case 'L2'
            x = max(0, 1-lambda/norm(v))*v;

        case 'L1'
            x = v; 
            for i = 1:size(v,1)
                if (v(i) >= lambda)
                    x(i) = x(i)-lambda;
                elseif (abs(v(i)) <= lambda) 
                    x(i) = 0;
                elseif (v(i) <= -lambda)
                    x(i) = x(i)+lambda;
                end
            end

        case 'Linfty'
            x = v; 
            for i = 1:size(v,1)
                if (v(i) > 1)
                    x(i) = 1;
                elseif (abs(v(i)) <= 1) 
                    x(i) = v(i);
                elseif (v(i) < -1)
                    x(i) = -1;
                end
            end

        otherwise 
            error('The proximal operator for the selected norm is not supported.');
    end

end