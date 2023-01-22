%% Project: 
% Date: 30/01/2022

%% Orthogonal Bernstein derivative %% 
% Function to compute the basis of the tangent space of the orthogonal Bernstein polynomials 

% Inputs:  - scalar degree, the cardinal of the polynomial basis
%          - vector tvec, the parametrization of the basis
%          - scalar order, the order of the derivative required

% Outputs: - array B, the basis of the required tangent space, of size n+1 x length(tau)

function [Phi] = OB_derivative(degree, tvec, order)
    % Find number of steps (time increments)
    steps = length(tvec);

    % Initialize variable for n-order curve
    Phi = zeros(degree+1,steps); 
    
    % Calculation of the orthogonal Bernstein polynomials 
    for i = 0:degree
        % Indexing parameter
        k = i+1;

        % compute the derivatives of the Bernstein basis
        for j = 0:i 
            % Switch the derivative order
            B = bernstein_derivative(degree-j, tvec, order);

            % Compute the scaling 
            num = nchoosek(2*degree+1-j,i-j) * nchoosek(i,j);
            den = nchoosek(degree-j,i-j);
            K = num/den;

            % Compute the orthonomal basis
            Phi(k,:) = Phi(k,:) + (-1)^j * K * B(i-j+1,:);
        end

        % Final scaling
        Phi(k,:) = sqrt(2*(degree-i)+1)*Phi(k,:);
    end
end