%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial trajectory guess given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - scalar tfapp, the approximated time of flight initial guess
%         - vector x0, the initial 6 by 1 state vector (position and velocity)
%         - vector xf, the final 6 by 1 state vector (position and velocity)
%         - scalar thetaf, the final optimal insertion phase
%         - cell array B, the polynomial basis to be used
%         - string basis, to select the polynomial basis to be used in the
%           approximation
%         - vector n, the order of the approximating polynomial function 
%         - array P0, the initial array of control points 

% Outputs: - array P, the updated boundary conditions control points, of dimensions size(x,1) x n+1 

function [P] = boundary_conditions(Problem, beta, t0, tf, B, basis, n, P0)
    % Boundary conditions
    [s0, sf] = feval(@(initial, final, beta, t0, tf)Problem.BoundaryConditions(initial, final, beta, t0, tf), Problem.initial, Problem.final, beta, t0, tf);

    % Constants 
    m = Problem.StateDim;       % State dimension
    L = Problem.DerDeg;         % Order of the dynamics (maximum derivative order)

    if (size(s0,1) ~= L * m)
        error('Cauchy boundary conditions cannot be imposed. Dimensions of the final boundary conditions are not compatible with the dynamics.');
    end

    if (size(sf,1) ~= L * m)
        error('Cauchy boundary conditions cannot be imposed. Dimensions of the final boundary conditions are not compatible with the dynamics.');
    end

    P = P0;         % Initialization

    % Dimensionalizing of the generalized velocities 
    if (L > 1)
        for i = 1:L-1
            s0(1+m*i:m*(i+1),1) = (tf - t0)^i * s0(1+m*i:m*(i+1),1);
            sf(1+m*i:m*(i+1),1) = (tf - t0)^i * sf(1+m*i:m*(i+1),1);
        end
    end

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = s0(1:m);
            P(:,2) = s0(1:m) + (L > 1) * s0(m+1:end)./n;
            for i = 1:length(n)
                P(i,n(i)) = sf(i) - (L > 1) * sf(m+i)/n(i);
                P(i,n(i)+1) = sf(i);
            end

        case 'Chebyshev'
            % Symbolic regression 
            for i = 1:size(P,1)
                l = n(i)-2;
                X0(1,1) = s0(i)-sum(P0(i,3:n(i)-1).*(-1).^(2:l),2); 
                X0(2,1) = sf(i)-sum(P0(i,3:n(i)-1),2);

                if (L > 1)
                    X0(3,1) = s0(size(P,1)+i)-sum((-1).^(1:(l-1)).*P0(i,3:n(i)-1).*(2:l).^2,2);
                    X0(4,1) = sf(size(P,1)+i)-sum(P0(i,3:n(i)-1).*(2:l).^2,2);
                    A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); ...
                         1 1 1 1; ... 
                         0 1 (-1)^(n(i)-2)*(n(i)-1)^2 (-1)^(n(i)-1)*n(i)^2; ...
                         0 1 (n(i)-1)^2 n(i)^2];
                else
                    A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); 1 1 1 1]; 
                end

                sol = A\X0;
                P(i,[1 2 n(i) n(i)+1]) = sol;
            end

        case 'Legendre'
            % Symbolic regression 
            for i = 1:size(P,1)
                l = n(i)-2;
                X0(1,1) = s0(i)-sum(P0(i,3:n(i)-1).*(-1).^(2:l),2); 
                X0(2,1) = sf(i)-sum(P0(i,3:n(i)-1),2);

                if (L > 1)
                    X0(3,1) = s0(size(P,1)+i)-sum((-1).^(1:(l-1)).*P0(i,3:n(i)-1).*(2:l).*(3:l+1)/2,2);
                    X0(4,1) = sf(size(P,1)+i)-sum(P0(i,3:n(i)-1).*(2:l).*(3:l+1)/2,2);
                    A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); ...
                         1 1 1 1; ... 
                         0 1 (-1)^(n(i)-2)*(n(i)-1)*n(i)/2 (-1)^(n(i)-1)*n(i)*(n(i)+1)/2; ...
                         0 1 (n(i)-1)*n(i)/2 n(i)*(n(i)+1)/2];
                else
                    A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); 1 1 1 1];
                end

                sol = A\X0;
                P(i,[1 2 n(i) n(i)+1]) = sol;
            end

        otherwise
            C = zeros(L*m, 2);     % Compute the partial state evolution                

            for i = 1:length(n)
                C(i,:) = P0(i,3:n(i)-1)*B{i}(3:n(i)-1,[1 end]);

                if (L > 1)
                    C(length(n)+i,:) = P0(i,3:n(i)-1)*B{i}(n(i)+4:2*n(i),[1 end]);
                end
            end

            % Assemble the linear system 
            C = [s0-C(:,1) sf-C(:,end)];

            if (L > 1)
                C = [C(1:size(P,1),:) C(size(P,1)+1:2*size(P,1),:)];
            else
                C = C(1:size(P,1),:);
            end
            
            % Control points 
            for i = 1:length(n)
                index = [1 2 n(i) n(i)+1 n(i)+2 n(i)+3 2*n(i)+1 2*(n(i)+1)];
                A = B{i}(index,[1,end]);

                if (L > 1)
                    A = [A(1:4,:) A(5:8,:)];
                else
                    A = A(1:4,:);
                end

                P(i,[1 2 n(i) n(i)+1]) = C(i,:)*A^(-1);
            end
     end
end