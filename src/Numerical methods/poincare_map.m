%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: poincare_map.m 
% Issue: 0 
% Validated: 30/05/21

%% Poincare map %%
% This script contains the function to compute several Poincare surfaces of
% section

% Inputs: - structure Su, a portion of the unstable manifold
%         - structure SS, a portion of the stable manifold
%         - string map, selecting the Poincare map to compute
% Output: - array S0, the initial conditions on the manifolds
%         - vector TOF, containing the time of flight on each manifold

% New versions: 

function [S0, TOF] = poincare_map(Su, Ss, map)
    %Constants 
    tol = 1e-3;     
    
    %Branch the different Poincare maps
    switch (map)
        case 'Homoclinic rendezvous'
            %Preallocation of the intersection points 
            Iu = zeros(size(Su.Trajectory,1), size(Su.Trajectory,3));
            Is = zeros(size(Ss.Trajectory,1), size(Ss.Trajectory,3));
            
            %Extract the intersection points 
            for i = 1:size(Iu,1)
                %Intersection of the stable manifold with the Poincare section
                Iu(i,:) = shiftdim(Su.Trajectory(i,Su.ArcLength(i),:));    
            end
            
            for i = 1:size(Is,1)
                %Intersection of the stable manifold with the Poincare section
                if (size(Is,1) == 1)
                    Is(1,:) = shiftdim(Ss.Trajectory(1,Ss.ArcLength,:));
                else
                    Is(i,:) = shiftdim(Ss.Trajectory(i,Ss.ArcLength(i),:));   
                end
            end
            
            %Select the connection point in the configuration space       
            figure
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1), Iu(i,4), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1), Is(i,4), 'go');
            end
            hold off
            grid on; 
            title('Poincare map in the x tangent space')
            xlabel('Synodic x coordinate')
            ylabel('Synodic $V_x$ coordinate')
            
            figure 
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,3), Iu(i,6), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,3), Is(i,6), 'go');
            end
            hold off
            grid on; 
            title('Poincare map in the z tangent space')
            xlabel('Synodic z coordinate')
            ylabel('Synodic $V_z$ coordinate')
           
            figure 
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1), Iu(i,3), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1), Is(i,3), 'go');
            end
            hold off
            grid on; 
            title('Poincare map in the XZ plane')
            xlabel('Synodic x coordinate')
            ylabel('Synodic z coordinate')
            [x0, z0] = getpts(); 
            
            %Configuration space initital conditions
            X0 = [x0 0 z0];                         
               
            %Selected initial conditions on the unstable manifold
            k = 1; 
            Ius = zeros(size(Iu));
            indexu = zeros(size(Iu,1),2);
            for i = 1:size(Iu,1) 
                if (norm(Iu(i,1:3)-X0(1:3)) < tol)
                    Ius(k,:) = Iu(i,:);  
                    indexu(k) = i;
                    k = k+1;
                end
            end
            
            %Save the intersection results
            if (k == 1)
                error('The unstable manifold does not intersect the Poincare section at the selected point');
            else
                Ius = Ius(1:k-1,:);         %Selected intersection point
                indexu = indexu(1:k-1);     %Selected intersection point
            end
            
            %Selected initial conditions on the stable manifold
            k = 1;
            Iss = zeros(size(Is));
            indexs = zeros(size(Is,1),1);
            for i = 1:size(Is,1)
                flag = norm(Is(i,1:3)-X0(1:3)) < tol;     
                if (flag)
                    Iss(k,:) = Is(i,:);
                    indexs(k) = i;
                    k = k+1;
                end
            end
            
            %Save the intersection results
            if (k == 1)
                error('The stable manifold does not intersect the Poincare section at the selected point');
            else
                Iss = Iss(1:k-1,:);         %Selected intersection point
                indexs = indexs(1:k-1);     %Selected intersection point
            end

            %Determine the size of the impulsive maneuver
            dV = zeros(size(Ius,1), size(Iss,1)); 
            for i = 1:size(Ius,1)
                for j = 1:size(Iss,1)
                    dV(i,j) = norm(Ius(i,4:6)-Iss(j,4:6));    %Needed velocity change
                end
            end
          
            %Find the minimum maneuver size
            [Vr, Ir] = min(dV); 
            [~, Ic] = min(Vr);
            
            %Output
            S0.Unstable = shiftdim(Su.Trajectory(indexu(Ir(Ic)),1,:));
            S0.Stable = shiftdim(Ss.Trajectory(indexs(Ic),1,:));
            TOF(1) = Su.TOF(indexu(Ir(Ic))); 
            TOF(2) = Ss.TOF(indexs(Ic));

        case 'X crossing'
            %Preallocation of the intersection points 
            Iu = zeros(size(Su.Trajectory,1), size(Su.Trajectory,3));
            Is = zeros(size(Ss.Trajectory,1), size(Ss.Trajectory,3));
            
            %Extract the intersection points 
            for i = 1:size(Iu,1)
                %Intersection of the stable manifold with the Poincare section
                Iu(i,:) = shiftdim(Su.Trajectory(i,Su.ArcLength(i),:));    
            end
            
            for i = 1:size(Is,1)
                %Intersection of the stable manifold with the Poincare section
                if (size(Is,1) == 1)
                    Is(1,:) = shiftdim(Ss.Trajectory(1,Ss.ArcLength,:));
                else
                    Is(i,:) = shiftdim(Ss.Trajectory(i,Ss.ArcLength(i),:));   
                end
            end
            
            %Select the connection point in the configuration space       
            figure
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1), Iu(i,4), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1), Is(i,4), 'ko');
            end
            hold off
            grid on; 
            title('Poincare map in the x tangent space')
            xlabel('Synodic x coordinate')
            ylabel('Synodic $V_x$ coordinate')
            
            figure
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1), Iu(i,6), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1), Is(i,6), 'ko');
            end
            hold off
            grid on; 
            title('Poincare map in the z tangent space')
            xlabel('Synodic z coordinate')
            ylabel('Synodic $V_z$ coordinate')
           
            figure 
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1), Iu(i,3), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1), Is(i,3), 'go');
            end
            hold off
            grid on; 
            title('Poincare map in the XZ plane')
            xlabel('Synodic x coordinate')
            ylabel('Synodic z coordinate')
            [x0, z0] = getpts(); 
            
            %Configuration space initital conditions
            X0 = [x0 0 z0];                         
               
            %Selected initial conditions on the unstable manifold
            k = 1; 
            Ius = zeros(size(Iu));
            indexu = zeros(size(Iu,1),2);
            for i = 1:size(Iu,1) 
                if (norm(Iu(i,1:3)-X0(1:3)) < tol)
                    Ius(k,:) = Iu(i,:);  
                    indexu(k) = i;
                    k = k+1;
                end
            end
            
            %Save the intersection results
            if (k == 1)
                error('The unstable manifold does not intersect the Poincare section at the selected point');
            else
                Ius = Ius(1:k-1,:);         %Selected intersection point
                indexu = indexu(1:k-1);     %Selected intersection point
            end
            
            %Selected initial conditions on the stable manifold
            k = 1;
            Iss = zeros(size(Is));
            indexs = zeros(size(Is,1),1);
            for i = 1:size(Is,1)
                flag = norm(Is(i,1:3)-X0(1:3)) < tol;     
                if (flag)
                    Iss(k,:) = Is(i,:);
                    indexs(k) = i;
                    k = k+1;
                end
            end
            
            %Save the intersection results
            %Save the intersection results
            if (k == 1)
                error('The stable manifold does not intersect the Poincare section at the selected point');
            else
                Iss = Iss(1:k-1,:);         %Selected intersection point
                indexs = indexs(1:k-1);     %Selected intersection point
            end
            
            %Determine the size of the impulsive maneuver
            dV = zeros(size(Ius,1), size(Iss,1)); 
            for i = 1:size(Ius,1)
                for j = 1:size(Iss,1)
                    dV(i,j) = norm(Ius(i,4:6)-Iss(j,4:6));    %Needed velocity change
                end
            end
          
            %Find the minimum maneuver size
            [Vr, Ir] = min(dV); 
            [~, Ic] = min(Vr);
            
            %Output
            S0.Unstable = shiftdim(Su.Trajectory(indexu(Ir(Ic)),1,:));
            S0.Stable = shiftdim(Ss.Trajectory(indexs(Ic),1,:));
            TOF(1) = Su.TOF(indexu(Ir(Ic))); 
            TOF(2) = Ss.TOF(indexs(Ic));
            
        otherwise 
            error('No valid map was selected');
    end
end