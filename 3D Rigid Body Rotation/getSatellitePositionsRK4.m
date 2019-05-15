function [Storage] = getSatellitePositionsRK4(timestep,tf,Moments_of_Inertia,Angular_Velocity)
%Gets a quaternion representation for the orientation of a rigid body with
%respect to the world reference frame in order to model the 3D rotation of
%a satellite of given dimensions spinning in orbit under negligible external torque. 
%Function Call:
%getSatellitePosition(timestep,tf,[Ixx Iyy Izz],[Wx Wy Wz])
%Example: 
%[Storage]=getSatellitePositionsRK4(0.05,100,[9361 200000 300000],[0.1 -0.01 0.05])
%Inputs:
%[Ixx,Iyy,Izz] are the moments of inertia in kg.m^2 about each respective axis
%[Wx,Wy,Wz] are the angular velocities in rad/s about the principal axes at time zero
%timestep is the increment in time (seconds) between each numerical integration step
%tf is the final time in seconds at which the position of the satellite will be
%determined by a 4th order Runge-Kutta method (RK4).
%Outputs:
%Storage array containing changing angular velocity and orientation of tumbling
%satellite in quaternion form at time tf and constant angular momentum.

%$Date:27/04/2016 $Colum Crowe $Revision:5

%%

Tol=0.5; %tolerance for numerical error

%defining variables
Ixx = Moments_of_Inertia(1);
Iyy = Moments_of_Inertia(2);
Izz = Moments_of_Inertia(3);
Ib = [Ixx 0 0;0 Iyy 0; 0 0 Izz]; %creating inertia matrix in body plane

Wx = Angular_Velocity(1);
Wy = Angular_Velocity(2);
Wz = Angular_Velocity(3);
Wb=[Wx; Wy; Wz];    %angular velocity in the body frame

%initialising the orientation as the unit quaternion in the world frame
Q=[1;0;0;0];

Qmag=Q(1)^2+Q(2)^2+Q(3)^2+Q(4)^2;   %magnitude of unit quaternion should equal 1

t = 0;              %time t0
dt = timestep;      %step size
N = tf*(1/dt);      %no. of iterations

R = quatToRotMat3( Q );         %convert q to 3x3 rotation matrix
    
Ww=R*Wb;                    %determine angular velocity in world frame
Hw=R*(Ib*Wb);               %determine angular momentum in world frame

%creating storage array to hold variables
Storage = zeros(N+1,15);
Storage(1,1) = t;
Storage(1,2:4) = Wb;
Storage(1,5:8) = Q;
Storage(1,9) = Qmag;
Storage(1,10:12) = Ww;
Storage(1,13:15) = Hw;

%%

for count=1:N

%--------------------------------------------------------------------------
%Fourth Order Runge-Kutta Algorithm
%Integrating angular velocity in the body frame and thus quaternion
%representing orientation of the rigid body in the world frame

    kw1 = dt * solveEulers(Wb,Ib);
    kq1 = dt * quatDiff(Wb,Q);
    z1 = Wb + 0.5*kw1;
    z2 = Q + 0.5*kq1;
    
    kw2 = dt * solveEulers( z1, Ib );
    kq2 = dt * quatDiff( z1, z2 );
    z1 = Wb + 0.5*kw2;
    z2 = Q + 0.5*kq2;
    
    kw3 = dt * solveEulers( z1, Ib );
    kq3 = dt * quatDiff( z1, z2 );
    z1 = Wb + kw3;
    z2 = Q+kq3;
    
    kw4 = dt * solveEulers( z1, Ib );
    kq4 = dt * quatDiff( z1, z2 );
    
    %new angular velocity in the body reference frame
    Wb = Wb + (1/6) * (kw1+(2*kw2)+(2*kw3)+kw4);
    
    %quaternion representating new orientation of the rigid body in world frame
    Q = Q + (1/6) * (kq1+(2*kq2)+(2*kq3)+kq4);
    
    %----------------------------------------------------------------------
    %check if quaternion is still of unit length and if not renormalize
    %to prevent problem of quaternion drift
    if Qmag~=1
        
        Qmag=Q(1)^2+Q(2)^2+Q(3)^2+Q(4)^2;
        
        Q = Q ./((Q(1)^2+Q(2)^2+Q(3)^2+Q(4)^2)^(0.5)); 
        
    end
    %----------------------------------------------------------------------
    
    R = quatToRotMat3( Q );     %convert q to 3x3 rotation matrix
    
    Ww=R*Wb;                    %determine angular velocity in world plane
    Hw=R*(Ib*Wb);               %determine angular momentum in world plane
    
    t = t + dt ;    %increment by time step
    

    %update storage array
    Storage(count+1,1:15) = [t Wb(1) Wb(2) Wb(3) Q(1) Q(2) Q(3) Q(4) Qmag...
                               Ww(1) Ww(2) Ww(3) Hw(1) Hw(2) Hw(3)]; 
    
    %----------------------------------------------------------------------
    %numerical error check
    %angular momentum must be conserved
    
    if abs(Storage(count+1,13)-Storage(1,13))>Tol %if H not constant
        error('Warning: Numerical error. Inaccurate results. Please use a smaller timestep') %error
    elseif abs(Storage(count+1,14)-Storage(1,14))>Tol
        error('Warning: Numerical error. Inaccurate results. Please use a smaller timestep')
    elseif abs(Storage(count+1,15)-Storage(1,15))>Tol
        error('Warning: Numerical error. Inaccurate results. Please use a smaller timestep')
    end
    %----------------------------------------------------------------------

end


end


