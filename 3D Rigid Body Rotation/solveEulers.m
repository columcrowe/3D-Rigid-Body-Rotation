function [Wd] = solveEulers( W, I )
%solves for the angular velocity in the body frame
%Output: Wdot
%Inputs:
%Angular Velocity W in the body reference frame
%Inertia Tensor

%$Date:20/02/2016 $Colum Crowe $Revision:2

%Components of inertia matrix about principal axes
Ixx=I(1,1);
Iyy=I(2,2);
Izz=I(3,3);

%Angular velocity about each axis
Wx=W(1);
Wy=W(2);
Wz=W(3);

%Euler's equations (with no applied external torque) https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
Wdx=((Iyy-Izz)/Ixx)*Wy*Wz;
Wdy=((Izz-Ixx)/Iyy)*Wz*Wx;
Wdz=((Ixx-Iyy)/Izz)*Wx*Wy;

%Derivative of angular velocity in the body frame
Wd=[Wdx; Wdy; Wdz];

end