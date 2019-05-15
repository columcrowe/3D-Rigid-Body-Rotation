function [Qd] = quatDiff(wb, q1)
%function to perform quaternion differentiation
%https://fgiesen.wordpress.com/2012/08/24/quaternion-differentiation/
%formula is for constant or near constant angular velocity therefore
%a relatively small time step (dt) is required so that equation is valid

%$Date:20/02/2016 $Colum Crowe $Revision:2

%convert quaternion representing orientation in world frame to 3x3 rotation matrix
R = quatToRotMat3(q1);

%rotating angular velocity into the world frame
ww = R*wb;

%quaternion representing angular veolcity in the world frame
q2=[0; ww(1); ww(2); ww(3)];

%quaternion multiplication of the the angular velocity in the world frame
%and the current orientation of the rigid body gives the change in orientation
%of the rigid body in the world reference frame

Qd = 0.5*quatMult(q2,q1); % = 0.5*quatMult(q1,[0,wb(1),wb(2),wb(3)]); %http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf

end