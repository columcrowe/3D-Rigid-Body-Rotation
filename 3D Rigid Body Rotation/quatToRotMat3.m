function R = quatToRotMat3(q1)
%function algebraically manipulates a quaternion rotation p'=qpq^-1 into a matrix rotation p' = Rp
%Quaternion-derived rotation matrix - https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
%Input: quaternion q represented by [w,xi,yj,zk]
%Output: 3x3 rotation matrix R 

%$Date:7/02/2016 $Colum Crowe $Revision:1

w = q1(1);
x = q1(2);
y = q1(3);
z = q1(4);

R = [1 - 2*(y^2 + z^2), 2*(x*y - z*w), 2*(x*z + y*w) ; 2*(x*y + z*w), 1 - 2*(x^2 + z^2), 2*(y*z - x*w ); 2*(x*z - y*w ), 2*(y*z + x*w ), 1 - 2 *(x^2 + y^2)];

end