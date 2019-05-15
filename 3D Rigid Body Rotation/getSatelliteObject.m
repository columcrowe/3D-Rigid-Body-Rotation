function [M] = getSatelliteObject(height,length,width,timestep,tf,Moment_of_Inertia,Angular_Velocity)
%Get frames to make a movie in order to visualise rotation of tumbling
%satellite under zero external torque but with changing angular velocity.

%defining variables
h=height; %x-axis
l=length; %y-axis
w=width; %z-axis

%axis size
axis_size=max(w,max(h,l))+10;

%%
%visualization

p1=[-h,-l,-w];
p2=[h,-l,-w];
p3=[h,l,-w];
p4=[-h,l,-w];
p5=[-h,-l,w];
p6=[h,-l,w];
p7=[h,l,w];
p8=[-h,l,w];

vertices_bottom=[p1;p2;p3;p4];
vertices_top=[p5;p6;p7;p8];
vertices = [vertices_bottom; vertices_top];

faces = [[1:4];[1:4]' [[2:4] 1]' [(4+[2:4]) 5]' ((4+[1:4])')*ones(1,1);4+[1:4]];

%%

figure;
hold on;
grid off
axis([-axis_size axis_size -axis_size axis_size -axis_size axis_size]);
axis off
view(3);

% Compare two methods - Runge-Kutta performs better
[Storage]=getSatellitePositionsRK4(timestep,tf,Moment_of_Inertia,Angular_Velocity);
% [Storage]=getSatellitePositionsABM(timestep,tf,Moment_of_Inertia,Angular_Velocity);

%%

N=size(Storage,1);

for count=1:N
    
    R=quatToRotMat3(Storage(count,5:8));
    vertices = [vertices_bottom; vertices_top];
    vertices = (R * vertices')';
    S1 = patch('Vertices',vertices, 'Faces',faces, 'FaceColor','r');

    if count>=2 && ( Storage(count,13)-Storage(count-1,13)>=1e-6 || Storage(count,14)-Storage(count-1,14)>=1e-6 || Storage(count,15)-Storage(count-1,15)>=1e-6 )
        h(1)=line([0,Storage(count,13)*axis_size],[0,Storage(count,14)*axis_size],[0,Storage(count,15)*axis_size],'color', 'r',...        %angular momentum H
    'LineStyle','--','linewidth', 1.5);
    else
        h(1)=line([0,Storage(count,13)*axis_size],[0,Storage(count,14)*axis_size],[0,Storage(count,15)*axis_size],'color', 'k',...        %angular momentum H
    'LineStyle','--','linewidth', 1.5);
    end

    h(2)=line([0,Storage(count,10)*axis_size],[0,Storage(count,11)*axis_size],[0,Storage(count,12)*axis_size],'color','b',...   %angular velocity Ww    
    'LineStyle',':','linewidth', 1);
 
%     h(3)=line([0,Storage(count,2)*axis_size],[0,Storage(count,3)*axis_size],[0,Storage(count,4)*axis_size],'color','c',...   %angular velocity Wb    
%     'LineStyle',':','linewidth', 1);
    
    M(count) = getframe();
    cla;
    count;
end

close

end