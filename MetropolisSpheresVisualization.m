%% Visualize Spheres
clear all;

M = dlmread('output');

%% Pull out anchors and spheres

ntTotal = M(end,1)+1;
NSphere = M(1,2);
sRadius = M(1,3);


for i=1:6
    rAnchor.x(:,i) = M(:,5+3*(i-1));
    rAnchor.y(:,i) = M(:,6+3*(i-1));
    rAnchor.z(:,i) = M(:,7+3*(i-1));
end
for i=1:NSphere
    rSphere.x(:,i) = M(:,5+6*3+0+3*(i-1));
    rSphere.y(:,i) = M(:,5+6*3+1+3*(i-1));
    rSphere.z(:,i) = M(:,5+6*3+2+3*(i-1));
    i
    M(1,5+6*3+0+3*(i-1))
    M(1,5+6*3+1+3*(i-1))
    M(1,5+6*3+2+3*(i-1))
end

%% Plot Spheres moving around

lw = 2;

for nt=1:1000:ntTotal
    figure(1); clf; hold on; view([-48,7]);
    for i=1:6
        plot3(rAnchor.x(nt,i),rAnchor.y(nt,i),rAnchor.z(nt,i),'xr','LineWidth',lw);
    end
    for i=1:NSphere
        [sphereX, sphereY, sphereZ] = ellipsoid(rSphere.x(nt,i),rSphere.y(nt,i),rSphere.z(nt,i),sRadius,sRadius,sRadius);
        sphere = surf(sphereX,sphereY,sphereZ);
    end
    xlim([-100,100]);
    ylim([-100,100]);
end


