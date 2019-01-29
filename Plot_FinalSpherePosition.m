%% Read Metropolis Spheres Summary Output
% Plot spheres

% Radius of sphere to read
i=11;

folder = '~/Documents/pub/lclemens/polymer-c_runs/20190125MetropolisSpheres';

M = dlmread(fullfile(folder,['MetropolisSpheres.',num2str(i)]));

%% Read File

NSphere = M(2,1);
sRadius = M(3,1);

if(sRadius ~= i)
    disp('Wrong file');
end

for j=1:6
    rAnchor.x(j) = M(3+j,1);
    rAnchor.y(j) = M(3+j,2);
    rAnchor.z(j) = M(3+j,3);
end

for j=1:NSphere
    rSphere.x(j) = M(3+6+j,1);
    rSphere.y(j) = M(3+6+j,2);
    rSphere.z(j) = M(3+6+j,3);
end

%% Plot

lw = 2;
ms = 10;
figure(1); clf; hold on;

% Anchors
for j=1:6
    plot3(rAnchor.x(j),rAnchor.y(j),rAnchor.z(j),'xr','MarkerSize',ms,'LineWidth',lw);
end

% Polymer between anchors 1-3 to ligands 1-3
for i=1:6
    switch i
        case {1,2,3}
            plot3([rAnchor.x(i),rSphere.x(i)],[rAnchor.y(i),rSphere.y(i)],[rAnchor.z(i),rSphere.z(i)],'-','LineWidth',lw);
        case 4
            plot3([rAnchor.x(i),rSphere.x(i+2)],[rAnchor.y(i),rSphere.y(i+2)],[rAnchor.z(i),rSphere.z(i+2)],'-','LineWidth',lw);
        case {5,6}
            plot3([rAnchor.x(i),rSphere.x(i+4)],[rAnchor.y(i),rSphere.y(i+4)],[rAnchor.z(i),rSphere.z(i+4)],'-','LineWidth',lw);
    end
end

% Ligands
for j=1:NSphere
    [ligandX,ligandY,ligandZ] = ellipsoid(rSphere.x(j),rSphere.y(j),rSphere.z(j),sRadius,sRadius,sRadius);
    surf(ligandX,ligandY,ligandZ);
end

% Polymer between anchors 1-3 to ligands 1-3
for i=[3 4 6 7]
    plot3([rSphere.x(i),rSphere.x(i+1)],[rSphere.y(i),rSphere.y(i+1)],[rSphere.z(i),rSphere.z(i+1)],'-','LineWidth',lw);
end

% Membrane
[membraneX, membraneY, membraneZ] = meshgrid(-100:100,-100:100,0);
membranePlot = surf(membraneX, membraneY, membraneZ,'EdgeColor','none','LineStyle','none');
set(membranePlot,'FaceAlpha',0.9);
set(membranePlot,'FaceColor',[0.6745 0.8588 0.9647]);
membranePlot.FaceLighting = 'none';





