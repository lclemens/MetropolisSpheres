%% Read Metropolis Spheres Summary Output
% Plot spheres

% Radius of sphere to read
for r=1:23

folder = '~/Documents/pub/lclemens/polymer-c_runs/20191010MetropolisSpheresTCRConfig';

M = dlmread(fullfile(folder,['MetropolisSpheres.',num2str(r)]));

%% Read File
ntTotal(r) = M(1,1);
NSphere = M(2,1);
sRadius = M(4,1);
E = M(5,1);

if(sRadius ~= r)
    disp('Wrong file');
end

for j=1:6
    rAnchor.x(j) = M(5+j,1);
    rAnchor.y(j) = M(5+j,2);
    rAnchor.z(j) = M(5+j,3);
end

for j=1:NSphere
    rSphere.x(j) = M(5+6+j,1);
    rSphere.y(j) = M(5+6+j,2);
    rSphere.z(j) = M(5+6+j,3);
end

for j=1:NSphere
    rPolymer.x(j) = M(5+6+NSphere+j,1);
    rPolymer.y(j) = M(5+6+NSphere+j,2);
    rPolymer.z(j) = M(5+6+NSphere+j,3);
end
%% Confirm constraints
passedTF=1;
contourLength = [42,29,27,39,31,27,39,31,29,42];

% Check sphere-polymer dist
spherePolymerDist = sqrt((rSphere.x(:)-rPolymer.x(:)).^2 + (rSphere.y(:)-rPolymer.y(:)).^2 + (rSphere.z(:)-rPolymer.z(:)).^2);
spherePolymerDistErr = spherePolymerDist-sRadius;

if(any(abs(spherePolymerDistErr)>0.05))
    disp('Sphere-Polymer Distance Failure');
    passedTF=0;
end

% Check polymer-anchor distances
polymerInd = [1 2 3 6 9 10]
for i=1:length(polymerInd)
    polyAnchDist = sqrt((rAnchor.x(i)-rPolymer.x(polymerInd(i))).^2 + (rAnchor.y(i)-rPolymer.y(polymerInd(i))).^2 + (rAnchor.z(i)-rPolymer.z(polymerInd(i))).^2);
    if(polyAnchDist > contourLength(polymerInd(i)))
        disp('Polymer-Anchor Distance Failure');
        disp(polymerInd(i));
        passedTF=0;
    end
end

% Check polymer-polymer distances
for i=[3 4 6 7]
    polyAnchDist = sqrt((rPolymer.x(i)-rPolymer.x(i+1)).^2 + (rPolymer.y(i)-rPolymer.y(i+1)).^2 + (rPolymer.z(i)-rPolymer.z(i+1)).^2);
    if(polyAnchDist > contourLength(i+1))
        disp('Polymer-Polymer Distance Failure');
        disp(i);
        passedTF=0;
    end
end

% Check sphere-sphere interactions
for ib=1:NSphere
    for ib2 = (ib+1):NSphere
        if( ((rSphere.x(ib)-rSphere.x(ib2)).^2 + (rSphere.y(ib)-rSphere.y(ib2)).^2 + (rSphere.z(ib)-rSphere.z(ib2)).^2) <= sRadius.^2)
            disp('Sphere-Sphere Failure');
            passedTF=0;
        end
    end
end


% Check sphere-membrane
if(any(rSphere.z(:)<sRadius))
    disp('Sphere-Membrane Failure');
    passedTF=0;
end


% Check polymer-membrane
if(any(rPolymer.z(:)<0))
    disp('Sphere-Membrane Failure');
    passedTF=0;
end

% Print passed if finished constraints
if(passedTF)
    disp('Radius');
    disp(r);
    disp('Passed constraints.');
else
    disp('Failed constraints.');
end

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
            plot3([rAnchor.x(i),rSphere.x(i)],[rAnchor.y(i),rSphere.y(i)],[rAnchor.z(i),rSphere.z(i)],'-r','LineWidth',lw);
        case 4
            plot3([rAnchor.x(i),rSphere.x(i+2)],[rAnchor.y(i),rSphere.y(i+2)],[rAnchor.z(i),rSphere.z(i+2)],'-r','LineWidth',lw);
        case {5,6}
            plot3([rAnchor.x(i),rSphere.x(i+4)],[rAnchor.y(i),rSphere.y(i+4)],[rAnchor.z(i),rSphere.z(i+4)],'-r','LineWidth',lw);
    end
end

% Ligands
for j=1:NSphere
    [ligandX,ligandY,ligandZ] = ellipsoid(rSphere.x(j),rSphere.y(j),rSphere.z(j),sRadius,sRadius,sRadius);
    ligand = surf(ligandX,ligandY,ligandZ);
    ligand.FaceAlpha = 0.2;
    ligand.FaceColor = [j*20/255 j*20/255 (255-j*20)/255]
end

% Polymer between anchors 1-3 to ligands 1-3
for i=[3 4 6 7]
    plot3([rSphere.x(i),rSphere.x(i+1)],[rSphere.y(i),rSphere.y(i+1)],[rSphere.z(i),rSphere.z(i+1)],'-r','LineWidth',lw);
end

% Membrane
[membraneX, membraneY, membraneZ] = meshgrid(-100:100,-100:100,0);
membranePlot = surf(membraneX, membraneY, membraneZ,'EdgeColor','none','LineStyle','none');
set(membranePlot,'FaceAlpha',0.9);
set(membranePlot,'FaceColor',[0.6745 0.8588 0.9647]);
membranePlot.FaceLighting = 'none';





