%% Read Metropolis Spheres Summary Output
% Plot spheres

savefolder = '/Volumes/GoogleDrive/My Drive/Papers/MultisiteDisorder/Data_Figures';
savesubfolder = '3.SimultaneousBinding/MetropolisSpheres';
saveTF = 1;
% Radius of sphere to read
for r=25

folder = '~/Documents/pub/lclemens/polymer-c_runs/20191010MetropolisSpheresTCRConfig';

M = dlmread(fullfile(folder,['MetropolisSpheres.',num2str(r)]));

%% Read File
ntTotal(r) = M(1,1);
NSphere = M(2,1);
sRadius = 0.3*M(4,1); % nm
E = M(5,1);

if(sRadius ~= 0.3*r)
    disp('Wrong file');
end

for j=1:6
    rAnchor.x(j) = 0.3*M(5+j,1); % nm
    rAnchor.y(j) = 0.3*M(5+j,2); % nm
    rAnchor.z(j) = 0.3*M(5+j,3); % nm
end

for j=1:NSphere
    rSphere.x(j) = 0.3*M(5+6+j,1); % nm
    rSphere.y(j) = 0.3*M(5+6+j,2); % nm
    rSphere.z(j) = 0.3*M(5+6+j,3); % nm
end

for j=1:NSphere
    rPolymer.x(j) = 0.3*M(5+6+NSphere+j,1); % nm
    rPolymer.y(j) = 0.3*M(5+6+NSphere+j,2); % nm
    rPolymer.z(j) = 0.3*M(5+6+NSphere+j,3); % nm
end
%% Confirm constraints
passedTF=1;
contourLength = 0.3*[42,29,27,39,31,27,39,31,29,42];  % nm

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
    disp(r); disp('Kuhn lengths');
    disp('Passed constraints.');
else
    disp('Failed constraints.');
end

end
%% Plot

lw = 3;
ms = 15;
colors_fil = [0.7 0 0; 0 0.5 0.8; 0 0.5 0; 0 0.8 0; 0.7 0 0.7; 1 0 0]
figure(1); clf; hold on;

% Anchors
for j=1:6
    pM = plot3(rAnchor.x(j),rAnchor.y(j),rAnchor.z(j),'x','MarkerSize',ms,'LineWidth',lw);
    pM.Color = colors_fil(j,:);
end

% Polymer between anchors 1-6 to ligands 1-6
for i=1:6
    switch i
        case {1,2,3}
            pL = plot3([rAnchor.x(i),rSphere.x(i)],[rAnchor.y(i),rSphere.y(i)],[rAnchor.z(i),rSphere.z(i)],'-','LineWidth',lw);
            pL.Color = colors_fil(i,:);
        case 4
            pL = plot3([rAnchor.x(i),rSphere.x(i+2)],[rAnchor.y(i),rSphere.y(i+2)],[rAnchor.z(i),rSphere.z(i+2)],'-','LineWidth',lw);
            pL.Color = colors_fil(i,:);
        case {5,6}
            pL = plot3([rAnchor.x(i),rSphere.x(i+4)],[rAnchor.y(i),rSphere.y(i+4)],[rAnchor.z(i),rSphere.z(i+4)],'-','LineWidth',lw);
            pL.Color = colors_fil(i,:);
    end
end

% Ligands
for j=1:NSphere
    [ligandX,ligandY,ligandZ] = ellipsoid(rSphere.x(j),rSphere.y(j),rSphere.z(j),sRadius,sRadius,sRadius);
    ligand = surf(ligandX,ligandY,ligandZ);
    colormap parula
    ligand.FaceAlpha = 0.4;
    %ligand.FaceColor = [j*20/255 j*20/255 (255-j*20)/255];
    ligand.EdgeAlpha = 0.1;
    %ligand.FaceColor = [(255-j*20)/255 (255-j*20)/255 (255-j*20)/255]; %b&w
    %ligand.FaceColor = [0 0 0];
end

% Polymer between anchors 1-3 to ligands 1-3
for i=[3 4 6 7]
    pL = plot3([rSphere.x(i),rSphere.x(i+1)],[rSphere.y(i),rSphere.y(i+1)],[rSphere.z(i),rSphere.z(i+1)],'-r','LineWidth',lw);
    if(i==3 || i==4)
        pL.Color = colors_fil(3,:);
    else
        pL.Color = colors_fil(4,:);
    end
end

% Membrane
[membraneX, membraneY, membraneZ] = meshgrid(0.3.*(-100:100),0.3.*(-100:150),0);
membranePlot = surf(membraneX, membraneY, membraneZ,'EdgeColor','none','LineStyle','none');
set(membranePlot,'FaceAlpha',0.9);
set(membranePlot,'FaceColor',[0.6745 0.8588 0.9647]); % blue membrane
set(membranePlot,'FaceColor',[1 0.98 0.747]); % yellow membrane
%colors_membrane = jet(20); 
%set(membranePlot,'FaceColor',colors_membrane(1,:)); % dark blue jet
%membrane
membranePlot.FaceLighting = 'none';

set(gca,'Color',[0.9 0.9 0.9]);

% Labels
xlabel('x (nm)','FontName','Arial','FontSize',18);
ylabel('y (nm)','FontName','Arial','FontSize',18);
zlabel('z (nm)','FontName','Arial','FontSize',18);

% set view
set(gca,'view',[131.2000 62.8000]);
set(gca,'view',[128.0000 60.4000]);
set(gca,'view',[-32.8000 62.0000]);
set(gca,'view',[106.4000 59.6000]);

% all bases are not intersecting

set(gca,'view',[106.8000 86.8000]);

set(gca,'view',[86.4000 64.4000]);

set(gca,'view',[-81.2000 66.0000]);

% possible good paired views, but still not sufficient to show all pairs
% don't intersect
set(gca,'view',[255.2000   69.2000]);
set(gca,'view',[165.2000   78.0000]);


set(gca,'view',[119.2000   76.4000]);
set(gca,'view',[-98.0000   70.8000]);
  
set(gca,'view',[66.0000   70.0000]);
if(saveTF)
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View1.fig'),'fig');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View1.pdf'),'pdf');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View1.eps'),'epsc');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View1.png'),'png');
end
set(gca,'view',[196.4000   60.4000]);
set(gca,'view',[-69.6000   62.8000]);

if(saveTF)
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View2.fig'),'fig');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View2.pdf'),'pdf');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View2.eps'),'epsc');
    saveas(gcf,fullfile(savefolder,savesubfolder,'MetSpheres25_View2.png'),'png');
end
