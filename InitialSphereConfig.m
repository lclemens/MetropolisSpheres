%% Initial Configuration 
sRadius = 2;
rSphere = zeros(20,3);
rAnchor = zeros(6,3);
baseSepDistance = 5;
    for i=0:5
        rAnchor(i+1,1) = sqrt(baseSepDistance*baseSepDistance + 2.5*2.5) * cos( floor(i/2)*(2*pi/3) + (-1)^(i+1)*atan( 2.5/baseSepDistance) );
        rAnchor(i+1,2) = sqrt(baseSepDistance*baseSepDistance + 2.5*2.5) * sin( floor(i/2)*(2*pi/3) + (-1)^(i+1)*atan( 2.5/baseSepDistance) );
        rAnchor(i+1,3) = 0;
    end
    

    %initalize spheres in vertical position oriented over their respective anchors
    
    % spheres 0 to 2 attached to anchors 0 to 2
    for i=0:2
        rSphere(i+1,1) = rAnchor(i+1,1);
        rSphere(i+1,2) = rAnchor(i+1,2);
        rSphere(i+1,3) = sRadius;
    end
    % spheres 3,4 attached to sphere 2
    for i=3:4
        rSphere(i+1,1) = rAnchor(3,1);
        rSphere(i+1,2) = rAnchor(3,2);
        rSphere(i+1,3) = 2*sRadius*(i-1)+i*0.0001;
    end
    % sphere 5 attached to anchor 3
    for i=5
        rSphere(i+1,1) = rAnchor(4,1);
        rSphere(i+1,2) = rAnchor(4,2);
        rSphere(i+1,3) = sRadius;
    end
    % spheres 6,7 attached to sphere 5
    for i=6:7
        rSphere(i+1,1) = rAnchor(4,1);
        rSphere(i+1,2) = rAnchor(4,2);
        rSphere(i+1,3) = 2*sRadius*(i-4)+i*0.0001;
    end
    
    % spheres 8, 9 attached to anchors 4, 5
    for i=8:9
        rSphere(i+1,1) = rAnchor(i-4+1,1);
        rSphere(i+1,2) = rAnchor(i-4+1,2);
        rSphere(i+1,3) = sRadius;
    end
        %% Polymer points
    % spheres 0 to 2 attached to anchors 0 to 2
    for i=0:2
        rSphere(i+11,1) = rAnchor(i+1,1);
        rSphere(i+11,2) = rAnchor(i+1,2)-sRadius;
        rSphere(i+11,3) = sRadius;
    end
    % spheres 3,4 attached to sphere 2
    for i=3:4
        rSphere(i+11,1) = rAnchor(3,1);
        rSphere(i+11,2) = rAnchor(3,2)-sRadius;
        rSphere(i+11,3) = 2*sRadius*(i-1)+i*0.0001;
    end
    % sphere 5 attached to anchor 3
    for i=5
        rSphere(i+11,1) = rAnchor(4,1);
        rSphere(i+11,2) = rAnchor(4,2)-sRadius;
        rSphere(i+11,3) = sRadius;
    end
    % spheres 6,7 attached to sphere 5
    for i=6:7
        rSphere(i+11,1) = rAnchor(4,1);
        rSphere(i+11,2) = rAnchor(4,2)-sRadius;
        rSphere(i+11,3) = 2*sRadius*(i-4)+i*0.0001;
    end
    % spheres 8, 9 attached to anchors 4, 5
    for i=8:9
        rSphere(i+11,1) = rAnchor(i-4+1,1);
        rSphere(i+11,2) = rAnchor(i-4+1,2)-sRadius;
        rSphere(i+11,3) = sRadius;
    end
    
    %% Plot
    lw = 10;
    figure(1); clf; hold on; box on;
    for i=1:6
        plot3(rAnchor(i,1),rAnchor(i,2),rAnchor(i,3),'xk','LineWidth',lw);
    end
    for i=11:20
        plot3(rSphere(i,1),rSphere(i,2),rSphere(i,3),'*r','LineWidth',lw);
    end
    for i=1:10
        [elliX, elliY, elliZ] = ellipsoid(rSphere(i,1),rSphere(i,2),rSphere(i,3),sRadius,sRadius,sRadius);
        ligand = surf(elliX,elliY,elliZ);
        ligand.FaceAlpha = 0.3;
        ligand.FaceColor = [i*20/255 i*20/255 i*20/255]
    end
    
    
    
    
    
    
    