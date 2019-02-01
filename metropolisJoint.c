/*** Allard Group jun.allard@uci.edu                    ***/

void metropolisJoint();


int energyPasses;
void metropolisJoint()
{
    energyPasses = 0;

    /************************************************************************************/
    /********************** INITIALIZE CONFIGURATION *****************************/
    /************************************************************************************/


    // Structure of rSphere:
    // Sphere 1: x,y,z
    // ...
    // Sphere 10: x,y,z
    // Polymer 1: x,y,z
    // ...
    // Polymer 10: x,y,z
    
    /************************************************************************************/
    // Establish sphere anchors in membrane
    int baseSepDistance = 5;
    for(i=0;i<6;i++)
    {
        rAnchor[i][0] = sqrt(baseSepDistance*baseSepDistance + 2.5*2.5) * cos( floor((double)i/2)*(2*PI/3) + pow(-1,(double)(i+1))*atan( 2.5/baseSepDistance) );
        rAnchor[i][1] = sqrt(baseSepDistance*baseSepDistance + 2.5*2.5) * sin( floor((double)i/2)*(2*PI/3) + pow(-1,(double)(i+1))*atan( 2.5/baseSepDistance) );
        rAnchor[i][2] = 0;
    }
    
    /************************************************************************************/
    //initalize spheres in vertical position oriented over their respective anchors
    
    // spheres 0 to 2 attached to anchors 0 to 2
    for(i=0;i<3;i++)
    {
        rSphere[i][0] = rAnchor[i][0];
        rSphere[i][1] = rAnchor[i][1];
        rSphere[i][2] = sRadius;
    }
    // spheres 3,4 attached to sphere 2
    for(i=3;i<5;i++)
    {
        rSphere[i][0] = rAnchor[2][0];
        rSphere[i][1] = rAnchor[2][1];
        rSphere[i][2] = sRadius*(i-1);
    }
    // sphere 5 attached to anchor 3
    for(i=5;i<6;i++)
    {
        rSphere[i][0] = rAnchor[3][0];
        rSphere[i][1] = rAnchor[3][1];
        rSphere[i][2] = sRadius;
    }
    // spheres 6,7 attached to sphere 5
    for(i=6;i<8;i++)
    {
        rSphere[i][0] = rAnchor[3][0];
        rSphere[i][1] = rAnchor[3][1];
        rSphere[i][2] = sRadius*(i-1);
    }
    // spheres 8, 9 attached to anchors 4, 5
    for(i=8;i<10;i++)
    {
        rSphere[i][0] = rAnchor[i-4][0];
        rSphere[i][1] = rAnchor[i-4][1];
        rSphere[i][2] = sRadius;
    }
    
    /************************************************************************************/
    // initialize polymer points
    for(i=0;i<3;i++)
    {
        rSphere[i+10][0] = rAnchor[i][0];
        rSphere[i+10][1] = rAnchor[i][1]-sRadius;
        rSphere[i+10][2] = sRadius;
    }
    // spheres 3,4 attached to sphere 2
    for(i=3;i<5;i++)
    {
        rSphere[i+10][0] = rAnchor[2][0];
        rSphere[i+10][1] = rAnchor[2][1]-sRadius;
        rSphere[i+10][2] = sRadius*(i-1);
    }
    // sphere 5 attached to anchor 3
    for(i=5;i<6;i++)
    {
        rSphere[i+10][0] = rAnchor[3][0];
        rSphere[i+10][1] = rAnchor[3][1]-sRadius;
        rSphere[i+10][2] = sRadius;
    }
    // spheres 6,7 attached to sphere 5
    for(i=6;i<8;i++)
    {
        rSphere[i+10][0] = rAnchor[3][0];
        rSphere[i+10][1] = rAnchor[3][1]-sRadius;
        rSphere[i+10][2] = sRadius*(i-1);
    }
    // spheres 8, 9 attached to anchors 4, 5
    for(i=8;i<10;i++)
    {
        rSphere[i+10][0] = rAnchor[i-4][0];
        rSphere[i+10][1] = rAnchor[i-4][1]-sRadius;
        rSphere[i+10][2] = sRadius;
    }
    
    /************************************************************************************/
    /********************** INITIALIZE LOOP PARAMETERS *****************************/
    /************************************************************************************/
    
    int failedMembrane=0;
    int failedSpheres=0;
    int failedContourLength=0;
    int failedRadiusCheck=0;
	endConstraintPassedTF=0;
	nt=0;
    E = INF;

    //THINK ABOUT THIS PART HERE!!!!!!!!!!!!!
    // Do we need a different dChi for each filament? or is one ok?
    // stuff needed for automatic perturbation size adaptation
	for(iParam=0;iParam<2;iParam++)
	{
		dChi[iParam] = DCHIINIT;
		proposals[iParam] = 0;
		accepts[iParam] = 0;
	}


    /********* BEGIN LOOP THROUGH ITERATIONS! *******************/
	while(!endConstraintPassedTF && nt < NTMAX) // Time loop!
    {

        // adapt step size
        if (!(nt % NTADAPT))
        {
            //THINK ABOUT THIS PART HERE!!!!!!!!!!!!!
            rate[0] = (double)accepts[0]/(double)proposals[0];
            //printf("accepts[0]: %ld\n",accepts[0]);
            //printf("proposals[0]: %ld\n",proposals[0]);
            //printf("rate[0]: %ld\n",rate[0]);
            rate[1] = (double)accepts[1]/(double)proposals[1];
            
            //printf("current accepts/proposals: %d/%d=%f, %d/%d=%f\n", accepts[0], proposals[0], rate[0], accepts[1], proposals[1], rate[1]);
            
            for(iParam=0;iParam<2;iParam++)
            {
                accepts[iParam] = 0;
                proposals[iParam] = 0;
                if (nt < NTCHECK)
                {
                    if (rate[iParam] > 0.5 || rate[iParam] < 0.3)
                    {
                        dChi[iParam] = dChi[iParam]*rate[iParam]/0.44;
                        if (dChi[iParam] > 1)
                            dChi[iParam] = 1;
                        if (dChi[iParam] < DCHIMIN)
                            dChi[iParam] = DCHIMIN;
                    }
                }
            }
            
        } // done adapting step

        /********* OUTLINE OF ALGORITHM *******************/
        // 1. We create a new proposal configuration and
        // then decide:
        // 2. whether the proposal configuration satisfies any constraints (e.g., floor?), and
        // 3. whether the proposal configuration is accepted by Metropolis test (always accept for a freely-jointed chain)
        // After the configuration update, we:
        // 4. collect any data of interest (e.g., ree; whether there is steric occlusion), and write it to file
        // 5. test for stationarity (convergence)
        // 6. increment the iteration and repeat
		
        /****************************************************************/
        /********* 1. Generate proposal configuration *******************/
        /****************************************************************/
		
		iPropose = floor(2*NSphere*TWISTER); // Initialize. This is the sphere we will adjust this time.
        
        //debugging
        if(0)
        {
            printf("Sphere: %ld\n",iPropose);
            fflush(stdout);
        }
        
        
        //THINK ABOUT THIS PART HERE!!!!!!!!!!!!!
		if(iPropose==0)
		{
			dChiHere = dChi[0];
			proposals[0]++;
		}
		else 
		{
			dChiHere = dChi[1];
			proposals[1]++;
		}
		
        constraintProposalsTotal = 0;
        long rounds = 0;
		constraintSatisfiedTF = 0;

		while(!constraintSatisfiedTF && constraintProposalsTotal < CPMAX) // keep generating proposal configurations until we get one that satisfies the constraint
        {
            // initialize proposal configuration
            for(i=0;i<2*NSphere;i++)
            {
                rSpherePropose[i][0] = rSphere[i][0];
                rSpherePropose[i][1] = rSphere[i][1];
                rSpherePropose[i][2] = rSphere[i][2];
            }
            
            iChiPropose = floor(3*TWISTER);
            
            //debugging
            if(0)
            {
                printf("iChiPropose: %ld\n",iChiPropose);
                fflush(stdout);
            }
            // propose new configuration angles
            // We can use normal (Gaussian) random variables using the Box-Muller method (see eg wikipedia)
            // or uniform random variable

            //phiPropose[nfPropose][iPropose]   = phiPropose[nfPropose][iPropose]   + dChiHere*sqrt(-2*log(TWISTER))*cos(2*PI*TWISTER);
            rSpherePropose[iPropose][iChiPropose]   = rSpherePropose[iPropose][iChiPropose]   + dChiHere*(2*TWISTER-1);
           
            /****************************************************************/
            /******************* 2. Test constraints ************************/
            /****************************************************************/
            
            constraintSatisfiedTF=1;
            
            // No hard constraints right now

            
            constraintProposalsTotal++; //count number of times proposals are rejected this time step
        } //finish constraint while loop
        
        if (constraintProposalsTotal >= CPMAX) //if number of proposals exceeds CPMAX, exit program
        {
            printf("Exceeded maximum proposals.\n nt: %ld\n",nt);
            fflush(stdout);
            
            exit(0);
        }
        
        /****************************************************************/
        /********************* 3. Metropolis test ***********************/
        /****************************************************************/
        // We now have a proposal configuration that passes the constraints.
        // Step 3 is to see if it passes our acceptance test (Metropolis test).
        
            // Compute energy
            ENew = 0;
        
            /****************************************************************/
            // Force attaching spheres to anchors
        
            /* Structure of rSphere:
               Sphere 1: x,y,z
               ...
               Sphere 10: x,y,z
               Polymer 1: x,y,z
               ...
               Polymer 10: x,y,z
               Sphere set up
             */
            // Location of spheres, polymers
            /*
             1 - first sphere, only on anchor
             2 - second sphere, only on anchor
             3 - third sphere, on anchor and other sphere
             4 - fourth sphere, between sphere 3 and 5
             5 - fifth sphere, attached to sphere 4
             6 - sixth sphere, on anchor and 7th sphere
             7 - seventh sphere, between sphere 6 and 8
             8 - eighth sphere, attached to sphere 7
             9 - ninth sphere, only on anchor
             10 - tenth sphere, only on anchor
             */
            // For indexing polymer points instead of spheres, add 10
        
        
            // add energy for distance between polymer points and sphere points
            for(i=0;i<NSphere;i++)
            {
                distCurrent = sqrt((rSpherePropose[i+10][0]-rSpherePropose[i][0])*(rSpherePropose[i+10][0]-rSpherePropose[i][0])+
                                         (rSpherePropose[i+10][1]-rSpherePropose[i][1])*(rSpherePropose[i+10][1]-rSpherePropose[i][1])+
                                         (rSpherePropose[i+10][2]-rSpherePropose[i][2])*(rSpherePropose[i+10][2]-rSpherePropose[i][2]));
                
                ENew += 0.5*kBound*(distCurrent-sRadius)*(distCurrent-sRadius);
            
            }
        
            // add energy for polymers attached to anchor points if polymer point is further away than allowed contour length
            for(i=0;i<3;i++)
            {
                    distCurrent = sqrt((rSpherePropose[i+10][0]-rAnchor[i][0])*(rSpherePropose[i+10][0]-rAnchor[i][0])+
                                             (rSpherePropose[i+10][1]-rAnchor[i][1])*(rSpherePropose[i+10][1]-rAnchor[i][1])+
                                            (rSpherePropose[i+10][2]-rAnchor[i][2])*(rSpherePropose[i+10][2]-rAnchor[i][2]));
                
                if(i==0)
                    contourLength = 42;
                if(i==1)
                    contourLength = 29;
                if(i==2)
                    contourLength = 27;
                
                if(distCurrent>contourLength)
                {
                    ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
                }
            }
            i=5;
            distCurrent = sqrt((rSpherePropose[i+10][0]-rAnchor[3][0])*(rSpherePropose[i+10][0]-rAnchor[3][0])+
                                 (rSpherePropose[i+10][1]-rAnchor[3][1])*(rSpherePropose[i+10][1]-rAnchor[3][1])+
                                 (rSpherePropose[i+10][2]-rAnchor[3][2])*(rSpherePropose[i+10][2]-rAnchor[3][2]));

                contourLength = 27;
        
                if(distCurrent>contourLength)
                {
                    ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
                }
        
            for(i=8;i<10;i++)
            {
                distCurrent = sqrt((rSpherePropose[i+10][0]-rAnchor[i-4][0])*(rSpherePropose[i+10][0]-rAnchor[i-4][0])+
                                         (rSpherePropose[i+10][1]-rAnchor[i-4][1])*(rSpherePropose[i+10][1]-rAnchor[i-4][1])+
                                         (rSpherePropose[i+10][2]-rAnchor[i-4][2])*(rSpherePropose[i+10][2]-rAnchor[i-4][2]));

                if(i==8)
                    contourLength = 29;
                if(i==9)
                    contourLength = 42;
                
                if(distCurrent>contourLength)
                {
                    ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
                }
            }
        
            // add energy for polymers attached to other polymers if polymer point is further away than allowed contour length
            // polymer 4 attached to polymer point 3
            i1=3+10;
            i2=2+10;
            distCurrent = sqrt((rSpherePropose[i1][0]-rSpherePropose[i2][0])*(rSpherePropose[i1][0]-rSpherePropose[i2][0])+
                                     (rSpherePropose[i1][1]-rSpherePropose[i2][1])*(rSpherePropose[i1][1]-rSpherePropose[i2][1])+
                                     (rSpherePropose[i1][2]-rSpherePropose[i2][2])*(rSpherePropose[i1][2]-rSpherePropose[i2][2]));
        
            contourLength = 39;
    
            if(distCurrent>contourLength)
            {
                ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
            }
        
            //polymer point 5 attached to polymer point 4
            i1=4+10;
            i2=3+10;
            distCurrent = sqrt((rSpherePropose[i1][0]-rSpherePropose[i2][0])*(rSpherePropose[i1][0]-rSpherePropose[i2][0])+
                                     (rSpherePropose[i1][1]-rSpherePropose[i2][1])*(rSpherePropose[i1][1]-rSpherePropose[i2][1])+
                                     (rSpherePropose[i1][2]-rSpherePropose[i2][2])*(rSpherePropose[i1][2]-rSpherePropose[i2][2]));
        
            contourLength = 31;
        
            if(distCurrent>contourLength)
            {
                ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
            }
            //polymer point 7 attached to polymer point 6
            i1=6+10;
            i2=5+10;
            distCurrent = sqrt((rSpherePropose[i1][0]-rSpherePropose[i2][0])*(rSpherePropose[i1][0]-rSpherePropose[i2][0])+
                                     (rSpherePropose[i1][1]-rSpherePropose[i2][1])*(rSpherePropose[i1][1]-rSpherePropose[i2][1])+
                                     (rSpherePropose[i1][2]-rSpherePropose[i2][2])*(rSpherePropose[i1][2]-rSpherePropose[i2][2]));
            contourLength = 39;
        
            if(distCurrent>contourLength)
            {
                ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
            }
        
            //polymer point 8 attached to polymer point 7
            i1=7+10;
            i2=6+10;
            distCurrent = sqrt((rSpherePropose[i1][0]-rSpherePropose[i2][0])*(rSpherePropose[i1][0]-rSpherePropose[i2][0])+
                                     (rSpherePropose[i1][1]-rSpherePropose[i2][1])*(rSpherePropose[i1][1]-rSpherePropose[i2][1])+
                                     (rSpherePropose[i1][2]-rSpherePropose[i2][2])*(rSpherePropose[i1][2]-rSpherePropose[i2][2]));
        
            contourLength = 31;
        
            if(distCurrent>contourLength)
            {
                ENew += 0.5*kBound*(distCurrent-contourLength)*(distCurrent-contourLength);
            }
            
            /****************************************************************/
            // Bound ligand energies

                /**********************/
                // energy between bound ligands and membrane
                for (i=0;i<NSphere;i++) //for each bound ligand on filament
                {
                    
                    if (MEMBRANE)
                    {
                        if(rSpherePropose[i][2]<sRadius) // if any bound ligands intersect with membrane
                        {
                            // add energy based on intersection distance
                            ENew += 0.5*kBound*(rSpherePropose[i][2]-sRadius)*(rSpherePropose[i][2]-sRadius);
                        }
                    }
                }
                /**********************/
                // energy between polymer points and membrane
                for (i=NSphere;i<2*NSphere;i++) //for each polymer point
                {
                    
                    if (MEMBRANE)
                    {
                        if(rSpherePropose[i][2]<0) // if any bound ligands intersect with membrane
                        {
                            // add energy based on intersection distance
                            ENew += 0.5*kBound*(rSpherePropose[i][2])*(rSpherePropose[i][2]);
                        }
                    }
                }

                /**********************/
                // energy between bound ligands and other bound ligands
                for (i=0;i<NSphere;i++)
                {
                    for (ib2=(i+1);ib2<NSphere;ib2++) //for each next ligand
                    {
                        
                        if ((rSpherePropose[i][0]-rSpherePropose[ib2][0])*(rSpherePropose[i][0]-rSpherePropose[ib2][0]) +
                            (rSpherePropose[i][1]-rSpherePropose[ib2][1])*(rSpherePropose[i][1]-rSpherePropose[ib2][1]) +
                            (rSpherePropose[i][2]-rSpherePropose[ib2][2])*(rSpherePropose[i][2]-rSpherePropose[ib2][2])<=
                            (2*sRadius)*(2*sRadius)) //if distance between centers is less than 2*brLigand, then ligands are intersecting
                        {
                            boundCentertoBoundDistance = sqrt(((rSpherePropose[i][0]-rSpherePropose[ib2][0])*(rSpherePropose[i][0]-rSpherePropose[ib2][0]) +
                                                              (rSpherePropose[i][1]-rSpherePropose[ib2][1])*(rSpherePropose[i][1]-rSpherePropose[ib2][1]) +
                                                              (rSpherePropose[i][2]-rSpherePropose[ib2][2])*(rSpherePropose[i][2]-rSpherePropose[ib2][2]))) - (2*sRadius);
                            
                            //ENew += 0.5*kBound*boundCentertoBoundDistance*boundCentertoBoundDistance;
                            
                            //Soft spring to allow spheres to pass through each other if necessary
                            ENew += 0.5*kSphere*boundCentertoBoundDistance*boundCentertoBoundDistance;

                        }
                    }
                }
            /****************************************************************/
            // Accept or reject based on energy
            // should this be <=? Do we reject normal probability for no force?
            if (  TWISTER < exp(E-ENew) ) //always accepts if ENew<E, accepts with normal (?) probability if ENew>E
            {
                E = ENew;
                
                if(E<0.0001)
                {
                    printf("E: %f\n",E);
                    fflush(stdout);

                   for(i=0;i<2*NSphere;i++)
                   {
                       printf("Sphere Position %ld: %f, %f, %f\n",i,rSphere[i][0],rSphere[i][1],rSphere[i][2]);
                   }
                }
                if(0)
                {
                    printf("E: %f\n",E);
                }
                
                // Make configuration into the proposal configuration
                for(i=0;i<2*NSphere;i++)
                {
                    rSphere[i][0] = rSpherePropose[i][0];
                    rSphere[i][1] = rSpherePropose[i][1];
                    rSphere[i][2] = rSpherePropose[i][2];
                 
                }
                
                if(iPropose==0)
                    accepts[0] ++;
                else 		
                    accepts[1] ++;
                
                if(0)
                {
                    for(i=0;i<2*NSphere;i++)
                    {
                        printf("Sphere Position %ld: %f, %f, %f\n",i,rSphere[i][0],rSphere[i][1],rSphere[i][2]);
                    }
                }
                energyPasses++;
                
            }
         // end configuration loop
        
        /**************************************************************************/
        /**************** 4. Data collection and output to file *******************/
        /**************************************************************************/
        
        // check if spheres are intersecting and if they are past FJC limits
        if (1)
        {
            
            endConstraintPassedTF = 1;
            
            // check if spheres are intersecting membrane
            if (MEMBRANE)
            {
                for(i=0;i<NSphere;i++)
                {
                    if(rSphere[i][2]<sRadius) // if any bound ligands intersect with membrane
                    {
                        if(0)
                        {
                            printf("Membrane\n");
                            printf("%f\n",rSphere[i][2]);
                        }
                        //if(nt>=100000)
                            failedMembrane++;
                        
                        endConstraintPassedTF = 0;
                        i=NSphere;
                    }
                }
                if(endConstraintPassedTF)
                {
                    for(i=NSphere;i<2*NSphere;i++)
                    {
                        if(rSphere[i][2]<0) // if any bound ligands intersect with membrane
                        {
                            if(0)
                            {
                                printf("Membrane\n");
                                printf("%f\n",rSphere[i][2]);
                            }
                            //if(nt>=100000)
                            failedMembrane++;
                            
                            endConstraintPassedTF = 0;
                            i=2*NSphere;
                        }
                    }
                }
            }
            
            // check if spheres are intersecting other spheres
            if(endConstraintPassedTF)
            {
                for (i=0;i<NSphere;i++) //for each bound ligand on filament
                {
                    for (ib2=(i+1);ib2<NSphere;ib2++) //for each next ligand
                    {
                        
                        if ((rSphere[i][0]-rSphere[ib2][0])*(rSphere[i][0]-rSphere[ib2][0]) +
                            (rSphere[i][1]-rSphere[ib2][1])*(rSphere[i][1]-rSphere[ib2][1]) +
                            (rSphere[i][2]-rSphere[ib2][2])*(rSphere[i][2]-rSphere[ib2][2])<=
                            (2*sRadius)*(2*sRadius)) //if distance between centers is less than 2*brLigand, then ligands are intersecting
                        {
                            if(0)
                            {
                                printf("Spheres\n");
                            }
                            //if(nt>=100000)
                                failedSpheres++;
                            
                            endConstraintPassedTF = 0;
                            ib2=NSphere;
                            i=NSphere;
                        }
                     }
                  }
            }
            // check if polymers are further apart than allowed by inextensible WLC
            if(endConstraintPassedTF)
            {
                for(i=0;i<3;i++)
                {
                    if(i==0)
                        contourLength = 42;
                    if(i==1)
                        contourLength = 29;
                    if(i==2)
                        contourLength = 27;
                    if(sqrt((rSphere[i+10][0]-rAnchor[i][0])*(rSphere[i+10][0]-rAnchor[i][0])+
                                             (rSphere[i+10][1]-rAnchor[i][1])*(rSphere[i+10][1]-rAnchor[i][1])+
                                             (rSphere[i+10][2]-rAnchor[i][2])*(rSphere[i+10][2]-rAnchor[i][2])) > contourLength)
                    {
                        if(0)
                        {
                            printf("WLC\n");
                            printf("%ld,%f,%f\n",i,contourLength,sqrt((rSphere[i+10][0]-rAnchor[i][0])*(rSphere[i+10][0]-rAnchor[i][0])+
                                                                      (rSphere[i+10][1]-rAnchor[i][1])*(rSphere[i+10][1]-rAnchor[i][1])+
                                                                      (rSphere[i+10][2]-rAnchor[i][2])*(rSphere[i+10][2]-rAnchor[i][2])));
                        }
                        //if(nt>=100000)
                            failedContourLength++;
                        
                        endConstraintPassedTF = 0;
                        i=3;
                    }
                    
                }
            }
            if(endConstraintPassedTF)
            {
                i=5;
                contourLength = 27;
                if(sqrt((rSphere[i+10][0]-rAnchor[3][0])*(rSphere[i+10][0]-rAnchor[3][0])+
                   (rSphere[i+10][1]-rAnchor[3][1])*(rSphere[i+10][1]-rAnchor[3][1])+
                   (rSphere[i+10][2]-rAnchor[3][2])*(rSphere[i+10][2]-rAnchor[3][2])) > contourLength)
                {
                    if(0)
                    {
                        printf("WLC\n");
                        printf("%ld,%f,%f\n",i,contourLength,sqrt((rSphere[i+10][0]-rAnchor[3][0])*(rSphere[i+10][0]-rAnchor[3][0])+
                                                                  (rSphere[i+10][1]-rAnchor[3][1])*(rSphere[i+10][1]-rAnchor[3][1])+
                                                                  (rSphere[i+10][2]-rAnchor[3][2])*(rSphere[i+10][2]-rAnchor[3][2])));
                    }
                    //if(nt>=100000)
                        failedContourLength++;
                    
                    endConstraintPassedTF = 0;
                }
            }
            if(endConstraintPassedTF)
            {
                for(i=8;i<10;i++)
                {
                    if(i==8)
                        contourLength = 29;
                    if(i==9)
                        contourLength = 42;
                    
                    if(sqrt((rSphere[i+10][0]-rAnchor[i-4][0])*(rSphere[i+10][0]-rAnchor[i-4][0])+
                       (rSphere[i+10][1]-rAnchor[i-4][1])*(rSphere[i+10][1]-rAnchor[i-4][1])+
                       (rSphere[i+10][2]-rAnchor[i-4][2])*(rSphere[i+10][2]-rAnchor[i-4][2])) > contourLength)
                    {
                        if(0)
                        {
                            printf("WLC\n");
                            printf("%ld,%f,%f\n",i,contourLength,sqrt((rSphere[i+10][0]-rAnchor[i-4][0])*(rSphere[i+10][0]-rAnchor[i-4][0])+
                               (rSphere[i+10][1]-rAnchor[i-4][1])*(rSphere[i+10][1]-rAnchor[i-4][1])+
                               (rSphere[i+10][2]-rAnchor[i-4][2])*(rSphere[i+10][2]-rAnchor[i-4][2])));
                        }
                        //if(nt>=100000)
                            failedContourLength++;
                        
                        endConstraintPassedTF = 0;
                        i=10;
                    }
                    
                }
            }
            
            // test if polymer points further apart than allowed by chain
            // polymer 4 attached to polymer 3
            if(endConstraintPassedTF)
            {
                i1=3+10;
                i2=2+10;
                
                contourLength = 39;
                if((sqrt(rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                                         (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                   (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])) > contourLength)
                {
                    if(0)
                    {
                        printf("WLC\n");
                        printf("%d,%d,%f,%f\n",i1,i2,contourLength,sqrt((rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                                                                        (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                                                                        (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])));
                    }
                    //if(nt>=100000)
                        failedContourLength++;
                    
                    endConstraintPassedTF = 0;
                }

            }
            // polymer 5 attached to polymer 4
            if(endConstraintPassedTF)
            {
                i1=4+10;
                i2=3+10;
                contourLength = 31;
                if((sqrt(rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                   (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                   (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])) > contourLength)
                {
                    if(0)
                    {
                        printf("WLC\n");
                        printf("%d,%d,%f,%f\n",i1,i2,contourLength,sqrt((rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                           (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                                                                        (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])));
                    }
                    //if(nt>=100000)
                        failedContourLength++;
                    
                    endConstraintPassedTF = 0;
                }
            }
            // polymer 7 attached to polymer 6
            if(endConstraintPassedTF)
            {
                i1=6+10;
                i2=5+10;
                contourLength = 39;
                if((sqrt(rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                   (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                   (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])) > contourLength)
                {
                    if(0)
                    {
                        printf("WLC\n");
                        printf("%d,%d,%f,%f\n",i1,i2,contourLength,sqrt((rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                                                                        (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                                                                        (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])));
                    }
                    //if(nt>=100000)
                        failedContourLength++;
                    
                    endConstraintPassedTF = 0;
                }
            }
            
            // polymer 8 attached to polymer 7
            if(endConstraintPassedTF)
            {
                i1=7+10;
                i2=6+10;
                contourLength = 31;
                if((sqrt(rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                   (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                   (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])) > contourLength)
                {
                    if(0)
                    {
                        printf("WLC\n");
                        printf("%d,%d,%f,%f\n",i1,i2,contourLength,sqrt((rSphere[i1][0]-rSphere[i2][0])*(rSphere[i1][0]-rSphere[i2][0])+
                                                                        (rSphere[i1][1]-rSphere[i2][1])*(rSphere[i1][1]-rSphere[i2][1])+
                                                                        (rSphere[i1][2]-rSphere[i2][2])*(rSphere[i1][2]-rSphere[i2][2])));
                    }
                    //if(nt>=100000)
                        failedContourLength++;
                    
                    endConstraintPassedTF = 0;
                }
            }
                               
            // test how far polymers are from sphere centers
           if(endConstraintPassedTF)
           {
               for(i=0;i<NSphere;i++)
               {
                   if(((sqrt(rSphere[i+10][0]-rSphere[i][0])*(rSphere[i+10][0]-rSphere[i][0])+
                       (rSphere[i+10][1]-rSphere[i][1])*(rSphere[i+10][1]-rSphere[i][1])+
                       (rSphere[i+10][2]-rSphere[i][2])*(rSphere[i+10][2]-rSphere[i][2]))-sRadius) > 0.1)
                   {
                       //if(nt>=100000)
                       failedRadiusCheck++;
                       
                       endConstraintPassedTF = 0;
                   }
               }
           }
                               
                               
            if(0)
            {
                for(i=0;i<2*NSphere;i++)
                {
                    printf("Sphere %ld: %f, %f, %f\n",i,rSphere[i][0],rSphere[i][1],rSphere[i][2]);
                }
            }
            
            
            
        } //
        
        /********* 6. Increment time *******************/
        dataRecording();
		nt++;
		
	} // done time loop
    
    //debugging
    printf("Energy Passes: %f\n",(float)energyPasses/(float)nt);
    fflush(stdout);
    
    printf("Energy Fails: %f\n",1-((float)energyPasses/(float)nt));
    fflush(stdout);
    
    printf("Failed Membrane: %f\n",(float)failedMembrane/(float)nt);
    fflush(stdout);
    
    printf("Failed Spheres: %f\n",(float)failedSpheres/(float)nt);
    fflush(stdout);
    
    printf("Failed Contour Length: %f\n",(float)failedContourLength/(float)nt);
    fflush(stdout);
    
    printf("Failed Contour Length: %f\n",(float)failedRadiusCheck/(float)nt);
    fflush(stdout);
    
    printf("E Final: %f\n",E);
    fflush(stdout);
    
    if(endConstraintPassedTF)
    {
       printf("Not sterically constrained.\n");
        for(i=0;i<2*NSphere;i++)
        {
            printf("Sphere Position %ld: %f, %f, %f\n",i,rSphere[i][0],rSphere[i][1],rSphere[i][2]);
        }
       outputSummary();
     }
                                                                      
    if(nt>=NTMAX)
       printf("Past max iterations.\n");
    

}

