/*** Allard Group jun.allard@uci.edu                    ***/

void outputSummary();
void dataRecording();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void outputSummary()
{
    if (!verboseTF)
    {
        fList = fopen(listName, "a");
        
        fprintf(fList, "%ld\n%d\n%f\n%f\n",
                nt,        // 1
                NSphere,   // 2
                sRadius,    // 3
                E);         // 4
        
        for(i=0;i<6;i++)
        {
            fprintf(fList, "%f %f %f\n",rAnchor[i][0],rAnchor[i][1],rAnchor[i][2]);
        }
        
        for(i=0;i<NSphere;i++)
        {
            fprintf(fList, "%f %f %f\n", rSphere[i][0],rSphere[i][1],rSphere[i][2]);
        }
        
        
        fclose(fList);
    }

}

void dataRecording()
{
    if (verboseTF)
    {
        
        if ( nt>=0 && (nt <= 200000) )
        {
            // output results to file
            fList = fopen(listName, "a");
            
            fprintf(fList, "%ld %d %f %f %ld ",
                    nt,                         // 1
                    NSphere,                    // 2
                    sRadius,                    // 3
                    E,                          // 4
                    constraintProposalsTotal);  // 5
            
            for(i=0;i<6;i++)
            {
                fprintf(fList, "%f %f %f ",rAnchor[i][0],rAnchor[i][1],rAnchor[i][2]);
            }
            
            
            for(i=0;i<NSphere;i++)
            {
                fprintf(fList, "%f %f %f ",rSphere[i][0],rSphere[i][1],rSphere[i][2]);
            }
            
            fprintf(fList, "\n");
            fclose(fList);
        }
    } // finished verbose output
    
    
}



