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
        
        for(i=0;i<NSphere;i++)
        {
            fprintf(fList, "%f %f %f\n", rSphere[i][0],rSphere[i][1],rSphere[i][2]);
        }

        fclose(fList);
    }

}

void dataRecording()
{
    
    
    
}



