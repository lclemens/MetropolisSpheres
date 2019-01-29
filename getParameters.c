/*** Allard Group jun.allard@uci.edu                    ***/

void getParameters();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/
char tmpString[100];

/********************************************************************************************************/
void getParameters()
{
    paramsFile = fopen(paramsFilename,"r");
    
    fscanf(paramsFile,"%s %s", tmpString, listName);
    if (TALKATIVE) printf("This is output file name: %s\n", listName);
    
    fscanf(paramsFile,"%s %d", tmpString, &NSphere);
    if (TALKATIVE) printf("This is number of spheres: %d\n", NSphere);
    
    fscanf(paramsFile,"%s %lf", tmpString, &sRadius);
    if (TALKATIVE) printf("This is sphere radius: %lf\n", sRadius);
    
    fscanf(paramsFile,"%s %lf", tmpString, &kBound);
    if (TALKATIVE) printf("This is bound ligand spring constant: %lf\n", kBound);
    
    fscanf(paramsFile,"%s %lf", tmpString, &kSphere);
    if (TALKATIVE) printf("This is sphere spring constant: %lf\n", kSphere);
    
    fscanf(paramsFile,"%s %d", tmpString, &verboseTF);
    if (TALKATIVE) printf("This is sphere spring constant: %d\n", verboseTF);
    
    fclose(paramsFile);
    
}

/********************************************************************************************************/

