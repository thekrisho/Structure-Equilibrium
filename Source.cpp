

#include "Support.h"

//---------------------------------------------------------------------------------------
int main()
{
   JOINT *Jf;              // stores joint data from file
   MEMBER *Mf;             // stores member data from file
   FORCE *Ef, *Rf;         // stores external forces, reactions from file
   char buff[MAXLINE];     // temp char buffer
   TRUSS_SIZES ts;         // contains number of joints\members\external forces\reactions
   FILE *f;                // data file stream pointer
   double **M, **Minv;     // M is system matrix, Minv is inverse of M
   double *E, *S;          // E is external force vector (equation 10), S is solution
   int N;                  // NxN is size of M.  indx is used internally in the inverter
   int i, r, c;            // temp variables.  r=row, c=column
   char strUnits[MAXLINE]; // stores units as a string

   printf("Enter the name of the data file: ");
   gets(buff);  // like scanf. request file name 
   f=fopen(buff, "r"); // read text file and f is the file pointer to the actual file
   if(f==NULL) // if the file is empty print error 
   {
      printf("Cannot open %s.  Press ENTER to end program...\n", buff);
      getchar();
      return 0;
   }

   // get the number of joints, members, external forces, and reactions
   // so can allocate memory for truss data structures
   ts=getTrussSizes(f); // give the function the file data 

   N=2*ts.NJ;
   if(N!=ts.NM+ts.NR)
   {
      printf("2*NJ!=NM+NR.  The truss is unsolvable.  Press ENTER to end program...\n");
      getchar();
      fclose(f);
      return 0;
   }

   Jf  = new JOINT[ts.NJ];
   Mf  = new MEMBER[ts.NM];
   Rf  = new FORCE[ts.NR];
   Ef = new FORCE[ts.NEF];

   M   = new double *[N];
   Minv  = new double *[N];
   E    = new double[N];
   S    = new double[N];

   for(i=0; i<N; i++)
   {
      M[i]  = new double[N];
      Minv[i] = new double[N];
   }

   // get all the truss data from file
   if(!getTrussData(f, ts, Jf, Mf, Ef, Rf, strUnits)) // fields  
   {
      printf("Truss data not read properly.\n");
      cleanup(f, Jf, Mf, Rf, Ef, M, Minv, E, S, N);
      return 0;
   }

   //********** develop the system matrix and external force vector here*************************************************************************
  
   // SOLVE LENGTH OF MEMBERS (Using equation from handout)
   for (i = 0; i < ts.NM; i++)   
   {
      Mf[i].L = sqrtf(pow((Jf[Mf[i].j[1]].p[0]) - (Jf[Mf[i].j[0]].p[0]), 2) + pow((Jf[Mf[i].j[1]].p[1]) - (Jf[Mf[i].j[0]].p[1]),2));
   }


   // CREATE M MATRIX

   // Fill entire matrix with zero (Defines matrix size)
   for (r = 0; r < (ts.NJ * 2); r++)      
   {
	   for (c = 0; c < (ts.NJ * 2); c++)  
		   M[r][c] =  0;
   }

   // Check members for joint connectivity (Identifies joint - member - joing connectevity and places a "1 flag" in matrix position)
   for (i = 0; i < ts.NM; i++) 
   {  // X Row connections
	   M[(Mf[i].j[0]) * 2][i] = 1;
	   M[(Mf[i].j[1]) * 2][i] = 1;
      // Y Row connections 
	   M[((Mf[i].j[0]) *2 ) + 1] [i] = 1;
	   M[((Mf[i].j[1]) *2 ) + 1] [i] = 1; 
      // Reaction force connections
      if (i < ts.NR)
      {
         M[((Rf[i].j) * 2)][i + ts.NM] = 1;
         M[((Rf[i].j) * 2) + 1][i + ts.NM] = 1;
      }
   }

   int axis = 0; // Y  and  X  coordinate flage
 
   // Place values in matrix (Checks for "1 flag" and proceeds with calculation and value placement only where needed)
   printf("\n\nM Matrix\n--------\n");
   for (r = 0; r < ts.NJ*2; r++)  // NJ * 2 for (x and y)  
   {
	   for (c = 0; c < ts.NM + ts.NR; c++)  
	   {
		   if (M[r][c] == 1)       // Check for "1 flag"
		   { 
            if (r % 2)           // Check for X or Y Row and set flag
              axis = 1;

            // solve and print angles (cos and sin part of the matrix)
            if (c >= ts.NM)  
            {
               if (axis == 1)    // Y
                  M[r][c] = sin( (Rf[c - ts.NM].theta) );

               else              // X
                  M[r][c] = cos( (Rf[c - ts.NM].theta) );
            }

            // solve and print normal values (forces in x and y)
            else   
            {
               if ((Mf[c].j[0]) == (r/2))       // seen in first col, therfore use col 1 - col 0 equation
                  {
                     M[r][c] = ((Jf[ Mf[c].j[1] ].p[axis]) - (Jf[Mf[c].j[0]].p[axis])) / (Mf[c].L);
                  }
  
               else if ((Mf[c].j[1]) == (r/2))  // seen in second col, therfore use col 0 - col 1 equation
                  {
                     M[r][c] = ((Jf[Mf[c].j[0]].p[axis]) - (Jf[Mf[c].j[1]].p[axis])) / (Mf[c].L);
                  }
            }
            axis = 0;  // Reset axis flag
		   }
         printf("%s%.2f ", (M[r][c] >= 0 ? "+":"") , M[r][c]);  // Print Value
	   }
      printf("\n");
   }
  

   // PRODUCE ARRAY OF EXTERNAL FORCES

   // Fill array with zero (define array size)
   for (r = 0; r < (ts.NJ); r++)      
   {
	   for (c = 0; c < (ts.NEF); c++)  
	   {
      E[r * 2] = 0;
      E[r * 2 + 1] = 0;
      }
   }

   // Place external forces 
   for (r = 0; r < (ts.NJ); r++)      
   {
	   for (c = 0; c < (ts.NEF); c++) 
	   {
		   if ((Ef[ c ].j) == r)    // (Checks if row (joint) has an external force)
		   {
			   if (Ef[c].idir == 1)  // idir flag 1: Place in Y Row position
				   E[r * 2 + 1] = Ef[ c ].F;	
			   else			          // idir flag 0: Place in X Row position
				   E[r * 2] = Ef[ c ].F;
		   }
      }	   
   }



   // INVERT THE M MATRIX (Provided by dave: nxn.cpp)
   InverseNxN(M, N, Minv);
   
   //********** solve system here**********************************************************************************

   // SOLVE AND PRINT SOLUTION

   // Fill solution array with zero (define size of array)
   for (r = 0; r < (ts.NJ*2); r++)      
   {
	    for (c = 0; c < (ts.NM + ts.NR) ; c++) 
          {
          S[r] = 0;
          }
   }

   // Use matrix multiplication to solve X equations and X unknowns
   printf("\n\nSolution\n");
   printf("--------\n");
   for (r = 0; r < (ts.NJ*2); r++)      
   {
	    for (c = 0; c < (ts.NM + ts.NR) ; c++) 
          {
          S[r] += (Minv[r][c] * E[c]); // Matrix multiplication (Fill in solutions)
          }
          
          // Printing Member Solutions
          if (r < ts.NM)
               printf("Member %02d: F= %.2f\t %s [%s]\n", r, (S[r] < 0 ? S[r]*-1:S[r]), strUnits, (S[r] < 0 ? "T":"C") );
          // Print Reaction Solutions
          else 
             {     
               printf("Reaction on joint %02d = %s%.2f\t %s", Rf[r-ts.NM].j, ((S[r]*-1) < 0 ? "":"+"), S[r]*-1, strUnits);
               printf("(Rx= %s%.02lf %s", ((S[r]*-1)*cos(Rf[r-ts.NM].theta) <= 0 ? "":"+"), (S[r]*-1)*cos((Rf[r-ts.NM].theta)) , strUnits); // Rx = R Cos theta
               printf(" Ry= %s%.02lf %s", ((S[r]*-1)*sin(Rf[r-ts.NM].theta) <= 0 ? "":"+"), (S[r]*-1)*sin((Rf[r-ts.NM].theta)) , strUnits); // Ry = R Sin theta
               printf(", theta= %.2lf)\n",((Rf[r-ts.NM].theta)/(DEG2RAD)) );    
             }
   }

   // CLEAN UP DATA BUFFERS
   cleanup(f, Jf, Mf, Rf, Ef, M, Minv, E, S, N);
   return 0;
}

//-------------------------------------FUNCTIONS--------------------------------------------------------------------------------------------------------------------------


void cleanup(FILE *f, JOINT* Jf, MEMBER* Mf, FORCE* Rf, FORCE* Ef,
   double **M, double **Minv, double *E, double *S, int N)
{
   int i;

   for(i=0; i<N; i++)
   {
      delete[] M[i];
      delete[] Minv[i];
   }
   delete[] M;
   delete[] Minv;
   delete[] E;
   delete[] S;

   delete[] Jf;
   delete[] Mf;
   delete[] Ef;
   delete[] Rf;

   fclose(f);

   printf("\npress ENTER to end the program...");
   getchar();
}



//---------------------------------------------------------------------------------------

TRUSS_SIZES getTrussSizes(FILE *f) // gets file pointer, returns a truss size structure 
{
   TRUSS_SIZES ts;
   char buff[MAXLINE]; 
   char *header[4]={JOINT_COORDINATE_HEADER, MEMBER_CONNECTIVITY_HEADER,
      REACTIONS_HEADER, EXTERNAL_FORCES_HEADER};
   int i, size[4];


   for(i=0; i<4; i++) // loop to get data under each header
   {
      size[i]=-1;                
      rewind(f);     // sets file pointer to begining of the file 
      while(fgets(buff, MAXLINE, f)!=0)
      {
         if(strstr(buff, header[i])!=NULL)
         {
            size[i]=0;
            while(fgets(buff, MAXLINE, f)!=0)
            {
               if(strstr(buff, "[")) break;
               if(isEmpty(buff)) continue;
               size[i]++;
            }
            break;
         }
      }
      if(size[i]<1)
      {
         if(size[i]==-1)
            printf("cannot find header: %s\n", header[i]);
         else
            printf("No data for: %s\n", header[i]);
      }
   }

   ts.NJ =size[0];
   ts.NM =size[1];
   ts.NR =size[2];
   ts.NEF=size[3];
   return ts;
}

//---------------------------------------------------------------------------------------

bool getTrussData(FILE *f, TRUSS_SIZES ts, JOINT* Jf, MEMBER* Mf, FORCE* Ef, FORCE* Rf, char *strUnits)
{
   char ch, buff[MAXLINE], buff2[MAXLINE], *seps="\t ,;:";
   int i, j;

   //---------------------- joint coordinate data
   rewind(f);
   while(fgets(buff, MAXLINE, f)!=NULL)
   {
      if(strstr(buff, JOINT_COORDINATE_HEADER)!=NULL)  //look for the data header
      {
         i=0;
         while(i<ts.NJ)  // read each line, ignore blanks
         {
            if(fgets(buff, MAXLINE, f)==NULL || strstr(buff, "[")!=NULL)
            {
               printf("%s data corrupted.\n");
               return false;
            }
            if(isEmpty(buff)) continue;  // skip blank lines
            // finally got to a data line
            j=atoi(strtok(buff, seps));  // this is joint index.  Maybe out of order
            Jf[j].p[X]=atof(strtok(NULL, seps));  // x-coordinate
            Jf[j].p[Y]=atof(strtok(NULL, seps));  // y-coordinate
            i++;
         }
         break;
      }
   }

   //---------------------- member joint connectivity data
   rewind(f);
   while(fgets(buff, MAXLINE, f)!=NULL)
   {
      if(strstr(buff, MEMBER_CONNECTIVITY_HEADER)!=NULL)  //look for the data header
      {
         i=0;
         while(i<ts.NM)  // read each line, ignore blanks
         {
            if(fgets(buff, MAXLINE, f)==NULL || strstr(buff, "[")!=NULL)
            {
               printf("%s data corrupted.\n");
               return false;
            }
            if(isEmpty(buff)) continue;  // skip blank lines
            // finally got to a data line
            j=atoi(strtok(buff, seps));  // this is member index.  Maybe out of order
            Mf[j].j[0]=atoi(strtok(NULL, seps));  // joint index of one end
            Mf[j].j[1]=atoi(strtok(NULL, seps));  // joint index of opposite end
            i++;
         }
         break;
      }
   }

   //---------------------- reaction force data
   rewind(f);
   while(fgets(buff, MAXLINE, f)!=NULL)
   {
      if(strstr(buff, REACTIONS_HEADER)!=NULL)  //look for the data header
      {
         i=0;
         while(i<ts.NR)  // read each line, ignore blanks
         {
            if(fgets(buff, MAXLINE, f)==NULL || strstr(buff, "[")!=NULL)
            {
               printf("%s data corrupted.\n");
               return false;
            }
            if(isEmpty(buff)) continue;  // skip blank lines
            // finally got to a data line
            j=atoi(strtok(buff, seps));        // this is reaction index.  Maybe out of order
            Rf[j].j=atoi(strtok(NULL, seps));   // joint index of reaction
            Rf[j].F=0.0;                       // computed as part of solution
            strcpy(buff2, strtok(NULL, seps));  // text should be X or Y (reaction direction);
            ch=buff2[0];
            if(ch=='X' || ch=='x')
               Rf[j].idir=X;
            else if(ch=='Y' || ch=='y')
               Rf[j].idir=Y;
            else
            {
               //printf("Angle detected in reaction force data\n");
               Rf[j].theta=DEG2RAD*atof(buff2);
            }
            i++;
         }
         break;
      }
   }

   //---------------------- external force data
   rewind(f);
   while(fgets(buff, MAXLINE, f)!=NULL)
   {
      if(strstr(buff, EXTERNAL_FORCES_HEADER)!=NULL)  //look for the data header
      {
         i=0;
         while(i<ts.NEF)  // read each line, ignore blanks
         {
            fgets(buff, MAXLINE, f);
            if(isEmpty(buff)) continue;  // skip blank lines
            // finally got to a data line
            j=atoi(strtok(buff, seps));        // this is external force index.  Maybe out of order
            Ef[j].j=atoi(strtok(NULL, seps));  // joint index of external force
            Ef[j].F=atof(strtok(NULL, seps));  // force magnitude (signed)
            strcpy(buff2, strtok(NULL, seps));  // text should be X or Y (reaction direction);
            ch=buff2[0];
            if(ch=='X' || ch=='x')
               Ef[j].idir=X;
            else if(ch=='Y' || ch=='y')
               Ef[j].idir=Y;
            else
            {
               printf("Angle detected in external force data\n");
               Ef[j].theta=DEG2RAD*atof(buff2);
            }
            i++;
         }
         break;
      }
   }

   //---------------------- UNITS ---------------
   rewind(f);
   while(fgets(buff, MAXLINE, f)!=NULL)
   {
      if(strstr(buff, FORCE_UNITS_HEADER)!=NULL)  //look for the data header
      {
         while(fgets(buff, MAXLINE, f)!=NULL)
         {
            if(isEmpty(buff))
               continue;  // skip blank lines
            else
               break;
         }
         // finally got to a data line
         strcpy(strUnits, buff);
         for(i=0; i<(int)strlen(strUnits); i++) // get rid of newline character
         {
            if(strUnits[i]=='\n') strUnits[i]='\0';
         }
         break;
      }
   }
   return true;
}

//---------------------------------------------------------------------------------------

bool isEmpty(char *str)
{
   int i;
   for(i=0; i<(int)strlen(str); i++)
   {
      if(!isspace(str[i])) return false;
   }
   return true;
}
