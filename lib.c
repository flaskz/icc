#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>

#include "lib.h"
#include "main.h"

//------------------------------------------------------------------------------//    
//                                                                              //
//   ██████╗ ██████╗ ███████╗███████╗ ██████╗ ██╗    ██╗   ██╗███████╗██████╗   //
//   ██╔══██╗██╔══██╗██╔════╝██╔════╝██╔═══██╗██║    ██║   ██║██╔════╝██╔══██╗  //
//   ██████╔╝██║  ██║█████╗  ███████╗██║   ██║██║    ██║   ██║█████╗  ██████╔╝  //
//   ██╔═══╝ ██║  ██║██╔══╝  ╚════██║██║   ██║██║    ╚██╗ ██╔╝██╔══╝  ██╔══██╗  //
//   ██║     ██████╔╝███████╗███████║╚██████╔╝███████╗╚████╔╝ ███████╗██║  ██║  //
//   ╚═╝     ╚═════╝ ╚══════╝╚══════╝ ╚═════╝ ╚══════╝ ╚═══╝  ╚══════╝╚═╝  ╚═╝  //
//                                                                              //
//   Lib - Version 2.8                                                          //
//                                                                              //
//------------------------------------------------------------------------------// 

/* ============================================================================
  			Funções Básicas
  ============================================================================*/
double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

void imprimeTempo(double start, double end, char *str){

        fprintf(file, "# %s: %1.8f\n", str, end-start);
        fprintf(stdout, "# %s: %1.8f\n", str, end-start);
}

void imprimeMat(long double *mat, int matWidth, int nx, int ny, double hx, double hy){

	fprintf(file, "set title \"Entrada: %s %s %s %s %s %s\"\n", ar1, ar2, ar3, ar4, ar5, ar6);
	fprintf(file, "# X   Y   Z   \nsplot \"-\" using 1:2:3 with lines\n");
	for(int i=0; i<= ny; ++i){
		for(int j=0; j<= nx; ++j){
			fprintf(file, "%lf   %lf   %1.8Lf   \n", i*hy, j*hx, mat[j * matWidth + i]); //troquei j por i
		}
		fprintf(file, "\n");
	}

}


void gaussSeidelSOR(long double *mat, int matWidth, long double *fxy, int nx, int ny, int nosInternos, double w, double *X, int itmax){

	 
	double aux, Raux, start=0.0, end=0.0, accStart=0.0, accEnd=0.0, norm_start, norm_end, accNormS=0.0, accNormE=0.0;
	int l=1;
//	long double matHolderRes=0.0;
//      int matHolder=0, fxyWidth=nx-1;

              

        long double res=0.0, normVec[itmax];
      
      
 	long double *matAux = (long double *) malloc(sizeof(long double)*nosInternos);      
        for(int i=0; i<nosInternos; i++)
                matAux[i]=0.0;
        

	for (l=0; l<itmax; l++){
        res = 0.0; 
        Raux = 0.0;
	//copia os nós internos de mat para matAux
        	for (int i=1; i<ny; ++i){
			for (int j=1; j<nx ; ++j){
				matAux[(i-1) * (nx-1) + (j-1)] = mat[i*matWidth+j];
			}
		}




        long double matHolderRes=0.0;
        int matHolder=0, fxyWidth=nx-1, matIndexHolder=0;;
	int indexIMHolder=0,  indexJMHolder=0; 
        double flexHolder= 1-w;
        //double x0=X[0],x1=X[1],x2=X[2],x3=X[3],x4=X[4];
        //--------------------------------------------------------//
                LIKWID_MARKER_START("SOR");
        //likiwid start marker
        //--------------------------------------------------------//

	//--------------------------------------------------------//
		start=timestamp();
	// Time Stamp START
	//--------------------------------------------------------//
	
 		for (int i=1; i<ny; ++i){
			for (int j=1; j<nx; ++j){
                                //matIndexHolder = matWidth+j;
                                //matHolder = i*matWidth+j;
                                matHolder = i*matWidth+j;
				indexIMHolder = i-1;

				//----------------Original--------------------------------------------//
				//aux = fxy[(i-1)*(nx-1)+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
	  			//mat[(i-1)*matWidth+j]*X[2] - mat[i*matWidth+(j+1)]*X[3] +
	  			//mat[i*matWidth+(j-1)]*X[4];

                                //----------------Opt-1-----------------------------------------------//
                                //aux = fxy[(i-1)*(nx-1)+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                //mat[(i-1)*matWidth+j]*X[2] - mat[matHolder+1]*X[3] +
                                //mat[matHolder-1]*X[4];

                                //----------------Opt-2-----------------------------------------------//				
                                aux = fxy[indexIMHolder*fxyWidth+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                mat[indexIMHolder*matWidth+j]*X[2] - mat[matHolder+1]*X[3] +
                                mat[matHolder-1]*X[4];

                                //----------------Opt-3-----------------------------------------------//                                
                                //aux = fxy[indexIMHolder*fxyWidth+(j-1)]*x0 - mat[(i+1)*matWidth+j]*x1 +
                                //mat[indexIMHolder*matWidth+j]*x2 - mat[matHolder+1]*x3 +
                                //mat[matHolder-1]*x4;


				//Raux= w*aux + (1-w)*mat[i*matWidth+j];
                                Raux= w*aux + (flexHolder)*mat[matHolder];
				// Atualiza matriz:
				mat[matHolder] = Raux;
				
      			}
    		}
	//--------------------------------------------------------//
		end=timestamp();
        // Time Stamp STOP
        //--------------------------------------------------------//
		accStart += start;
		accEnd += end;
		
        //--------------------------------------------------------//
                LIKWID_MARKER_STOP("SOR");
        //likiwid stop marker
        //--------------------------------------------------------//



        //--------------------------------------------------------//
                LIKWID_MARKER_START("RES");
        //likiwid start marker
        //--------------------------------------------------------//

	//Calcula o resíduo
        //--------------------------------------------------------//
                norm_start=timestamp();
        // Time Stamp START
        //--------------------------------------------------------//

                for (int i=1; i<ny; ++i){
                        for (int j=1; j<nx; ++j){
                                //matIndexHolder = matWidth+j;
                                matHolder = i*matWidth+j;  
				//matHolder = i*matIndexHolder;
                                indexIMHolder = i-1;
				indexJMHolder= j-1;
                                //aux = fxy[(i-1)*(nx-1)+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                //mat[(i-1)*matWidth+j]*X[2] - mat[i*matWidth+(j+1)]*X[3] +
                                //mat[i*matWidth+(j-1)]*X[4];

				//--------------Opt-1-----------------------------------------------//
                                //aux = fxy[(i-1)*(nx-1)+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                //mat[(i-1)*matWidth+j]*X[2] - mat[matHolder+1]*X[3] +
                                //mat[matHolder-1]*X[4];
                              
                                //--------------Opt-1-----------------------------------------------//
                                //aux = fxy[(i-1)*fxyWidth+(j-1)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                //mat[(i-1)*matWidth+j]*X[2] - mat[matHolder+1]*X[3] +
                                //mat[matHolder-1]*X[4];

                                //--------------Opt-2-----------------------------------------------//
                                aux = fxy[(indexIMHolder)*fxyWidth+(indexJMHolder)]*X[0] - mat[(i+1)*matWidth+j]*X[1] +
                                mat[(indexIMHolder)*matWidth+j]*X[2] - mat[matHolder+1]*X[3] +
                                mat[matHolder-1]*X[4];


                                //--------------Opt-3-----------------------------------------------//
                                //aux = fxy[(indexIMHolder)*fxyWidth+(indexJMHolder)]*x0 - mat[(i+1)*matWidth+j]*x1 +
                                //mat[(indexIMHolder)*matWidth+j]*x2 - mat[matHolder+1]*x3 +
                                //mat[matHolder-1]*x4;

                                //--------------------------------------------------------//
                                // Calculo do resíduo:

                                //matHolderRes = (matAux[(i-1)* (nx-1) + (j-1)] - aux);
                                //matHolderRes = (matAux[(i-1)* fxyWidth + (j-1)] - aux);
                                matHolderRes = (matAux[(indexIMHolder)* fxyWidth + (indexJMHolder)] - aux);
                                //res += (matAux[(i-1)* (nx-1) + (j-1)] - aux) * (matAux[(i-1)* (nx-1) + (j-1)] - aux) ;
                                res += (matHolderRes) * (matHolderRes) ;

                        }
              	}
 
		//Calcula a norma L^2 do resíduo ||res||;		
		
                normVec[l] = sqrt(res);

	//--------------------------------------------------------//
                norm_end=timestamp();
        // Time Stamp STOP
	//--------------------------------------------------------//
		accNormS += norm_start;
		accNormE += norm_end;

        //--------------------------------------------------------//
                LIKWID_MARKER_STOP("RES");
        //likiwid stop marker
        //--------------------------------------------------------//
  	}

	
	start = accStart/l+1;
	end = accEnd/l+1;
	fprintf(file, "###########\n");
	imprimeTempo(start,end,"Tempo Método SOR");
	
	norm_start = accNormS/l+1;
	norm_end = accNormE/l+1;
	imprimeTempo(norm_start,norm_end,"Tempo Resíduo");

	fprintf(file, "#\n# Norma do Resíduo\n");
        for(int i=0; i<itmax; ++i){
        	fprintf(file, "# i=%d: %1.8Lf\n", i+1, normVec[i]);
        }
        fprintf(file, "###########\n");

}

//double norma(long double *res, int nosInternos){
//	long double aux=0.0;
//	for(int i=0; i< nosInternos; ++i){
//		aux += res[i]*res[i]; 
//	}
//	return sqrt(aux);
//}


