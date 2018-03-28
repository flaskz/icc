#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
//   Main - Version 2.8                                                         //
//                                                                              //
//------------------------------------------------------------------------------//                       
int main(int argc, char *argv[]){
        //--------------------------------------------------------//
        	LIKWID_MARKER_INIT;
        //likiwid initialization marker
        //--------------------------------------------------------//



	if(argc < 9 ||  strcmp(argv[1], "-hx") ||  strcmp(argv[3], "-hy")||  strcmp(argv[5], "-i")||  strcmp(argv[7], "-o")  ){
		fputs("Uso correto: <pdeSolver> -hx <Hx> -hy <Hy> -i <maxIter> -o [arquivo_saida]\n",stderr);
		exit(-1);
	}

	file = fopen(argv[8], "w");
	if(file == NULL){
		fputs("Erro ao criar arquivo de saída\n",stderr);
		exit(-1);
	}

	//----------------------------------------------------
	//Calculo das constantes:
	//
	ar1 = argv[1]; ar2 = argv[2]; ar3 = argv[3];
	ar4 = argv[4]; ar5 = argv[5]; ar6 = argv[6];

	double hx = strtod(argv[2], NULL);//Transforma de string para double;
	double hy = strtod(argv[4], NULL);

	//if (hx < 0.0044 || hy < 0.0044){
	//  fputs("hx e hy devem ser maiores ou iguais a 0.0044\n",stderr);
	//  exit(-1);
	//}
	  
	
	if(hx <= 0.0 || hx > PI ||hy <= 0.0 || hy > PI){
		fputs("hx e hy devem pertencer ao intervalo ]0,PI] !\n",stderr);
		exit(-1);
	}

	int itmax = strtod(argv[6], NULL);

	int nx = round(PI/hx)+1; //Calcula nx e ny de acordo com a fórmla para discretizar o intervalo;
	int ny = round(PI/hy)+1;
	
	int nosInternos = (nx-1)*(ny-1); // quantidade de nós internos;
        int eleMatriz = (nx+1)*(ny+1); // elementos da matriz quadrada;

        fprintf(stdout,"Quantidade de pontos: %d \n", nx+1);

	hx = PI/nx; //faz ajustes do hx e hy caso o valor inserido pelo usuário seja inválido;
	hy = PI/ny;
	
		//Cálculo de todas as constantes da equação;
	double k5 = (4*hy*hy)+(4*hx*hx)+(8*PI*PI*hy*hy*hx*hx); //Alterado de 16 para 8
	double k1 = hy*hy*(hx-2)/k5;
	double k2 = hy*hy*(hx+2)/k5;
	double k3 = hx*hx*(hy-2)/k5;
	double k4 = hx*hx*(hy+2)/k5;
	double k6 = 2*hx*hx*hy*hy/k5;

	double X[5] = {k6, k1, k2, k3, k4}; //Vetor de constantes;
	//
	

	double w = 2-( (hx+hy)/2 ); //fator de relaxação;

	int matWidth = nx; //Largura da matriz //MUDEI AQUI! SOMEI 1 !!!
	//double *mat = calloc(eleMatriz, sizeof(double)); //Aloca memória e zera os valores;
        long double *mat = (long double *) malloc (eleMatriz*sizeof(long double));
        for(int i=0; i<eleMatriz; i++)
		mat[i]=0.0;

	// para fazer o acesso usar: mat[i * matWidth + j]; onde i percorre no eixo Y (linha) e j no eixo X (coluna);


	//----------------------------------------------------
	//Calculo das constantes das bordas envolvendo seno:
        //

		//---> calcula valores superior e inferior da matriz excluindo os cantos;
	for(int j=1; j< nx; ++j){ //Percorrendo as colunas na borda superior e inferior;
		mat[0 * matWidth + j] = sin(2.0*PI*(PI-(j*hx)))*sinh(PI*PI);
		mat[ny * matWidth + j] = sin(2.0*PI*j*hx)*sinh(PI*PI);       
	}

		//---> calcula os valores de F(x,y)*4*k6*PI^{2} para os nós internos
	long double *fxy = (long double *) malloc(sizeof(long double)*nosInternos);
        for(int i=0; i<nosInternos; i++)
                fxy[i]=0.0;        
        //fflush(stdout);
        //fprintf(stdout,"Nos internos: %d \n",nx); 
	//for(int i=0; i< ny-1; ++i){
        for(int i=0; i< ny-1; ++i){
                //for(int j=0; j< nx-1; ++j){
		for(int j=0; j< nx-1; ++j){// X é j : Y é i ---> U(j,i)
			//k6 removido do cálculo pq ele é feito no gauss-seidel
                        //fprintf(stdout,"valor dentro do fxy[]: %d \n", (i*(nx-1)+j));
      			fxy[i * (nx-1) + j] = 4*PI*PI*(sin(2*PI*(j+1)*hx)*sinh(PI*(i+1)*hy) + sin(2*PI*(PI-((j+1)*hx)))*sinh(PI*(PI-((i+1)*hy))));
    		}	
  	}
	//	


	//----------------------------------------------------
	//Executa método de Gauss-Seidel com Successive Over Relaxation:
	//
	gaussSeidelSOR(mat, matWidth,fxy, nx, ny, nosInternos, w, X,itmax);

	imprimeMat(mat,matWidth,nx,ny,hx,hy);
	fprintf(file, "end\npause -1\n");
	
	//--------------------------------------------------------//
        	LIKWID_MARKER_CLOSE;
        //likiwid initialization marker
        //--------------------------------------------------------//

	return 0;
}
