/*
############################################################################
# Author:
# Date:
# Purpose:
# 
# 
# Related files: MCMC.c, simulation.R, MCMCdynoNEMdocu.tex, functionsMCMC.h
#		 netlearn.h, LogfilesMCMCrun.txt, and ouputMCMCdynoNEM
# 		 directory
############################################################################
*/

#include "netlearn.h"

//#####################

//void dynoNEM_posteriorEGenePos(double* Psi_R, int* nsgenes, int* negenes, int* T, double* D_R, double* Egene_prior_R, int* type, int* nrep, double* alpha, double* beta, double* res_posteriorPos);

/*void dynoNEM_posteriorEGenePos(double* Psi_R, int* nsgenes, int* negenes, int* T, double* D_R, double* Egene_prior_R, int* type, int* nrep, double* alpha, double* beta, double* res_posteriorPos){
	double*** D = (double***) R_alloc(*T, sizeof(double**));
	double** Egene_prior = (double**) R_alloc(*negenes, sizeof(double*));
	double** Psi = (double**) R_alloc(*nsgenes, sizeof(double*));
    	int t, i, s, k;
    	double** likelihood_post;
	for(t = 0; t < *T; t++){
		D[t] = (double**) R_alloc(*negenes, sizeof(double*));
		for(i = 0; i < *negenes; i++){
		    D[t][i] = (double*) R_alloc(*nsgenes, sizeof(double));
		    Egene_prior[i] = (double*) R_alloc(*nsgenes, sizeof(double));
		    for(s = 0; s < *nsgenes; s++){
			D[t][i][s] = D_R[s * *T * *negenes + i * *T + t];
			Egene_prior[i][s] = Egene_prior_R[s * *nsgenes + i];
		    }
		}
	}
	for(s = 0; s < *nsgenes; s++){
        	Psi[s] = (double*) R_alloc(*nsgenes, sizeof(double));
	        for(k = 0; k < *nsgenes; k++){
            		Psi[s][k] = Psi_R[k * *nsgenes + s];
        	}
    	}
    	likelihood_post = posteriorEGenePos(Psi, *nsgenes, *negenes, *T, D, Egene_prior, *type, *nrep, *alpha, *beta);
	for(i = 0; i < *negenes; i++){
        	for(s = 0; s < *nsgenes; s++){
            		res_posteriorPos[s * *nsgenes + i] = likelihood_post[i][s];
        	}
		free(likelihood_post[i]);
    	}	
	free(likelihood_post);	
}*/

//double** MCMCrun(long SAMPLE, long BURNIN, int TIMELAG, double THRESHOLD, double** net, int nsgenes, int negenes, int T, double*** D, double **Egene_prior, int type, int nrep, double alpha, double beta, time_t* seed, double* DIC, double* allLikelihoods);

SEXP dynoNEM_getPerturbProb(SEXP Psi_R, SEXP T_R, SEXP nsgenes_R, SEXP k_R);
//#####################

SEXP dynoNEM_getPerturbProb(SEXP Psi_R, SEXP T_R, SEXP nsgenes_R, SEXP k_R){    
    int T = INTEGER(T_R)[0];
    int nsgenes = INTEGER(nsgenes_R)[0];
    int k = INTEGER(k_R)[0];
    double** Psi = (double**) R_alloc(nsgenes, sizeof(double*));
    double** perturb_prob = (double**) R_alloc(nsgenes, sizeof(double*));
    int s, q, t;
    for(s = 0; s < nsgenes; s++){
        Psi[s] = (double*) R_alloc(nsgenes, sizeof(double));
	
	perturb_prob[s] = (double*) R_alloc(nsgenes, sizeof(double));// calloc changed to R_alloc	
        for(q = 0; q < nsgenes; q++){
            Psi[s][q] = (double) REAL(Psi_R)[q * nsgenes + s];	    
	}
    }
    double** perturb_prob2 = getPerturbProb(Psi, T, nsgenes, k, perturb_prob);
    SEXP res;
    PROTECT(res = NEW_NUMERIC(T * nsgenes));
    for(s = 0; s < nsgenes; s++){
        for(t = 0; t < T; t++){
            NUMERIC_POINTER(res)[t * nsgenes + s] = perturb_prob2[s][t];
        }
    }
    UNPROTECT(1);
    return(res);
}

//##############

SEXP MCMCrunWrapper(SEXP SAMPLE, SEXP BURNIN, SEXP initial_R, SEXP nsgenes_R, SEXP negenes_R, SEXP T_R, SEXP D_R, SEXP networkPrior_R, SEXP Egene_prior_R, SEXP priorScale_R, SEXP theta_R, SEXP type_R, SEXP nrep_R, SEXP alpha_R, SEXP beta_R, SEXP seed_R);

//#####################

// consider new arguments

SEXP MCMCrunWrapper(SEXP SAMPLE, SEXP BURNIN, SEXP initial_R, SEXP nsgenes_R, SEXP negenes_R, SEXP T_R, SEXP D_R, SEXP networkPrior_R, SEXP Egene_prior_R, SEXP priorScale_R, SEXP theta_R, SEXP type_R, SEXP nrep_R, SEXP alpha_R, SEXP beta_R, SEXP seed_R){
//long* sample; long* burnin;
    long sample = (long) INTEGER(SAMPLE)[0];
    long burnin = (long) INTEGER(BURNIN)[0];
    int T = INTEGER(T_R)[0];
    int nsgenes = INTEGER(nsgenes_R)[0];
    int negenes = INTEGER(negenes_R)[0];
    double priorScale = REAL(priorScale_R)[0];    
    double theta = REAL(theta_R)[0];
    int type = INTEGER(type_R)[0];
    int nrep = INTEGER(nrep_R)[0];
    double alpha = REAL(alpha_R)[0];
    double beta = REAL(beta_R)[0];
    int seed = INTEGER(seed_R)[0];
    int useMCMC = (sample > 0 & burnin > 0);
    if(!useMCMC){
	sample = 0;
	burnin = 0;
    }
  
    double*** D = (double***) R_alloc(T, sizeof(double**));
    double** initial = (double**) R_alloc(nsgenes, sizeof(double*));
    double** sdMat = (double**) R_alloc(nsgenes, sizeof(double*));
    //
    double** mat = (double**) R_alloc(nsgenes, sizeof(double*));
    double* allLikelihoods = (double*) R_alloc(sample + burnin + 1, sizeof(double));
    //
    double** networkPrior = NULL;
    if(networkPrior_R != NULL)
	    networkPrior = (double**) R_alloc(nsgenes, sizeof(double*));
    double** Egene_prior = (double**) R_alloc(negenes, sizeof(double*));
    int t, i, s, k;
    //double likelihood;
    for(t = 0; t < T; t++){
        D[t] = (double**) R_alloc(negenes, sizeof(double*));
        for(i = 0; i < negenes; i++){
            D[t][i] = (double*) R_alloc(nsgenes, sizeof(double));
            Egene_prior[i] = (double*) R_alloc(nsgenes + 1, sizeof(double));
            for(s = 0; s < nsgenes; s++){
                D[t][i][s] = REAL(D_R)[i * T + t + s * T * negenes]; // To convert one dim vector of R to 3D array in C 		
                Egene_prior[i][s] = REAL(Egene_prior_R)[s * negenes + i];// To convert one dim vector of R to 2D array in C	(2. index mal #Zeilen + 1. index)	==> BUG FIXED	
            }
	    Egene_prior[i][nsgenes] = REAL(Egene_prior_R)[nsgenes * negenes + i];
        }
    }
    //
    
    //
    for(s = 0; s < nsgenes; s++){
        initial[s] = (double*) R_alloc(nsgenes, sizeof(double));
	networkPrior[s] = (double*) R_alloc(nsgenes, sizeof(double));
        //
	sdMat[s] = (double*) R_alloc(nsgenes, sizeof(double));
        mat[s] = (double*) R_alloc(nsgenes, sizeof(double));
        //
        for(k = 0; k < nsgenes; k++){
            initial[s][k] = REAL(initial_R)[k * nsgenes + s];
	    if(networkPrior_R != NULL)
	            networkPrior[s][k] = REAL(networkPrior_R)[k * nsgenes + s];// To convert one D vector of R to 2D array in C
	    sdMat[s][k] = 0;
	    mat[s][k] = 0;
        }
    }
    //
    Rprintf("(initial) prior scale parameter = %g\n", priorScale);
    double loglik = learn_network( T, nsgenes,negenes, D, initial, networkPrior, Egene_prior, priorScale, mat, type, nrep, alpha, beta);    
    if(useMCMC){ // if proper arguments are given, use MCMC		
		for(s = 0; s < nsgenes; s++){
			for(k = 0; k < nsgenes; k++)
				initial[s][k] = mat[s][k];
		}	   
		MCMCrun( sample, //1
			  burnin,  //2		 
			  initial, //5
			  nsgenes,  //6
			  negenes,  //7
			  T, //8
			  D, //9
			  networkPrior, // 10
			  Egene_prior,  //11
			  priorScale,
			  theta,
			  type,  //12
			  nrep,  //13
			  alpha,  //14
			  beta, //15
			  seed,  //16
			  allLikelihoods, 
			  sdMat, 
			  mat); 

	    Rprintf("Sampling finished\n");	    
    }
    else{
	allLikelihoods[0] = loglik;
    }
    SEXP res_network, res_networkSD, res_loglik, result, wnames;
    PROTECT(res_network = NEW_NUMERIC(nsgenes * nsgenes));
    PROTECT(res_networkSD = NEW_NUMERIC(nsgenes * nsgenes));
    for(s = 0; s < nsgenes; s++){
	for(k = 0; k < nsgenes; k++){
	    NUMERIC_POINTER(res_network)[k * nsgenes + s] = mat[s][k];
	    NUMERIC_POINTER(res_networkSD)[k * nsgenes + s]= sdMat[s][k];
	    //Rprintf("%lf #\n",network[s][k]);
	}
    }
    PROTECT(res_loglik = NEW_NUMERIC(burnin + sample + 1));
    for(long counter = 0; counter < burnin + sample + 1; counter++){
	NUMERIC_POINTER(res_loglik)[counter] = allLikelihoods[counter];
    }
    PROTECT(result = NEW_LIST(3));
    PROTECT(wnames = NEW_CHARACTER(3));    
    SET_STRING_ELT(wnames, 0, mkChar("net.res"));
    SET_STRING_ELT(wnames, 1, mkChar("sdmat"));
    SET_STRING_ELT(wnames, 2, mkChar("allLikelihoods"));
    SET_NAMES(result, wnames);
    SET_ELEMENT(result, 0, res_network);
    SET_ELEMENT(result, 1, res_networkSD);
    SET_ELEMENT(result, 2, res_loglik);
    UNPROTECT(5);
    return(result);
}// ### Wrapper OK passes the values

