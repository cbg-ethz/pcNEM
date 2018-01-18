/*
 * netlearn.h
 *
 *  Created on: Aug 2, 2010
 *      Author: frohlich
 */

#ifndef NETLEARN_H_
#define NETLEARN_H_

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <float.h>
#include <time.h>
#include <string.h>

#define PVAL_DENS 0
#define EFFECT_PROB 1
#define DISCRETE 2


double** getPerturbProb(double** Psi, int T, int nsgenes, int k, double** perturb_prob);
double updateFactor (double likLogOld, double logPriorOld, double logPriorScaleOld, double likLogNew, double logPriorNew, double logPriorScaleNew, int oldNeighborhood, int newNeighborhood);
int neighborhoodSize(double** net, int nsgenes, int T);
double network_likelihood (double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta, double*** perturb_prob, double* loglik0); 
double logPrior(int nsgenes, double** net, double** prior, double inv_nu);
double logPriorLambda(double inv_nu, double theta);
void copyNet(int nsgenes, double** net, double** netCopy);//
int alterNet(double** net, int nsgenes, int T, double** temp1);
void MCMCrun(long sample, long burnin, double** net, int nsgenes, int negenes, int T, double*** D, double** networkPrior,  double **Egene_prior, double priorScale, double theta, int type, int nrep, double alpha, double beta, int seed, double* allLikelihoods, double** sdMat, double** matrix_r); // Egene_prior ????

double learn_network(int T, int nsgenes,int negenes, double*** D, double** initial, double** network_prior, double** Egene_prior, double prior_scale, double **net, int type, int nrep, double alpha, double beta);
void print_network(double** net, int nsgenes);

//double** posteriorEGenePos(double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta);
#endif /* NETLEARN_H_ */
