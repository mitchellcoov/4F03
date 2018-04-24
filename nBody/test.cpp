#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <omp.h>
#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

double G = 0.0000000000667384;

vec3 forceKonQ(vec3 posq, vec3 posk, double mk);

typedef double vect_t[DIM];  /* Vector type for position, etc. */

struct particle_s {
   double m;  /* Mass     */
   vect_t s;  /* Position */
   vect_t v;  /* Velocity */
};

int main(int argc, char* argv[]){

	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = 0;
	int numParticleMedium = 0;
	int numParticleHeavy = 0;

	int numSteps = 0;
	int subSteps = 0;
	double timeSubStep;

	int width, height;

	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){
		numParticlesLight = stoi(argv[1]);
		numParticlesMedium = stoi(argv[2]);
		numParticlesHeavy = stoi(argv[3]);

		numSteps = stoi(argv[4]);
		subSteps = stoi(argv[5]);
		timeSubStep = stoi(argv[6]);

		width = stoi(argv[7]);
		height = stoi(argv[8]);
		int depth = 10;

		printf("numParticlesLight: %d\nnumParticlesMedium: %d\nnumParticlesLarge: %d\nnumSteps: %d\nsubSteps: %d\ntimeSubStep: %f\nWidth: %d\nHeight: %d\n",numParticlesLight,numParticleMedium,numParticleHeavy,numSteps,numSubSteps,timeSubStep, width,height);


		int totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;


		Body particles[totalParticles];
		int * w = (int *) malloc(sizeof(int) * numParticlesTotal);
		double * s_x = (double *) malloc(sizeof(double) * numParticlesTotal);
 		double * s_y = (double *) malloc(sizeof(double) * numParticlesTotal);
 		double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal);
 		double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal);

		#pragma omp parallel num_threads(3)
		{
			switch (omp_get_thread_num()) {
				case 0:
					for (int i = 0; i < numParticlesLight; i++) {
						// w[i] =rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
						// s_x[i]=drand48()*(width-1);
						// s_y[i]=drand48()*(height-1);
						// v_x[i]=rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						// v_y[i]=rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						particles[i].m = rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
						particles[i].p[X] = rand48()*(width-1);
						particles[i].p[Y] = rand48()*(height-1);
						particles[i].p[Z] = rand48()*(depth);
						particles[i].v[X] = rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						if (part % 2 == 0){
							particles[i].v[Y] =rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						}else{
							particles[i].v[Y] =(-1)*rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						}
						particles[i].m = rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
					}
				case 1:
					for (int i = numParticlesLight; i < numParticlesMedium + numParticlesLight; i++) {
						// w[i] =rangeRand((double)drand48(), (double)properties::massMediumMin, (double)properties::massMediumMax);
						// s_x[i]=drand48()*(width-1);
						// s_y[i]=drand48()*(height-1);
						// v_x[i]=rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						// v_y[i]=rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						particles[i].m = rangeRand((double)drand48(), (double)properties::massMediumMin, (double)properties::massMediumMax);
						particles[i].p[X] = rand48()*(width-1);
						particles[i].p[Y] = rand48()*(height-1);
						particles[i].p[Z] = rand48()*(depth);
						particles[i].v[X] = rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						if (part % 2 == 0){
							particles[i].v[Y] =rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						}else{
							particles[i].v[Y] =(-1)*rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						}
					}
				case 2:
					for (int i = numParticleMedium + numParticlesLight; i < totalParticles; i++) {
						// w[i] =rangeRand((double)drand48(), (double)properties::massHeavyMin, (double)properties::massHeavyMax);
						// s_x[i]=drand48()*(width-1);
						// s_y[i]=drand48()*(height-1);
						// v_x[i]=rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						// v_y[i]=rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						particles[i].m = rangeRand((double)drand48(), (double)properties::massHeavyMin, (double)properties::massHeavyMax);
						particles[i].p[X] = rand48()*(width-1);
						particles[i].p[Y] = rand48()*(height-1);
						particles[i].p[Z] = rand48()*(depth);
						particles[i].v[X] = rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						if (part % 2 == 0){
							particles[i].v[Y] =rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						}else{
							particles[i].v[Y] =(-1)*rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						}
					}
			}
		}


			//almost done, just save the image
		saveBMP(argv[9], image, width, height);
	}
	//all other nodes do this
	else{

	}

	free(image);

	MPI_Finalize();
	return 0;
}

//The force of particle k on particle q
vec3 forceKonQ(vec3 posq, vec3 posk, double mk) {
	vec3 diff = posq - posk;
	double mag = diff.Magnitude;
	return diff*(mk/mag);
}


//Random number in range
double rangeRand(double rand, double start, double end) {
	double frand = (rand / RAND_MAX)*(end - start) + start;
	return frand;
}
