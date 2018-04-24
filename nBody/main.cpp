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

typedef struct Body{
	//position
	vec3 p;
	//velocity
	vec3 v;
	//net force
	vec3 f;
	//mass
	double m;
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
		xpos = (double *)malloc(N * sizeof(double));
  	ypos = (double *)malloc(N * sizeof(double));
		zpos = (double *)malloc(N * sizeof(double));

		printf("numParticlesLight: %d\nnumParticlesMedium: %d\nnumParticlesLarge: %d\nnumSteps: %d\nsubSteps: %d\ntimeSubStep: %f\nWidth: %d\nHeight: %d\n",numParticlesLight,numParticleMedium,numParticleHeavy,numSteps,numSubSteps,timeSubStep, width,height);


		int totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;

		//Generate black image to put stuff on
		for(int a=0; a<(3*frameSize);a++){
				image[a]=(unsigned char) 0;
		}


		Body particles[totalParticles];
		#pragma omp parallel num_threads(3)
		{
			switch (omp_get_thread_num()) {
				case 0:
					for (int i = 0; i < numParticlesLight; i++) {
						// xpos[i] = drand48()*(width-1);
						// ypos[i] = drand48()*(height-1);
						// zpos[i] = depth;
						// particles[i].p = vec3(xpos[i],ypos[i],zpos[i],vec3 properties::colourLight);
						// particles[i].v = rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						// particles[i].m = rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
						particles[i].p = vec3(xpos[i],ypos[i],zpos[i],vec3 properties::colourLight);
						particles[i].v = rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						particles[i].m = rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
					}
				case 1:
					for (int i = numParticlesLight; i < numParticlesMedium + numParticlesLight; i++) {
						// xpos[i] = drand48()*(width-1);
						// ypos[i] = drand48()*(height-1);
						// zpos[i] = depth;
						// particles[i].p = vec3(xpos[i],ypos[i],zpos[i],vec3 properties::colourMedium);
						// particles[i].v = rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						// particles[i].m = rangeRand((double)drand48(), (double)properties::massMediumMin, (double)properties::massMediumMax);
					}
				case 2:
					for (int i = numParticleMedium + numParticlesLight; i < totalParticles; i++) {
						// xpos[i] = drand48()*(width-1);
						// ypos[i] = drand48()*(height-1);
						// zpos[i] = depth;
						// particles[i].p = vec3(xpos[i],ypos[i],zpos[i],vec3 properties::colourHeavy);
						// particles[i].v = rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						// particles[i].m = rangeRand((double)drand48(), (double)properties::massHeavyMin, (double)properties::massHeavyMax);
					}
			}
		}

		unsigned char* image = (unsigned char *) calloc(3*width*height, sizeof(unsigned char));

		for(int j=0; g < totalParticles; g++){
			image[j] = particles[j];
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

vec3 randVelocity(double min, double max) {
	double
}

//Random number in range
double rangeRand(double rand, double start, double end) {
	double frand = (rand / RAND_MAX)*(end - start) + start;
	return frand;
}
