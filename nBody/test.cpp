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
#define DIM 2  /* Two-dimensional system */
#define X 0    /* x-coordinate subscript */
#define Y 1    /* y-coordinate subscript */
#define Z 2
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

		image = (unsigned char *)malloc(3*width*height*sizeof(unsigned char));
		//Create black image
		for (int b = 0; b < (3*width*height); b++)
		{
			image[b] = unsigned char 0;
		}

		particles_s particles[totalParticles];


		#pragma omp parallel num_threads(3)
		{
			switch (omp_get_thread_num()) {
				case 0:
					for (int i = 0; i < numParticlesLight; i++) {
					  particles[i].m = rangeRand((double)drand48(), (double)properties::massLightMin, (double)properties::massLightMax);
						particles[i].p[X] = drand48()*(width-1);
						particles[i].p[Y] = drand48()*(height-1);
						particles[i].p[Z] = drand48()*(depth);
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
										particles[i].m = rangeRand((double)drand48(), (double)properties::massMediumMin, (double)properties::massMediumMax);
						particles[i].p[X] = drand48()*(width-1);
						particles[i].p[Y] = drand48()*(height-1);
						particles[i].p[Z] = drand48()*(depth);
						particles[i].v[X] = rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						if (part % 2 == 0){
							particles[i].v[Y] =rangeRand((double)drand48(), (double)properties::velocityMediumMin, (double)properties::velocityMediumMax);
						}else{
							particles[i].v[Y] =(-1)*rangeRand((double)drand48(), (double)properties::velocityLightMin, (double)properties::velocityLightMax);
						}
					}
				case 2:
					for (int i = numParticleMedium + numParticlesLight; i < totalParticles; i++) {
						particles[i].m = rangeRand((double)drand48(), (double)properties::massHeavyMin, (double)properties::massHeavyMax);
						particles[i].p[X] = drand48()*(width-1);
						particles[i].p[Y] = drand48()*(height-1);
						particles[i].p[Z] = drand48()*(depth);
						particles[i].v[X] = rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						if (part % 2 == 0){
							particles[i].v[Y] =rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						}else{
							particles[i].v[Y] =(-1)*rangeRand((double)drand48(), (double)properties::velocityHeavyMin, (double)properties::velocityHeavyMax);
						}
					}
			}
		}

		for (int j = 0; j < totalParticles; g++){
			int Xval = (int) particles[j].p[X];
			int Yval = (int) particles[j].p[Y];
			image[(Xval*3)+Yval*width*3] = (unsigned char) particles[j];
			image[(Xval*3)+Yval*width*3+1] = (unsigned char) particles[j];
			image[(Xval*3)+Yval*width*3+2] = (unsigned char) particles[j];
		}

		char file[20];
		strcpy(file, argv[9]);
		strcat(file,"_000000.bmp");
		const char* filename = file;
		const unsigned char* result = (image);
		saveBMP (filename, result, width, height);

		//Computing forces here

		// for (step = 1; step <= n_steps; step++) {
		// 	 t = step*substeps;
		// 	 memset(forces, 0, n*sizeof(vect_t));
		// 	 for (int particle = 0; part < totalParticles-1; part++)
		// 		 for(int k= particle+1; k<totalParticles;k++){
		//
		// 		 }
		// 	 	for (int particle = 0; part < totalParticles; part++)
		// 			Update_part(part, forces, curr, n, delta_t);
		//
		// }
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
// vec3 forceKonQ(vec3 posq, vec3 posk, double mk) {
// 	vec3 diff = posq - posk;
// 	double mag = diff.Magnitude;
// 	return diff*(mk/mag);
// }
//Compute force of particles. Taken from Sample code from book
// void Compute_force(int particle, vect_t forces[], struct particle_s curr[],int n) {
//    int k;
//    double mg;
//    vect_t f_part_k;
//    double len, len_3, fact;
//
//    for (k = part+1; k < n; k++) {
//       /* Compute force on part due to k */
//       f_part_k[X] = curr[part].s[X] - curr[k].s[X];
//       f_part_k[Y] = curr[part].s[Y] - curr[k].s[Y];
//       len = sqrt(f_part_k[X]*f_part_k[X] + f_part_k[Y]*f_part_k[Y]);
//       len_3 = len*len*len;
//       mg = -G*curr[part].m*curr[k].m;
//       fact = mg/len_3;
//       f_part_k[X] *= fact;
//       f_part_k[Y] *= fact;
//
//       /* Add force in to total forces */
//       forces[part][X] += f_part_k[X];
//       forces[part][Y] += f_part_k[Y];
//       forces[k][X] -= f_part_k[X];
//       forces[k][Y] -= f_part_k[Y];
//    }
// }  /* Compute_force */
// //Random number in range
double rangeRand(double rand, double start, double end) {
	double frand = (rand / RAND_MAX)*(end - start) + start;
	return frand;
}
