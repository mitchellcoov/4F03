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
double rangeRand(double rand, double start, double end);

int main(int argc, char* argv[]){

	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep width height imageFilenamePrex\n", argv[0]);
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
		int numParticlesLight = std::stoi(argv[1]);
		int numParticlesMedium = std::stoi(argv[2]);
		int numParticlesHeavy = std::stoi(argv[3]);

		int numSteps = std::stoi(argv[4]);
		int subSteps = std::stoi(argv[5]);
		int timeSubStep = std::stoi(argv[6]);

		int width = std::stoi(argv[7]);
		int height = std::stoi(argv[8]);
		int totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;

		int * w = (int *) malloc(sizeof(int) * totalParticles);
 		int * s_x = (int *) malloc(sizeof(int) * totalParticles);
 		int * s_y = (int *) malloc(sizeof(int) * totalParticles);
 		double * v_x = (double *) malloc(sizeof(double) * totalParticles);
 		double * v_y = (double *) malloc(sizeof(double) * totalParticles);


		image = (unsigned char *)malloc(3*width*height*sizeof(unsigned char));
		//Create black image
		for (int b = 0; b < (3*width*height); b++){
			image[b] = (unsigned char) 0;
		}

		//GENERATE STUFF HERE
		#pragma omp parallel num_threads(3)
		{
			switch (omp_get_thread_num()) {
				case 0:
					for (int i = 0; i < numParticlesLight; i++) {
							//w[i] = drand48() * (massLightMax-massLightMin+1) + massLightMin;
							w[i] = (double) rangeRand(drand48(),massLightMin,massLightMax);
			 				s_x[i] = drand48() * width;
			 				s_y[i] = drand48() * height;

	 				if(i%2 ==0){
		 					v_x[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
		 					v_y[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
	 				} else{
		 					v_x[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
		 					v_y[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
	 				}
				}
				case 1:
					for (int i = numParticlesLight; i < numParticlesMedium + numParticlesLight; i++) {
						//w[i] = drand48() * (massMediumMax-massMediumMin+1) + massMediumMin;
							w[i] = (double) rangeRand(drand48(),massMediumMin,massMediumMax);
			  			s_x[i] = drand48()*width;
			 				s_y[i] = drand48()*height;
 				if(i%2 ==0){
		 					v_x[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
		 					v_y[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 				} else{
		 					v_x[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
		 					v_y[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 				}
					}
				case 2:
					for (int i = numParticleMedium + numParticlesLight; i < totalParticles; i++) {
						//w[i] = drand48() * (massHeavyMax-massHeavyMin+1) + massHeavyMin;
						w[i] = (double) rangeRand(drand48(),massHeavyMin,massHeavyMax);
		 				s_x[i] = drand48()*width;
		 				s_y[i] = drand48()*height;
		 				if(i%2 ==0){
		 					v_x[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
		 					v_y[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
		 				} else{
		 					v_x[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
		 					v_y[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
		 				}
					}
			}
		}
		unsigned char* image = (unsigned char *) calloc(3*width*height, sizeof(unsigned char));

		for(int i = 0; i < totalParticles; i++){

			int index = (s_y[i] * width + s_x[i])*3;
			if (index < (sizeof(unsigned char) *3*width*height) && index >= 0){
				if(w[i] >= 1 && w[i] <= 5){
					image[index] = 0;
					image[index+1] = 0;
					image[index+2] = 255;
					printf("derp1 \n");
				} else if(w[i] >= 6 && w[i] <= 10){
					image[index] =0;
					image[index+1] = 255;
					image[index+2] = 0;
					printf("derp2 \n");
				} else if(w[i] >= 11 && w[i]<=15) {
					image[index] = 255;
					image[index+1] = 0;
					image[index+2] = 0;
					printf("derp3 \n");
				}
			}
		}

		for (int step = 1; step <= totalParticles, step++){

			for (int particle = 0; i < totalParticles; particle++){

				for (int j=i+1; j<totalParticles; j++){

					forces = computeForce(s_x[i],s_y[i],w[i],s_x[j],s_y[j],w[j],0);

				}
			}

		}



			//almost done, just save the image
		saveBMP(argv[9], image, width, height);

	}

	else{

	}

	free(image);

	//MPI_Finalize();
	return 0;
}
//
// //The force of particle k on particle q
// vec3 forceKonQ(vec3 posq, vec3 posk, double mk) {
// 	vec3 diff = posq - posk;
// 	double mag = diff.Magnitude;
// 	return diff*(mk/mag);
// }
//
double rangeRand(double rand, double start, double end) {
	double frand = (rand / RAND_MAX)*(end - start) + start;
	return frand;
}
