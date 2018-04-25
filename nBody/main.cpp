#include <stdio.h>
#include <stdlib.h>

#include <string.h>

//#include <omp.h>
#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

double G = 0.0000000000667384;

vec3 forceKonQ(vec3 posq, vec3 posk, double mk);

struct Body{
	//position
	vec3 p;
	//velocity
	vec3 v;
	//net acceleration
	vec3 a;
	//mass
	double m;
};

//The force of particle k on particle q
/*vec3 forceKonQ(vec3 posq, vec3 posk, double mk) {
	vec3 diff = posq - posk;
	double mag = diff.Magnitude;
	return diff*(mk/mag);
}*/

/*vec3 randVelocity(double min, double max) {
	double
}*/

//Random number in range
double genRand(double seed, double min, double max) {
	double r = (seed)*(max - min) + min;
	return r;
}

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
	int numParticlesMedium = 0;
	int numParticlesHeavy = 0;

	int numSteps = 0;
	int subSteps = 0;
	double timeSubStep;

	int width, height;

	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){
		numParticlesLight = atoi(argv[1]);
		numParticlesMedium = atoi(argv[2]);
		numParticlesHeavy = atoi(argv[3]);

		numSteps = atoi(argv[4]);
		subSteps = atoi(argv[5]);
		timeSubStep = atoi(argv[6]);

		width = atoi(argv[7]);
		height = atoi(argv[8]);

		int zDepth = 100;

		Body *particles;
		//Body *temp;
		vec3 *accel;

		image = (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
		//initialize image (all pixels to 0)
		for(int i = 0; i < (width * height * 3); i++){
			image[i] = 0;
		}

		//initialize light particles (randomized)
		int totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;
		particles = (Body *)malloc(totalParticles * sizeof(Body));
		//temp = (Body *)malloc(totalParticles * sizeof(Body));
		accel = (vec3 *)malloc(totalParticles * sizeof(vec3));

		for(int i = 0; i < numParticlesLight; i++){
			//mass
			particles[i].m = genRand(drand48(), massLightMin, massLightMax);
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));
			//velocity
			if(i%2 == 0){
				particles[i].v.x = genRand(drand48(), velocityLightMin, velocityLightMax);
				particles[i].v.y = genRand(drand48(), velocityLightMin, velocityLightMax);
				particles[i].v.z = genRand(drand48(), velocityLightMin, velocityLightMax);
			}else{
				particles[i].v.x = -1 * genRand(drand48(), velocityLightMin, velocityLightMax);
				particles[i].v.y = -1 * genRand(drand48(), velocityLightMin, velocityLightMax);
				particles[i].v.z = -1 * genRand(drand48(), velocityLightMin, velocityLightMax);
			}	
		}
		//initialize medium particles (randomized)
		for(int i = numParticlesLight; i < (numParticlesLight + numParticlesMedium); i++){
			//mass
			particles[i].m = genRand(drand48(), massMediumMin, massMediumMax);
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));
			//velocity
			if(i%2 == 0){
				particles[i].v.x = genRand(drand48(), velocityMediumMin, velocityMediumMax);
				particles[i].v.y = genRand(drand48(), velocityMediumMin, velocityMediumMax);
				particles[i].v.z = genRand(drand48(), velocityMediumMin, velocityMediumMax);
			}else{
				particles[i].v.x = -1 * genRand(drand48(), velocityMediumMin, velocityMediumMax);
				particles[i].v.y = -1 * genRand(drand48(), velocityMediumMin, velocityMediumMax);
				particles[i].v.z = -1 * genRand(drand48(), velocityMediumMin, velocityMediumMax);
			}
		}
		//initialize heavy particles (randomized)
		for(int i = (numParticlesLight + numParticlesMedium); i < totalParticles; i++){
			//mass
			particles[i].m = genRand(drand48(), massHeavyMin, massHeavyMax);
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));
			//velocity
			if(i%2 == 0){
				particles[i].v.x = genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
				particles[i].v.y = genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
				particles[i].v.z = genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
			}else{
				particles[i].v.x = -1 * genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
				particles[i].v.y = -1 * genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
				particles[i].v.z = -1 * genRand(drand48(), velocityHeavyMin, velocityHeavyMax);
			}
		} 
		//frame generation (per step)
		for(int frame = 0; frame < numSteps; frame++){
			//reset image pixels
			for(int i = 0; i < (width * height * 3); i++){
				image[i] = 0;
			}
			//sub step calculations
			for(int step = 0; step < subSteps; step++){
				//reset acceleration for particles
				for(int i = 0; i < totalParticles; i++){
					accel[i].x = 0;
					accel[i].y = 0;
					accel[i].z = 0;
				}
				
				//calculate acceleration from force of k exerted on q
				for(int i = 0; i < totalParticles; i++){
					Body q = particles[i];
					for(int j = 0; j < totalParticles; j++){
						if(i == j) continue;
						Body k = particles[j];

						vec3 dp = q.p - k.p;
						double dist = dp.Magnitude();

						double magQ = k.m / (dist * dist * dist);
						accel[i].x -= G * magQ * dp.x;
						accel[i].y -= G * magQ * dp.y;
						accel[i].z -= G * magQ * dp.z;

						double magK = q.m / (dist * dist * dist);
						accel[j].x += G * magK * dp.x;
						accel[j].y += G * magK * dp.y;
						accel[j].z += G * magK * dp.z;
					}
				}
				//update position & velocity for each particle based on acceleration
				for(int i = 0; i < totalParticles; i++){
					particles[i].v.x += accel[i].x * timeSubStep;
					particles[i].v.y += accel[i].y * timeSubStep;
					particles[i].v.z += accel[i].z * timeSubStep;

					particles[i].p.x += particles[i].v.x * timeSubStep;
					particles[i].p.y += particles[i].v.y * timeSubStep;
					particles[i].p.z += particles[i].v.z * timeSubStep;
				}
			}

			//light(red), medium(green), heavy(blue)
			int offset;
			for(int i = 0; i < totalParticles; i++){
				if(i < numParticlesLight){
					offset = (width * (int)particles[i].p.y + (int)particles[i].p.x) * 3 + 2;
				}else if(i >= numParticlesLight && i < (numParticlesLight + numParticlesMedium)){
					offset = (width * (int)particles[i].p.y + (int)particles[i].p.x) * 3 + 1;
				}else{
					offset = (width * (int)particles[i].p.y + (int)particles[i].p.x) * 3 + 0;
				}
				image[offset] = 255;
			}

			//almost done, just save the image
			std::string filename = std::string(argv[9]) + '_';
			std::string frameNum = std::to_string(frame);
			filename.append(5 - frameNum.length(), '0');
			filename = filename + frameNum + std::string(".bmp");
			printf("%s", filename.c_str());
			saveBMP(filename.c_str(), image, width, height);
		}
			
	}
	//all other nodes do this
	else{

	}

	free(image);

	MPI_Finalize();
	return 0;
}
