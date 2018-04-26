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

struct Body{
	//position
	vec3 p;
	//velocity
	vec3 v;
	//mass
	double m;
};

//Random number in range
double genRand(double seed, double min, double max) {
	double r = (seed)*(max - min) + min;
	return r;
}

void computeAccel(vec3 *p, double *m, vec3 *accelBlock, int blockSize){
	
	for(int i = 0; i < blockSize; i++){
		//Body q = particles[i];
		for(int j = i + 1; j < blockSize; j++){
			//Body k = particles[j];

			vec3 dp = p[i] - p[j];
			double dist = dp.Magnitude();

			double magQ = m[j] / (dist * dist * dist);
			accelBlock[i].x -= G * magQ * dp.x;
			accelBlock[i].y -= G * magQ * dp.y;
			accelBlock[i].z -= G * magQ * dp.z;

			double magK = m[i] / (dist * dist * dist);
			accelBlock[j].x += G * magK * dp.x;
			accelBlock[j].y += G * magK * dp.y;
			accelBlock[j].z += G * magK * dp.z;
		}
	}
}

int main(int argc, char* argv[]){

	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
	}

	int numParticlesLight = atoi(argv[1]);
	int numParticlesMedium = atoi(argv[2]);
	int numParticlesHeavy = atoi(argv[3]);

	int numSteps = atoi(argv[4]);
	int subSteps = atoi(argv[5]);
	double timeSubStep = atoi(argv[6]);

	int width = atoi(argv[7]);
	int height = atoi(argv[8]);
	int totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;

	int p, my_rank;

	//variables

	Body particles[totalParticles];
	vec3 accel[totalParticles];
	int zDepth;

	//temporary (mpi) arrays
	vec3 pos[totalParticles];
	double mas[totalParticles];

	unsigned char image[width*height*3];	

	double start, end;
	double timing[(numSteps*subSteps)];
	int t = 0;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	MPI_Datatype mpi_vec3;
	MPI_Type_contiguous(3, MPI_DOUBLE, &mpi_vec3);
	MPI_Type_commit(&mpi_vec3);

	int blockSize = totalParticles / p;
	vec3 accelBlock[blockSize];
	

	//root node stuff goes here
	if(my_rank == 0){


		zDepth = 100;

		//initialize image (all pixels to 0)
		for(int i = 0; i < (width * height * 3); i++){
			image[i] = 0;
		}

		//initialize light particles (randomized)		

		for(int i = 0; i < numParticlesLight; i++){
			//mass
			particles[i].m = genRand(drand48(), massLightMin, massLightMax);
			//temp
			mas[i] = particles[i].m;
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));

			//temp
			pos[i].x = particles[i].p.x;
			pos[i].y = particles[i].p.y;
			pos[i].z = particles[i].p.z;
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
			//temp
			mas[i] = particles[i].m;
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));
			
			//temp
			pos[i].x = particles[i].p.x;
			pos[i].y = particles[i].p.y;
			pos[i].z = particles[i].p.z;
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
			//temp
			mas[i] = particles[i].m;
			//position
			particles[i].p.x = genRand(drand48(), (double)0, (double)(width - 1));
			particles[i].p.y = genRand(drand48(), (double)0, (double)(height - 1));
			particles[i].p.z = genRand(drand48(), (double)0, (double)(zDepth - 1));
			
			//temp
			pos[i].x = particles[i].p.x;
			pos[i].y = particles[i].p.y;
			pos[i].z = particles[i].p.z;
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

	}
	
	printf("total particles  and rank %d, %d\n", totalParticles, my_rank);

	MPI_Bcast(mas, totalParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(pos, totalParticles, mpi_vec3, 0, MPI_COMM_WORLD);
	//frame generation (per step)
	for(int frame = 0; frame < numSteps; frame++){
		
		if(my_rank == 0){
			//reset image pixels
			for(int i = 0; i < (width * height * 3); i++){
				image[i] = 0;
			}
		}
		//sub step calculations
		for(int step = 0; step < subSteps; step++){
			
			//reset acceleration for particles
			for(int i = 0; i < totalParticles; i++){
				accel[i].x = 0;
				accel[i].y = 0;
				accel[i].z = 0;
			}

			MPI_Barrier(MPI_COMM_WORLD);
			
			MPI_Scatter(accel, blockSize, mpi_vec3, accelBlock, blockSize, mpi_vec3, 0, MPI_COMM_WORLD);

			//printf("Rank %d has %f\n", my_rank, accelBlock[0].x);

			start = MPI_Wtime();
			computeAccel(pos, mas, accelBlock, blockSize);
			end = MPI_Wtime();
			timing[t] = end - start;
			t++;

			//printf("COMPUTES rank %d\n", my_rank);

			MPI_Gather(accelBlock, blockSize, mpi_vec3, accel, blockSize, mpi_vec3, 0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);

			if(my_rank == 0){
				//update position & velocity for each particle based on acceleration
				for(int i = 0; i < totalParticles; i++){
					int pm = particles[i].m;

					particles[i].v.x += accel[i].x * timeSubStep;
					particles[i].v.y += accel[i].y * timeSubStep;
					particles[i].v.z += accel[i].z * timeSubStep;
					//clamp (restrict) velocity values to given range
					if(pm >= 11 && pm <= 15){
						clamp(particles[i].v.x, velocityHeavyMin, velocityHeavyMax);
						clamp(particles[i].v.y, velocityHeavyMin, velocityHeavyMax);
						clamp(particles[i].v.z, velocityHeavyMin, velocityHeavyMax);
					}else if(pm >= 6 && pm <= 10){
						clamp(particles[i].v.x, velocityMediumMin, velocityMediumMax);
						clamp(particles[i].v.y, velocityMediumMin, velocityMediumMax);
						clamp(particles[i].v.z, velocityMediumMin, velocityMediumMax);
					}else if(pm >= 1 && pm <= 5){
						clamp(particles[i].v.x, velocityLightMin, velocityLightMax);
						clamp(particles[i].v.y, velocityLightMin, velocityLightMax);
						clamp(particles[i].v.z, velocityLightMin, velocityLightMax);
					}

					particles[i].p.x += particles[i].v.x * timeSubStep;
					particles[i].p.y += particles[i].v.y * timeSubStep;
					particles[i].p.z += particles[i].v.z * timeSubStep;
					//handle boundary position conditions
					if(particles[i].p.x < 0 || particles[i].p.x > width){
						particles[i].v.x *= -1;
					}else if(particles[i].p.y < 0 || particles[i].p.y > height){
						particles[i].v.y *= -1;
					}else if(particles[i].p.z < 0 || particles[i].p.z > zDepth){
						particles[i].v.z *= -1;
					}
				}
			}
		}
		if(my_rank == 0){
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
			//printf("%s", filename.c_str());
			saveBMP(filename.c_str(), image, width, height);
		}
	}
	
/*		
	//all other nodes do this
	else{

	}
	*/
	if(my_rank == 0){
		//timing data processing
		double min = timing[0];
		double max = timing[0];
		double sum = 0.0;
		double avg;
		int total = numSteps * subSteps;
		for(int i = 0; i < total; i++){
			if(timing[i] < min){
				min = timing[i];
			}else if(timing[i] > max){
				max = timing[i];
			}
			sum += timing[i];
		}
		avg = sum / (double)total;
		printf("%f %f %f\n", min, max, avg);
	}
	
//	free(image);
	printf("Done\n");
	MPI_Finalize();
	return 0;
}
