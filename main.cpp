vel_y#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

#define g_constant (6.673e-11)

double * pos_x;
double * pos_y;
double * vel_x;
double * vel_y;
double * mass;


double ForceCalc(int x1, int y1, int x2, int y2, int m1, int m2, int xory){
	double total_force = 0;
	int dx = x1-x2;
	int dy = y1-y2;
	int dist = sqrt(dx*dx+dy*dy);
	int 3dist = dist*dist*dist;

	/* COMPUTE FORCE */
	if(xory == 0){
		total_force = -g*w1*w1/3dist*dx;
	} else if (xory == 1){
		total_force = -g*w1*w1/3dist*dy;
	}
	return total_force;
}
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
	int numParticlesLight = stoi(argv[1]);
	int numParticleMedium = stoi(argv[2]);
	int numParticleHeavy = stoi(argv[3]);
	int width = stoi(argv[7]);
	int height = stoi(argv[8]);

	int totalParticles = lightParticles + mediumParticles + heavyParticles;

	int numSteps = 0;
	int subSteps = 0;
	double timeSubStep;
	int particles_torecv, particle_perproc = numParticlesTotal/p, particle_left = numParticlesTotal%p;
	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){

		pos_x = (double *) malloc(sizeof(double) * numParticlesTotal);
		pos_y = (double *) malloc(sizeof(double) * numParticlesTotal);
		vel_x = (double *) malloc(sizeof(double) * numParticlesTotal);
		vel_y = (double *) malloc(sizeof(double) * numParticlesTotal);
		mass = (double *) malloc(sizeof(double) * numParticlesTotal);

		for(i = 0; i < numParticlesTotal; i++){
			if(i < numParticlesLight){
				mass[i] = 1;
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			} else if(i >= numParticlesLight && i < (numParticlesLight+numParticlesMedium)){
				mass[i] = 2;
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			} else{
				mass[i] = 3;
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			}
		}


		for (int outp = 0; outp < p; outp++){

		//particles for each processor
			if (outp < particle_left) {
				particles_torecv = particle_perproc + 1;
			}
			else {
				particles_torecv = particle_perproc;
			}

			local_mass = (int *) malloc(sizeof(int) * particles_torecv);
			loc_arr_src_x = (int *) malloc(sizeof(int) * particles_torecv);
			loc_arr_dst_x = (double *) malloc(sizeof(double) * particles_torecv);
			loc_arr_src_y = (int *) malloc(sizeof(int) * particles_torecv);
			loc_arr_dst_y = (double *) malloc(sizeof(double) * particles_torecv);
			loc_arr_pointer = (int *) malloc(sizeof(int) * particles_torecv);

			temp_mass = (int *) malloc(sizeof(int) * particles_torecv);
			temp_arr_src_x = (int *) malloc(sizeof(int) * particles_torecv);
			temp_arr_dst_x = (double *) malloc(sizeof(double) * particles_torecv);
			temp_arr_src_y = (int *) malloc(sizeof(int) * particles_torecv);
			temp_arr_dst_y = (double *) malloc(sizeof(double) * particles_torecv);
			temp_arr_pointer = (int *) malloc(sizeof(int) * particles_torecv);


			compute_pos_x = (int *) malloc(sizeof(int) * particles_torecv);
			compute_pos_y = (int *) malloc(sizeof(int) * particles_torecv);



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
