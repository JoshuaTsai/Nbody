#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222
#define G (6.673e-11)
#define DataType float

DataType * pos_x;
DataType * pos_y;
DataType * vel_x;
DataType * vel_y;
DataType * mass;

using namespace std;

DataType ForceCalc(int x1, int y1, int x2, int y2, int m1, int m2, int xory){
	DataType total_force = 0;
	int dx = x1-x2;
	int dy = y1-y2;
	int dist = sqrt(dx*dx+dy*dy);
	int dist3 = dist*dist*dist;

	/* COMPUTE FORCE */
	if(xory == 0){
		total_force = -G*m1*m2/dist3*dx;
	} else if (xory == 1){
		total_force = -G*m1*m2/dist3*dy;
	}
	return total_force;
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

	int numParticlesTotal = numParticlesLight + numParticleMedium + numParticleHeavy;

	int numSteps = 0;
	int subSteps = 0;
	DataType timeSubStep;
	int particles_torecv, particle_perproc = numParticlesTotal/p, particle_left = numParticlesTotal%p;
	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){

		pos_x = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		pos_y = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		vel_x = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		vel_y = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		mass = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);

		for(int i = 0; i < numParticlesTotal; i++){
			if(i < numParticlesLight){
				mass[i] = 1;
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			} else if(i >= numParticlesLight && i < (numParticlesLight+numParticleMedium)){
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

			DataType* local_mass = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* loc_arr_src_x = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* loc_arr_dst_x = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* loc_arr_src_y = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* loc_arr_dst_y = (DataType *) malloc(sizeof(DataType) * particles_torecv);

			DataType* temp_mass = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* temp_arr_src_x = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* temp_arr_dst_x = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* temp_arr_src_y = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* temp_arr_dst_y = (DataType *) malloc(sizeof(DataType) * particles_torecv);


			DataType* compute_pos_x = (DataType *) malloc(sizeof(DataType) * particles_torecv);
			DataType* compute_pos_y = (DataType *) malloc(sizeof(DataType) * particles_torecv);
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
