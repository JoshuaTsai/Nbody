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
	if(xory == 0){ //x == 0
		total_force = -G*m1*m2/dist3*dx;
	} else if (xory == 1){
		total_force = -G*m1*m2/dist3*dy;
	}
	return total_force;
}

int main(int argc, char* argv[])
{

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
	int particle_perproc = numParticlesTotal/p;
	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){

		//------------------------------------------ Init
		pos_x = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		pos_y = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		vel_x = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		vel_y = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);
		mass = (DataType *) malloc(sizeof(DataType) * numParticlesTotal);

		for(int i = 0; i < numParticlesTotal; i++){
			if(i < numParticlesLight){
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			} else if(i >= numParticlesLight && i < (numParticlesLight+numParticleMedium)){
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			} else{
				pos_x[i] = drand48()*width;
				pos_y[i] = drand48()*height;
				vel_x[i] = drand48();
				vel_y[i] = drand48();
				numParticlesTotal--;
			}
		}

		DataType* compute_pos_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* compute_pos_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* compute_vel_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* compute_vel_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		//DataType* compute_mass = (DataType *) malloc(sizeof(DataType) * particle_perproc);

		for (int k = 1 ; k < p +1 ; k+=p ) { //distribute data
			for(int j = 0; j < particle_perproc ; j+=p){
				compute_pos_x[j] = pos_x[j*(k-1)];
				compute_pos_y[j] = pos_y[j*(k-1)];
			}
			MPI_Send(&(compute_pos_x[0]), particle_perproc, MPI_INT, k, 0, MPI_COMM_WORLD);
			MPI_Send(&(compute_pos_y[0]), particle_perproc, MPI_INT, k, 0, MPI_COMM_WORLD);
			MPI_Send(&(compute_vel_x[0]), particle_perproc, MPI_INT, k, 0, MPI_COMM_WORLD);
			MPI_Send(&(compute_vel_y[0]), particle_perproc, MPI_INT, k, 0, MPI_COMM_WORLD);
		}

		free(pos_x);
		free(pos_y);
		free(vel_x);
		free(vel_y);

		free(compute_pos_x);
		free(compute_pos_y);
		free(compute_vel_x);
		free(compute_vel_y);
	}

	else{ //all other nodes do this

		MPI_Status* stat[5];
		DataType* loc_mass = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_pos_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_pos_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_vel_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_vel_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_force_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* loc_force_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);

		DataType* temp_mass = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* temp_pos_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* temp_pos_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* temp_force_x = (DataType *) malloc(sizeof(DataType) * particle_perproc);
		DataType* temp_force_y = (DataType *) malloc(sizeof(DataType) * particle_perproc);

		MPI_Recv(&(temp_mass[0]),  particle_perproc, MPI_INT, 0, 0, MPI_COMM_WORLD, stat[0]);
		MPI_Recv(&(temp_pos_x[0]),  particle_perproc, MPI_INT, 0, 0, MPI_COMM_WORLD, stat[1]);
		MPI_Recv(&(temp_pos_x[0]),  particle_perproc, MPI_INT, 0, 0, MPI_COMM_WORLD, stat[2]);
		MPI_Recv(&(temp_force_x[0]),  particle_perproc, MPI_INT, 0, 0, MPI_COMM_WORLD, stat[3]);
		MPI_Recv(&(temp_force_y[0]),  particle_perproc, MPI_INT, 0, 0, MPI_COMM_WORLD, stat[4]);

		for (int frames = 0; frames < numSteps ; frames++){
			for (int sub = 0; sub < subSteps ; sub++){
				for(int ringcomm = 0; ringcomm < p - 1; ringcomm++){
					if (my_rank != 0){
						int rcv_frm = (my_rank)-1 > 0 ? (my_rank)-1 : p-1;
						int send_to = (my_rank%(p-1))+1;

						//compute FOrcesssssssss and update local_arrays
					}
				}
				if(my_rank != 0){
				//update position and velocity
				//MPI barrier
				}
			}
			if(my_rank != 0){
				//everyone sends to master
			}
			else{
				//master receives all and combines and then calls save BMP
				//saveBMP(argv[9], image, width, height); //almost done, just save the image
			}
			
		}
	}

	free(image);
	MPI_Finalize();
	return 0;
}
