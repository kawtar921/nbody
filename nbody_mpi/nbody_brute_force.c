/*
** nbody_brute_force.c - nbody simulation using the brute-force algorithm (O(n*n))
**
**/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>

#ifdef DISPLAY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#include "ui.h"
#include "nbody.h"
#include "nbody_tools.h"

FILE* f_out=NULL;

int nparticles=10;      /* number of particles */
float T_FINAL=1.0;     /* simulation end time */
particle_t*particles;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

int nbre_proc = 3;

/*MPI*/
int rank , size, tag = 44;
MPI_Status status;
double* positions_part;
double* masses_part;
MPI_Datatype positionstype;
MPI_Datatype massestype;
double* max_global;

void init() {
  
    
}

void fragmentation_mpi(int nbre_proc, int nbre_part, int* part_proc,int* count_part)
{
    int part_particules = nbre_part/nbre_proc;
    int reste_particules = nbre_part%nbre_proc;
   
    int i;
    part_proc[0] = 0;
    for(i = 0 ; i < size; i++)
    {
        
        count_part[i] = part_particules + (i < reste_particules ? 1 : 0);
        part_proc[i+1] = part_proc[i] + count_part[i];
        
    }
}

#ifdef DISPLAY
Display *theDisplay;  /* These three variables are required to open the */
GC theGC;             /* particle plotting window.  They are externally */
Window theMain;       /* declared in ui.h but are also required here.   */
#endif

/* compute the force that a particle with position (x_pos, y_pos) and mass 'mass'
 * applies to particle p
 */
void compute_force(particle_t*p, double x_pos, double y_pos, double mass) {
  double x_sep, y_sep, dist_sq, grav_base;

  x_sep = x_pos - p->x_pos;
  y_sep = y_pos - p->y_pos;
  dist_sq = MAX((x_sep*x_sep) + (y_sep*y_sep), 0.01);

  /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
  grav_base = GRAV_CONSTANT*(p->mass)*(mass)/dist_sq;

  p->x_force += grav_base*x_sep;
  p->y_force += grav_base*y_sep;
}

void compute_local_force(int i, int j, int diff, double*positions, double*pos_local, double*masses, double* local_force)
{
    //compute local force
    double x_sep, y_sep, dist_sq, grav_base;
    
    x_sep = positions[2*j+0] - pos_local[2*diff+0];
    y_sep = positions[2*j+1] - pos_local[2*diff+1];
    dist_sq = MAX((x_sep*x_sep) + (y_sep*y_sep), 0.01);
    
    grav_base = GRAV_CONSTANT*(masses[i])*(masses[j])/dist_sq;

    local_force[2*diff+0] += grav_base*x_sep;
    local_force[2*diff+1] += grav_base*y_sep;
}

/* compute the new position/velocity */
void move_particle(double step, double* pos_local, double* vel_local, double masse, double* force) {

    
    *pos_local += (*vel_local)*step;
    *(pos_local+1) += (*(vel_local+1))*step;
    double x_acc = *force/(masse);
    double y_acc = *(force+1)/(masse);
    *vel_local += x_acc*step;
    *(vel_local+1) += y_acc*step;
    
    //compute statistics
    double cur_acc = (x_acc*x_acc + y_acc*y_acc);
    cur_acc = sqrt(cur_acc);
    double speed_sq = (*vel_local)*(*vel_local) + (*(vel_local+1))*(*(vel_local+1));
    double cur_speed = sqrt(speed_sq);
    
    sum_speed_sq += speed_sq;
    max_acc = MAX(max_acc, cur_acc);
    max_speed = MAX(max_speed, cur_speed);
    
}

void calcul_max(double*buffer)
{
    double* max = malloc(2*sizeof(double));
    max[0]=0.0;
    max[1]=0.0;
    int i;
    for(i=0;i<size;i++)
    {
        max[0] = buffer[2*i] > max[0] ? buffer[2*i] : max[0];
        
    }
    
    for(i=0;i<size;i++)
    {
        max[1] = buffer[2*i+1] > max[1] ? buffer[2*i+1] : max[1];
    }
    
    max_global[0]=max[0];
    max_global[1]=max[1];
}


/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step,double* local_force, int start, int end, double* positions, double* pos_local, double* masses, double* vel_local, int*part_proc, int*count_part)
{
  /* First calculate force for particles. */
  
    /*compute force*/
    int i,j;
    int nlocal = end - start;
        for(i=0; i<nlocal;i++)
    {
        local_force[2*i+0] = 0;
        local_force[2*i+1] = 0;
    }
    for(i = start; i < end; i++)
    {
        int diff = i - start;
        for(j = 0; j < nparticles; j++)
        {
            if(i!=j)
                compute_local_force(i, j, diff, positions, pos_local,masses,local_force);
        }
    }
    
    /*compute positions*/
    for(i=0; i<nlocal; i++) {
        move_particle(step,pos_local+2*i, vel_local+2*i, *(masses+start+i),local_force+2*i);
    }
    
    /*send all positions to the combined data*/
    MPI_Allgatherv(pos_local, 2*nlocal,MPI_DOUBLE,positions_part,count_part, part_proc, positionstype,MPI_COMM_WORLD);

}

/* display all the particles */
void draw_all_particles() {
  int i;
  for(i=0; i<nparticles; i++) {
    int x = POS_TO_SCREEN(particles[i].x_pos);
    int y = POS_TO_SCREEN(particles[i].y_pos);
    draw_point (x,y);
  }
}

void run_simulation(double* local_force, int start, int end, double* positions, double* pos_local, double* masses, double* vel_local,int*part_proc, int*count_part) {
    double t = 0.0, dt = 0.01;
    
    while (t < T_FINAL && nparticles>0) {
        
    /* Update time. */
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt,local_force,start, end, positions, pos_local, masses, vel_local,part_proc,count_part);
    
        double* max_send = malloc(2*sizeof(double));
        max_send[0]=max_speed;
        max_send[1]=max_acc;
        MPI_Send(max_send,2,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
    
    if(rank==0)
    {
        int proc;
        double* recu = malloc(2*size*sizeof(double));
        for(proc = 0;proc < size; proc++)
        {
            MPI_Recv(recu+2*proc,2,MPI_DOUBLE,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
        }
        
        calcul_max(recu);
    }
    
    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */
        
    MPI_Bcast(max_global,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    dt = 0.1*max_global[0]/max_global[1];
        
    
    /* Plot the movement of the particle */
#if DISPLAY
    clear_display();
    draw_all_particles();
    flush_display();
#endif
  }
    

}

/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char**argv)
{
  if(argc >= 2) {
    nparticles = atoi(argv[1]);
  }
  if(argc == 3) {
    T_FINAL = atof(argv[2]);
  }

  /*Variables locales relatives à chaque proc*/
    int local_nparticles;//nombre de particules à la charge du proc
    double* local_pos;//Tableau des positions des particules
    double* local_vel;//Tableau des vitesses des particules
    double* local_force;//Tableau des forces
    int* part_proc;
    int* count_part;
    
  //init mpi
    if(MPI_Init( &argc, &argv ) != MPI_SUCCESS)
    {
        printf("Can't initialize MPI");
        return 1;
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_Type_vector(1, 2, 1, MPI_DOUBLE, &positionstype);
    MPI_Type_commit(&positionstype);
    
    MPI_Type_vector(1, 2, 1, MPI_DOUBLE, &massestype);
    MPI_Type_commit(&massestype);
    
    part_proc = malloc((size+1)*sizeof(int));
    count_part = malloc(size*sizeof(int));
    
  /* Allocate global shared arrays for the particles data set. */
  particles = malloc(sizeof(particle_t)*nparticles);
  positions_part = malloc(sizeof(double)*2*nparticles);
  masses_part = malloc(sizeof(double)*nparticles);
  if(rank==0)
  {
      all_init_particles(nparticles,positions_part, masses_part);
   
  }
    
  /*Broadcast informations to all processes*/
    MPI_Bcast(positions_part,2*nparticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(masses_part,nparticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    
    
  /*Partager les particules entre les processus*/
    fragmentation_mpi(size,nparticles,part_proc,count_part);
    local_nparticles = count_part[rank];
    
    /*Structures de données locales*/
    local_pos = malloc(sizeof(double)*2*local_nparticles);
    local_vel = malloc(sizeof(double)*2*local_nparticles);
    local_force = malloc(sizeof(double)*2*local_nparticles);
    memcpy(local_pos, positions_part+2*part_proc[rank], 2*local_nparticles*sizeof(double));
    local_init_particles_v(local_nparticles,local_vel,local_pos);
    max_global = malloc(sizeof(double)*2);
    max_global[0]=0.0;
    max_global[1]=0.0;
    
    
  /* Initialize thread data structures */
#ifdef DISPLAY
  /* Open an X window to display the particles */
  simple_init (100,100,DISPLAY_SIZE, DISPLAY_SIZE);
#endif

  struct timeval t1, t2;
  gettimeofday(&t1, NULL);
   
  /* Main thread starts simulation ... */
  run_simulation(local_force,part_proc[rank], part_proc[rank + 1], positions_part, local_pos, masses_part,local_vel,part_proc,count_part);
    
  gettimeofday(&t2, NULL);

  double duration = (t2.tv_sec -t1.tv_sec)+((t2.tv_usec-t1.tv_usec)/1e6);

#ifdef DUMP_RESULT
  FILE* f_out = fopen("particles.log", "w");
  assert(f_out);
  print_particles(f_out, root);
  fclose(f_out);
#endif

    
  printf("--------------rank %d---------------\n",rank);
  printf("nparticles: %d\n", nparticles);
  printf("T_FINAL: %f\n", T_FINAL);
  printf("-----------------------------\n");
  printf("Simulation took %lf s to complete\n", duration);
  MPI_Finalize();

#ifdef DISPLAY
  clear_display();
  draw_all_particles();
  flush_display();

  printf("Hit return to close the window.");

  getchar();
  /* Close the X window used to display the particles */
  XCloseDisplay(theDisplay);
#endif
  
    free(part_proc);
    free(count_part);
    free(positions_part);
    free(masses_part);
    free(local_pos);
    free(local_vel);
    free(local_force);
  return 0;
}
