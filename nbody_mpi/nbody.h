#ifndef NBODY_H
#define NBODY_H

struct node;
struct particle;

/*
  This structure holds information for a single particle,
  including position, velocity, and mass.
*/
typedef struct particle{
  double x_pos, y_pos;
  double x_vel, y_vel;
  double x_force, y_force;
  double mass;
  struct node* node;
} particle_t;

typedef struct node{
  struct node *parent;
  struct node *children;
  particle_t *particle;
  int n_particles; //number of particles in this node and its sub-nodes
  double mass; // mass of the node (ie. sum of its particles mass)
  double x_center, y_center; // center of the mass
  int depth;
  int owner;
  double x_min, x_max;
  double y_min, y_max;
} node_t;


extern int nparticles;		/* number of particles to simulate */

//#define DUMP_RESULT 1
//#define DRAW_BOXES 1

#define DISPLAY_SIZE       512      /* pixel size of display window */
#define SCALE               0.03    /* sets the magnification at the origin */
                                    /* smaller #'s zoom in */
#define XMIN (-1/SCALE)
#define XMAX (1/SCALE)
#define YMIN (-1/SCALE)
#define YMAX (1/SCALE)

#define DISPLAY_RANGE       20      /* display range of fish space */
#define STEPS_PER_DISPLAY   10      /* time steps between display of fish */
#define GRAV_CONSTANT       0.01    /* proportionality constant of
                                       gravitational interaction */

#define POS_TO_SCREEN(pos)   ((int) ((pos/SCALE + DISPLAY_SIZE)/2))

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))  /* utility function */


#endif
