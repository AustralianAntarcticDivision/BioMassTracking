#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Utils.h>

/**************************************************************************
 *
 * This is a library of functions used in the R package "BiomassTracking",
 * developed by Thorsten Lenser as part of his final project for the course
 * "MRes in Mathematics in the Living Environment while he was a guest
 * of the Australian Antarctic Division in Kingston, Tasmania, in 
 * May - August 2005.
 *
 * Author: Thorsten Lenser (thorsten.lenser@gmx.net)
 *
 * Some notes:
 * - Some of the source is is quite hard to read, because the three-
 * dimensional tensors coming from the R scripts had to be translated
 * into 1-D arrays. This is indicated by comments.
 * - Although greatest care has been taken in developing this code, it was
 * written in a relatively tight timeframe because it formed part of a MRes
 * project. I therefore advise the user not to take any results of this code
 * for granted and at least to compare it with observational data. If any 
 * problems arise, please contact me under the email address given above.
 *
 **************************************************************************/

/* Functions called from the R scripts */

// Carry out the particle tracking for the estimation
void particle_tracking(double*, double*, int*, int*, int*, int*,
                       double*, int*, int*, double*, double*, 
                       double*, int*, int*, int*, double*,
                       int*, double*, int*);
                    
// Particle tracking to estimate the movement of "native" biomass
void particle_tracking_native_biomass(double*,double*,int*,int*,int*,int*,
                       double*,int*,int*,double*,double*,double*,
                       int*,int*,double*);
 
// Carry out the particle tracking for comparison
void particle_tracking_compare_inflow(double*,double*,int*,int*,int*, int*,double*, int*,
                       double*, double*, double*, int*,int*,int*, double*);
                       
// Comparison with a given start setup of biomass distribution
void particle_tracking_compare_start_setup(double*,double*,int*,int*,int*, 
                       int*,double*,int*,double*,double*,double*,double*,
                       int*,double*);                       
                       
// Estimate a movement tensor
void estimate_P(double*, int*, int*, int*, double*, int*, int*, int*, double*);

// Calculate the biomass distribution from the movement kernel
void biomass_tracking_inflow(double*, double*, double*, int*, int*, int*,
                       int*, double*, int*, double*, int*);

// Calculate the biomass distributions starting from a given setup
void biomass_tracking_start_setup(double*, double*, double*, 
                       int*, int*, int*, int*, double*, double*);
                       
                       
                       
/* Local helper functions */
                       
// Return the local neighbourhood number of a polygon relative to another polygon
int local_neighbour(int, int, int*, int*, int);

// One advection time step using a 4th order Runge-Kutta scheme
void runge_kutta(double*, double*, double*, double*, int*, int*, double*);

// Get the position of a single particle
int get_single_position(double, double, int*, int, int*, int*);

// Get the polygon position of all particles
void get_position(int*, double*, double*, int*, int*, int, int*, int*);

// Move the particles for one time step
void move_particles(double*,double*,double*,double*,int*,int*,int*,int*,int,double*,
                    int,int,double*,double,double);
                    
// Find the infl probabilities into the boundary cells
double* get_infl_prob(double*, double*, int*, int*, int*);
     

                       

/*****************************************************************************
 *
 * "particle_tracking"
 *
 * The implementation of a simple 4th-order Runge-Kutta 2-D integration
 * scheme to follow the paths of particles through a given flow field.
 * After a phase to settle down, the moves of the particles from one
 * polygon to another are recorded and reported back to the calling R function.
 * These can then be used for statistics to estimate movement characteristics
 * in the given flow field.
 * 
 * Note that the integration scheme uses a fixed time step which leaves it to 
 * clumsy for serious numerical integration. However, it is meant to be used 
 * for demonstration purposes only and can easily be extended to a more precise
 * algorithm when needed.
 *
 * Be careful: R converts matrices to arrays on a "by column" basis, i.e. 
 * different than the notation I used in some R scripts!
 *
 * rnorm is a vector of length n*m +1, containing random numbers drawn 
 * according to a normal distribution with mean = 0 and sd = ?
 *
 * If logyesno == 1, the particle positions will be logged in 'log',
 * which holds successively the x and y position of each particle
 * in each time step. This should only be used for a small number of
 * particles.
 *
 *****************************************************************************/ 

void particle_tracking(double *U,double *V,int *n,int *m,int *num_part, int *S,
                       double *t_step, int *condition_window, int *end_t_counter,
                       double *x, double *y, double *result_t, int *result_S1, 
                       int* result_S2, int *result_S3, double *diff_sd,
                       int* logyesno, double* log, int* from_boundary){

	int t,i;
	int k;
	int pos[*num_part];
	int max_S = 0;
	int num_poly;
	int a,b;
	double rand_no;
	int left_end, right_end;
	double mean_abs = 0;// For the random diffusion.
	int not_into_right_poly = 1; // Boolean
	int continue_estimation = 1;
	int* num_moves = (int*) R_alloc(*condition_window, sizeof(int));
	for(i=0; i < *condition_window; i++){
        num_moves[i] = 0;
    }
	int current_moves = 0;
	int sum_moves = 0;
	
	// Create an array of counters for the number of moves a particle has made
	int* move_counter = (int*) R_alloc(*num_part, sizeof(int));
	for(i=0; i < *num_part; i++){
        move_counter[i] = 0;
    }
    int move_tmp = 0;

	for(i=0; i < (*n)*(*m); i++){
		if(S[i] > max_S){
			max_S = S[i];
		}
		mean_abs += fabs(U[i]);
		mean_abs += fabs(V[i]);
	}
	mean_abs /= 2* *n * *m;
	
	num_poly = max_S +1;
	
	// Start Random Number Generator
	GetRNGstate();
	
	// Get inflow probabilities from outside the boundary
	double* infl_prob = get_infl_prob(U,V,S,m,n);
	
    // Let the particles enter from the outside
    if(*from_boundary == 1){
    for (i=0; i < *num_part; i++){
        // Make sure particle is inserted into polygon 1
        not_into_right_poly = 1;
        while(not_into_right_poly == 1){               
            rand_no = unif_rand();
            // Search for entry point using divide-and-conquer strategy
            b = 1;
            left_end = 0;
            right_end = (2*(*m) + 2*(*n))-1;
            while(b==1){
                k = (right_end + left_end)/2;
                if ((k > 0) && (infl_prob[k-1] >= rand_no)){
                    right_end = k;
                }
                else if((infl_prob[k] <= rand_no) && (k < ((2*(*m) * 2*(*n))-1))){
                    left_end = k+1;
                }
                else b = 0;
                if(k == ((2*(*m) * 2*(*n))-1)) b = 0;
                if(k==0) b = 0;
                if(left_end == right_end) b = 0;
            }
              
            // Get particle position at entry point
            if(k < (*m)){
                x[i] = 0.5001;
                y[i] = k +1 + 0.999 * unif_rand() - 0.499;
            }
            else if(k < 2*(*m)){
                x[i] = *n + 0.4999;
                y[i] = k - (*m) +1 + 0.999 * unif_rand() - 0.499;
            }
            else if(k < (2*(*m) + (*n))){
                x[i] = k - 2*(*m) +1 + 0.999 * unif_rand() - 0.499;
                y[i] = 0.5001;
            }
            else{
                x[i] = k - 2*(*m) - (*n) +1 + 0.999 * unif_rand() - 0.499;
                y[i] = *m + 0.4999;
            }
            if(get_single_position(x[i],y[i],S,num_poly,m,n) != 0) not_into_right_poly = 0; 
        }
    }
    } // if(*from_boundary == 1)
    
    // Get current position of particles when move-recording begins
    // All particles came from the outside, so origin = num_poly
    get_position(pos,x,y,num_part,S,num_poly,m,n);
    for(i=0; i < *num_part; i++){
        result_S2[i] = pos[i];
        result_S1[i] = num_poly;
    }
                   
    t = 1;
    // Estimate until the number of moves in the last condition_window time steps is 0.
    while((continue_estimation == 1) && ((*end_t_counter == 0) || (t <= *end_t_counter))){
        current_moves = 0;
              
        // Check whether an interruption is requested by R
        R_CheckUserInterrupt();
        	
	 	for(i=0;i < *num_part; i++){
            if(*logyesno == 1){ // Log the particle tracks
                log[(t-1)*2*(*num_part) + 2*i] = x[i];
                log[(t-1)*2*(*num_part) + 2*i +1] = y[i];
            }
                  
            if(result_S3[(num_poly-1) * (*num_part) +i] == 0){      // Particle has not already done its work
                // Advection
                runge_kutta(&(x[i]), &(y[i]), U, V, m, n, t_step);
                // Diffusion
               	x[i] = x[i] + (mean_abs * *diff_sd * norm_rand() * *t_step);
               	y[i] = y[i] + (mean_abs * *diff_sd * norm_rand() * *t_step);
    			
			    // Give a particle position -10,-10 if it has left the domain
                // Then get position in S and see if particle has moved
                if((x[i] >= (*n + 0.5)) || (x[i] <= 0.5) 
    			|| (y[i] >= (*m + 0.5)) || (y[i] <= 0.5)){
    				x[i] = -10;
    				y[i] = -10;
                                
    				pos[i] = num_poly;
    			}
    			else{
    				a = (int) fround(x[i],0);
    				b = (int) fround(y[i],0);
    				pos[i] = S[(a-1)*(*m) + b-1];
    			}
    			
    			// Has the particle hit land?
    			if((pos[i] == 0) && (x[i] != -10)){
                    // Move it outside
                    x[i] = -10;
                    y[i] = -10;
                    pos[i] = num_poly;
                }
                         
                // Record the move if the particle has moved
                for(move_tmp=0; move_tmp < num_poly; move_tmp++){
                    if(result_S3[move_tmp * (*num_part) + i] == 0){            
                        if(result_S2[move_tmp * (*num_part) + i] != pos[i]){
                            current_moves++;
                            if(move_tmp < (num_poly-1)) result_S1[(move_tmp+1) * (*num_part) + i] = result_S2[move_tmp * (*num_part) + i];
                            result_t[move_tmp * (*num_part) +i] = t * (*t_step);
                            result_S3[move_tmp * (*num_part) +i] = pos[i];
                            if(move_tmp < (num_poly-1)) result_S2[(move_tmp+1) * (*num_part) + i] = pos[i];
                        }
                        break;
                    }
                }
                				              
                // Insert the particles again which have left the domain
                if((pos[i] == num_poly) && (result_S3[(num_poly-1) * (*num_part) +i] == 0)){
                    not_into_right_poly = 1;
                    while(not_into_right_poly == 1){
                        rand_no = unif_rand();
                        // Search for entry point using divide-and-conquer strategy
                        b = 1;
                        left_end = 0;
                        right_end = (2*(*m) + 2*(*n))-1;
                        while(b==1){
                            k = (right_end + left_end)/2;
                            if ((k > 0) && (infl_prob[k-1] >= rand_no)){
                                right_end = k;
                            }
                            else if((infl_prob[k] <= rand_no) && (k < ((2*(*m) * 2*(*n))-1))){
                                left_end = k+1;
                            }
                            else b = 0;
                            if(k == ((2*(*m) * 2*(*n))-1)) b = 0;
                            if(k==0) b = 0;
                            if(left_end == right_end) b = 0;
                        }
                        
                        // Get particle position at entry point, 0 <= k < 2*m + 2*n
                        if(k < (*m)){
                             x[i] = 0.501;
                             y[i] = k + 1 + 0.999*unif_rand() - 0.499;
                        }
                        else if(k < 2*(*m)){
                             x[i] = *n + 0.499;
                             y[i] = k - (*m) +1 + 0.999*unif_rand() - 0.499;
                        }
                        else if(k < (2*(*m) + (*n))){
                             x[i] = k - 2*(*m) +1 + 0.999*unif_rand() - 0.499;
                             y[i] = 0.501;
                        }
                        else{
                             x[i] = k - 2*(*m) - (*n) +1 + 0.999*unif_rand() - 0.499;
                             y[i] = *m + 0.499;
                        }
                        
                        // if(t > (*end_t_counter / 2)){
                        a = (int) fround(x[i],0);
        				b = (int) fround(y[i],0);
        				for(move_tmp=1; move_tmp < num_poly; move_tmp++){
                            if(result_S3[move_tmp * (*num_part) +i] == 0){
                                result_S1[move_tmp * (*num_part) +i] = num_poly;
                                result_S2[move_tmp * (*num_part) +i] = S[(a-1)*(*m) + b-1];
                                pos[i] = S[(a-1)*(*m) + b-1];
                                break;
                            }
                        }
                        
                        // Repeat if particle has entered on land
                        if(pos[i] != 0){
                            not_into_right_poly = 0;
                        }
                    }
                    // Do not log the move across the domain
                    if(*logyesno == 1){ // Log the particle tracks
                        log[(t-1)*2*(*num_part) + 2*i] = -10;
                        log[(t-1)*2*(*num_part) + 2*i +1] = -10;
                    }
                }
            } // if(particle has not done its last move yet)
            else{
                x[i] = -10;
                y[i] = -10;
                pos[i] = num_poly;
            }
		}

		if(*condition_window > 0){
            sum_moves = 0;
            for(i = (*condition_window-1); i > 0; i--){
                num_moves[i] = num_moves[i-1];
                sum_moves += num_moves[i];
            }
            num_moves[0] = current_moves;
            sum_moves += num_moves[0];
            if((sum_moves == 0) && (t > *condition_window)){
                continue_estimation = 0;
            }
        }

        t++;
	}
	
	// Take care of all the particles which got stuck
    for(i=0; i < *num_part; i++){
        for(move_tmp=0; move_tmp < num_poly; move_tmp++){
            if(result_S3[move_tmp * (*num_part) +i] == 0){
                if(*condition_window > 0){
                    result_t[move_tmp * (*num_part) +i] = (t - *condition_window + 10) * *t_step;
                }
                else result_t[move_tmp * (*num_part) +i] = t * *t_step;
                break;
            }
        }     
    }

	// Stop Random Number Generator
	PutRNGstate();
}

/*****************************************************************************
 *
 * "particle_tracking_native_biomass"
 *
 * Simpler variant of "particle_tracking, estimates the movement of "native" 
 * biomass, i.e. biomass that is uniformly distributed over the area and 
 * makes only one move, namely out of its starting polygon. Further movement
 * is then governed by the movement description estimated in the original
 * function.
 *
 *****************************************************************************/
void particle_tracking_native_biomass(double *U,double *V,int *n,int *m,
                                      int *num_part, int *S, double *t_step, 
                                      int *condition_window, int *end_t_counter,
                                      double *x, double *y, 
                                      double *result_t,
                                      int* result_S2, int *result_S3,
                                      double *diff_sd){

	int t,i;
	int pos[*num_part];
	int max_S = 0;
	int num_poly;
	int a,b;
	double mean_abs = 0;     // For the random diffusion.
	int continue_estimation = 1;
	int* num_moves = (int*) R_alloc(*condition_window, sizeof(int));
	for(i=0; i < *condition_window; i++){
        num_moves[i] = 0;
    }
	int current_moves = 0;
	int sum_moves = 0;
	
	for(i=0; i < (*n)*(*m); i++){
		if(S[i] > max_S){
			max_S = S[i];
		}
		mean_abs += fabs(U[i]);
		mean_abs += fabs(V[i]);
	}
	mean_abs /= 2* *n * *m;
	
	num_poly = max_S +1;
	
	// Start Random Number Generator
	GetRNGstate();
	
    // Get current position of particles when move-recording begins
    get_position(pos,x,y,num_part,S,num_poly,m,n);
    for(i=0; i < *num_part; i++){
        result_S2[i] = pos[i];
    }
                   
    t = 1;
    // Estimate until the number of moves in the last no_move_window time steps is 0.
    while((continue_estimation == 1) && ((*end_t_counter == 0) || (t <= *end_t_counter))){
        current_moves = 0;
              
        // Check whether an interruption is requested by R
        R_CheckUserInterrupt();
        	
	 	for(i=0;i < *num_part; i++){
            if(result_S3[i] == 0){      // Particle has not already done its work
                // Advection
                runge_kutta(&(x[i]), &(y[i]), U, V, m, n, t_step);
                // Diffusion
               	x[i] = x[i] + (mean_abs * *diff_sd * norm_rand() * *t_step);
               	y[i] = y[i] + (mean_abs * *diff_sd * norm_rand() * *t_step);
    			
			    // Give a particle position -10,-10 if it has left the domain
                // Then get position in S and see if particle has moved
                if((x[i] >= (*n + 0.5)) || (x[i] <= 0.5) 
    			|| (y[i] >= (*m + 0.5)) || (y[i] <= 0.5)){
    				x[i] = -10;
    				y[i] = -10;
                                
    				pos[i] = num_poly;
    			}
    			else{
    				a = (int) fround(x[i],0);
    				b = (int) fround(y[i],0);
    				pos[i] = S[(a-1)*(*m) + b-1];
    			}
    			
    			// Has the particle hit land?
    			if((pos[i] == 0) && (x[i] != -10)){
                    // Move it outside
                    x[i] = -10;
                    y[i] = -10;
                    pos[i] = num_poly;
                }
                         
                // Record the move if the particle has moved
                if(result_S2[i] != pos[i]){
                    current_moves++;
                    result_t[i] = t * (*t_step);
                    result_S3[i] = pos[i];
                }  				              
            } 
		}

		if(*condition_window > 0){
            sum_moves = 0;
            for(i = (*condition_window-1); i > 0; i--){
                num_moves[i] = num_moves[i-1];
                sum_moves += num_moves[i];
            }
            num_moves[0] = current_moves;
            sum_moves += num_moves[0];
            if((sum_moves == 0)  && (t > *condition_window)){
                continue_estimation = 0;
            }
        }

        t++;
	}
	
	// Take care of all the particles which got stuck
    for(i=0; i < *num_part; i++){
        if(result_S3[i] == 0){
            result_t[i] = (t - *condition_window +10) * *t_step;
            break;
        }
    }

	// Stop Random Number Generator
	PutRNGstate();
}




/*****************************************************************************
 *
 * "particle_tracking_compare"
 *
 * This version of "particle_tracking" is used to compare the result of the 
 * original function used in "biomass_tracking" against the result of a 
 * purely particle-oriented simulation, i.e. where particles are advected 
 * and their locations are counted.
 *
 * Particles are inserted in "insert_at_boundary" steps into polygon 
 * "infl_poly", according to the inflow probabilities at the boundary of 
 * the domain.
 *
 * Be careful: R converts matrices to arrays on a "by column" basis, i.e. 
 * different than the notation I used in some R scripts!
 *
 *****************************************************************************/

void particle_tracking_compare_inflow(double *U,double *V,int *n,int *m,
                               int *num_part, int *S, double *t_step, 
                               int *end_t_counter, double *x, double *y, 
                               double *N, int *insert_at_boundary,
                               int *infl_poly, int *subdiv, double *diff_sd){

	int t,i;
	int k;
	int pos[*num_part];
	int max_S = 0;
	int num_poly;
	int b; // Used as boolean
	int left_end, right_end;
	double rand_no;
	int part_count = 0; // Count the number of particles already inserted
	int not_into_right_poly = 1; // Used as boolean variable
	double mean_abs = 0;
	double dbl_tmp;

	for(i=0; i < (*n)*(*m); i++){
		if(S[i] > max_S){
			max_S = S[i];
		}
		mean_abs += fabs(U[i]);
		mean_abs += fabs(V[i]);
	}
	mean_abs /= 2* *n * *m;
	
	num_poly = max_S +1;
	
	// Start Random Number Generator
	GetRNGstate();
	
	// Get inflow probabilities into the boundary cells
	double* infl_prob = get_infl_prob(U,V,S,m,n);
	
	// Check whether there is any inflow into the inflow polygon
    // Otherwise the particle insertion will not work.
    dbl_tmp = 0;
    for(i=0;i < *m;i++){
        if(S[i] == *infl_poly) dbl_tmp += infl_prob[i];
        if(S[(*n -1)*(*m) + i] == *infl_poly) dbl_tmp += infl_prob[*m +i];
    }
    for(i=0;i < *n;i++){
        if(S[i * (*m)] == *infl_poly) dbl_tmp += infl_prob[2*(*m) + i];
        if(S[i * (*m) + *m -1] == *infl_poly) dbl_tmp += infl_prob[2*(*m) + (*n) +i];
    }

    if(dbl_tmp <= 0){
        error("There is no inflow into the inflow polygon!\n Program exits\n");
        exit(-1);
    }
    
	for (t=1;t <= *end_t_counter;t++){
        // Check whether an interruption is requested by R
        R_CheckUserInterrupt();
        
        move_particles(x,y,U,V,m,n,&(part_count),S,num_poly,t_step,*subdiv,0,infl_prob,mean_abs,*diff_sd);
        // Insert the particles	
        if(t <= *insert_at_boundary){
            for (i=0; i < ((*num_part) / (*insert_at_boundary)); i++){
                // Make sure particle is inserted into polygon 1
                not_into_right_poly = 1;
                while(not_into_right_poly == 1){               
                    rand_no = unif_rand();
                    // Search for entry point using divide-and-conquer strategy
                    b = 1;
                    left_end = 0;
                    right_end = (2*(*m) + 2*(*n))-1;
                    while(b==1){
                        k = (right_end + left_end)/2;
                        if ((k > 0) && (infl_prob[k-1] >= rand_no)){
                            right_end = k;
                        }
                        else if((infl_prob[k] <= rand_no) && (k < ((2*(*m) * 2*(*n))-1))){
                             left_end = k+1;
                        }
                        else b = 0;
                        if(k == ((2*(*m) * 2*(*n))-1)) b = 0;
                        if(k==0) b = 0;
                        if(left_end == right_end) b = 0;
                    }
              
                    // Get particle position at entry point
                    if(k < (*m)){
                         x[part_count + i] = 0.5001;
                         //x[part_count + i] = 0.5001 + 4*unif_rand();
                         y[part_count + i] = k +1 + 0.999 * unif_rand() - 0.499;
                    }
                    else if(k < 2*(*m)){
                         x[part_count + i] = *n + 0.4999;
                         y[part_count + i] = k - (*m) +1 + 0.999 * unif_rand() - 0.499;
                    }
                    else if(k < (2*(*m) + (*n))){
                         x[part_count + i] = k - 2*(*m) +1 + 0.999 * unif_rand() - 0.499;
                         y[part_count + i] = 0.5001;
                    }
                    else{
                         x[part_count + i] = k - 2*(*m) - (*n) +1 + 0.999 * unif_rand() - 0.499;
                         y[part_count + i] = *m + 0.4999;
                    }
                    if(get_single_position(x[part_count+i],y[part_count+i],S,num_poly,m,n) == *infl_poly) not_into_right_poly = 0; 
                }
            }
            part_count += i; // i is increased before the condition yields FALSE
        } 
        
        get_position(pos,x,y,&(part_count),S,num_poly,m,n);
        
        // Estimate biomass distribution:
        for(i=0;i < part_count;i++){
       	    N[t*(num_poly)+pos[i]-1] += 1;
        }                 
        // Norm to 100 in R!
	}
	
	// Stop Random Number Generator
	PutRNGstate();
}


/***************************************************************************
 * 
 * "particle_tracking_compare_start_setup"
 *
 * This function is very similar to "particle_tracking_compare", with the
 * difference that the biomass is not insert into a specific polygon, but
 * has a given starting distribution in the polygons.
 *
 ***************************************************************************/

void particle_tracking_compare_start_setup(double *U,double *V,int *n,int *m,int *num_part, 
                       int *S, double *t_step, int *end_t_counter,
                       double *x, double *y, double *N, double *start_setup,
                       int *subdiv, double *diff_sd){

	int t,i;
	int pos[*num_part];
	int max_S = 0;
	int num_poly = 0;
	double mean_abs = 0;
	
	for(i=0; i < (*n)*(*m); i++){
		if(S[i] > max_S){
			max_S = S[i];
		}
		mean_abs += fabs(U[i]);
		mean_abs += fabs(V[i]);
	}
	mean_abs /= 2* *n * *m;
	
	num_poly = max_S +1;
	
	// Start Random Number Generator
	GetRNGstate();       
	
	// Count how many particles have been inserted per polygon
	double *dist_count = (double*) R_alloc(num_poly,sizeof(double));
	for(i=0; i < num_poly; i++){
        dist_count[i] = 0;
    }
    int pos_tmp;
    int not_in_right_poly = 1;
    int counter = 0;
    
    // Insert the particles according to the start setup
	for(i=0; i < *num_part; i++){
        counter = 0;
        while(not_in_right_poly == 1){
            counter++;
            x[i] = unif_rand() * (*n - 0.0002) + 0.5001;
            y[i] = unif_rand() * (*m - 0.0002) + 0.5001;
            pos_tmp = get_single_position(x[i],y[i],S,num_poly,m,n);
            if(dist_count[pos_tmp-1] < start_setup[pos_tmp-1]){
                not_in_right_poly = 0;
                dist_count[pos_tmp-1] += ((double)100 / *num_part);
            }
            // Stop after 1000 tries, hopefully distribution will be nearly reached by then.
            if(counter > 1000){
                i = *num_part;
                break;
            }
        }
        not_in_right_poly = 1;
    }        
            
    get_position(pos,x,y,num_part,S,num_poly,m,n);
    for(i=0;i < *num_part;i++){
        N[pos[i]-1] += 1;
    }
    
	for (t=1;t <= *end_t_counter;t++){
        // Check whether an interruption is requested by R
        R_CheckUserInterrupt();
        
        move_particles(x,y,U,V,m,n,num_part,S,num_poly,t_step,*subdiv,0,NULL,mean_abs,*diff_sd);
        
        get_position(pos,x,y,num_part,S,num_poly,m,n);
        
        // Estimate biomass distribution:
        for(i=0;i < *num_part;i++){
       	    N[t*(num_poly)+pos[i]-1] += 1;
        }                 
        // Norm to 100 in R!
	}
	
	// Stop Random Number Generator
	PutRNGstate();
}



/****************************************************************************
 *
 * "estimate_P"
 *
 * The function performs the estimation of the movement tensor after the 
 * particle tracking has been done. Its methodology is described in the 
 * report this works belongs to.
 *
 ****************************************************************************/
 
void estimate_P(double *P, int *nk, int *origin, int *dest, double *times,
                int *no_age_classes, int *no_particles, int *no_neighbours, 
                double *delta){
    // P is assumed to be a 3-D tensor of dimensions no_age_classes x 
    // no_origin x no_dest.
    // It is also assumed that only records relating to one polygon have 
    // been passed to this function, for which the estimation is returned.
    int age,i,j;
    int count, abs_count,k; 
                  // Count counts the particles fulfilling a certain condition
                  // abs_count counts all particles in a group
                  // k runs through all the records

    // Leave plus class out
    for(age = 1; age < *no_age_classes; age++){
        for(i = 0; i < *no_neighbours; i++){
            for(j = 0; j < *no_neighbours; j++){
                count = 0;
                abs_count = 0;
                for(k = 0; k < *no_particles; k++){
                    if((origin[k] == nk[i]) && (times[k] >= (age-1)*(*delta))){ 
                        abs_count++;
                        if((times[k] < age * (*delta)) && (dest[k] == nk[j])){
                            count++;
                        }
                    }
                }
                // when converting an array to a vector, 'i' is the fastest variable,
                // 'j' the second and 'age' the slowest.
                P[i +  (j * *no_neighbours) + (age-1)*(*no_neighbours * *no_neighbours)] 
                    = ((double) count) / ((double) imax2(1,abs_count));      
            }
        }
    }

    // Treat plus class differently, this class retains all biomass
    for(i = 0; i < *no_neighbours; i++){
        for(j = 0; j < *no_neighbours; j++){
            // when converting an array to a vector, 'i' is the fastest variable,
            // 'j' the second and 'age' the slowest.
            P[i +  (j * *no_neighbours) + (*no_age_classes-1)*(*no_neighbours * *no_neighbours)] 
            = 0;
        }
    }
}



/************************************************************************
 *
 * "biomass_tracking_inflow"
 * 
 * This function does the hard work for the biomass_tracking.R script,
 * calculating the flow of biomass through the domain on the basis
 * of a movement kernel.
 *
 * Since we cannot use lists here, the addressing of different elements
 * in the movement kernel is a bit more complicated than in R. It works
 * like this:
 * The first information available is the number of polygons including
 * the outside, 'num_poly' and the number of age classes. With 'num_poly',
 * the number of neighbours per polygon can be accessed, stored in 
 * the array 'num_neighbours'. With all this together, finally,
 * the movement tensors can be accessed, which are all stored together
 * in the double array 'P'. 'P' contains the tensors in their natural order,
 * with origin ('i') as the fastest variable, destination ('j') the second
 * and the age class ('age') as the slowest.
 *
 * Note that while R addresses array in 1..n style, C uses 0..n-1 style.
 * I kept this consisted, so that a lot of '-1's appear in the code below.
 * => The polygons can be thought of as labelled 0..num_poly-1.
 *
 * The results are stored in N_delta, which is called that way because it
 * gives the biomass distribution in delta-time steps. The is a vector of
 * ent_t_counter blocks of num_poly doubles, in intrinsic order.
 *
 **************************************************************************/
 
void biomass_tracking_inflow(double *N_delta, double* max_times, double *P,
                       int *num_neighbours, int *num_poly, int *nk,
                       int *num_age_classes, double *delta,
                       int *insert_time, double *insert_t_step, 
                       int *infl_poly){
                       
    int i,j,k;
    int age;
   
    double dbl_sum = 0; 

    // The inflow polygon is actually infl_poly -1 in C terminology
    int infl_poly_c = *infl_poly -1;                       
    
    // N holds the current biomass distribution. Its first index is the
    // polygon, second index is origin and third index age class.
    // It does not include the outside!
    double ***N = (double***) R_alloc(*num_poly,sizeof(double**));
    for(i = 0; i < (*num_poly -1); i++){
        N[i] = (double**) R_alloc(num_neighbours[i],sizeof(double*));
        for(j = 0; j < num_neighbours[i]; j++){
            N[i][j] = (double*) R_alloc(*num_age_classes,sizeof(double));
            for(age = 0; age < *num_age_classes; age++){
                N[i][j][age] = 0;
            }
        }
    }
    double outside = 0;
    
    // M[i][j] holds the biomass to be moved from i to j
    // Put an additional row into M containing the inflow from the outside
   	double **M = (double**) R_alloc(*num_poly,sizeof(double*));
   	for(i=0; i < *num_poly; i++){
        M[i] = (double*) R_alloc(*num_poly,sizeof(double));
        for(j=0; j < *num_poly; j++){
            M[i][j] = 0;
        }
    }

    if((*infl_poly < 1) || (*infl_poly >= *num_poly)){
        error("Wrong inflow polygon given to 'biomass_tracking1'!\n Remember that the outside is not a valid inflow polygon\n");
    }
    
    
    // Calculate the starting points of the polygons in P and nk
    int *start_P = (int*) R_alloc(*num_poly,sizeof(int));
    int *start_nk = (int*) R_alloc(*num_poly,sizeof(int));
    start_P[0] = 0;
    start_nk[0] = 0;
    for(i = 1; i < *num_poly; i++){
        start_P[i] = start_P[i-1] 
            + num_neighbours[i-1] * num_neighbours[i-1] * *num_age_classes;
        start_nk[i] = start_nk[i-1] + num_neighbours[i-1];
    }
    
     
    double delta_time = 0;
    int counter = 0;
    int insert_counter = 0;
    // Convert insert_time into the number of delta-time steps in which 
    // biomass is inserted. Note that insert_delta_time is a non-integer 
    // number, so some biomass has to be substracted at the end of the 
    // insertion because the last delta-step with insertion should be a 
    // full insertion step. I will therefore floor insert_delta_time
    // and use the remainder for a final partial insertion step. 
    // This introduces a slight error in the timing of the insertion, 
    // which is tolerable.
    double insert_delta_time_dbl = *insert_time * *insert_t_step / *delta;
    double remainder = insert_delta_time_dbl - ftrunc(insert_delta_time_dbl);
    int insert_delta_time = (int)(ftrunc(insert_delta_time_dbl));

    while(delta_time < (*max_times + *delta)){
    	// Insert biomass from the outside
    	if((insert_counter < insert_delta_time) && (delta_time > 0)){
    		for(i = 0; i < num_neighbours[infl_poly_c]; i++){
    			if(nk[start_nk[infl_poly_c] + i] == *num_poly){
    				N[infl_poly_c][i][0] += (*delta / *insert_t_step)*((double)100 / *insert_time);
    			}
    		}
    		insert_counter += 1;
    	}
 
    	if((insert_counter == insert_delta_time) && (delta_time > 0)){
    		// Last, partial insertion step
    		for(i = 0; i < num_neighbours[infl_poly_c]; i++){
    			if(nk[start_nk[infl_poly_c] + i] == *num_poly){
    				N[infl_poly_c][i][0] += (*delta / *insert_t_step)*((double)100 / *insert_time) * remainder;
    			}
    		}
    		insert_counter += 1;
        }   	                    

    	// Calculate the distribution
    	for(i = 0; i < (*num_poly-1); i++){
            dbl_sum = 0;
            for(j = 0; j < num_neighbours[i]; j++){
                for(age = 0; age < *num_age_classes; age++){          
                    dbl_sum += N[i][j][age];
                }
            }
    		N_delta[counter * *num_poly + i] = dbl_sum;
    	}
    	N_delta[counter * *num_poly + *num_poly -1] = outside;


    	// Determine biomass to move
    	for(i=0; i < (*num_poly -1); i++){
    		for(j=0; j < num_neighbours[i]; j++){
                dbl_sum = 0;
                for(k=0; k < num_neighbours[i]; k++){
                    for(age=0; age < *num_age_classes; age++){
                        dbl_sum += N[i][k][age] 
                          * P[start_P[i] + k +  (j * num_neighbours[i]) 
                          + age*(num_neighbours[i] * num_neighbours[i])];
                    }
                }
                M[i][nk[start_nk[i] + j] -1] = dbl_sum;
    		}
    	}	

	
    	// Move biomass between the age classes:
    	for(i=0; i < (*num_poly -1); i++){
            for(j=0; j < num_neighbours[i]; j++){
                for(age=0; age < *num_age_classes; age++){
                    dbl_sum = 0;
                    for(k=0; k < num_neighbours[i]; k++){    
                        dbl_sum += P[start_P[i] + j +  (k * num_neighbours[i])
                            + age*(num_neighbours[i] * num_neighbours[i])];
                    }
                    N[i][j][age] -= N[i][j][age] * dbl_sum;
                }
            }
    		for(age = (*num_age_classes-1); age > 0; age--){
                for(j=0; j < num_neighbours[i]; j++){
                    N[i][j][age] += N[i][j][age-1];
    			    N[i][j][age-1] = 0;
                }
    		}
    	}    
    	
    	// Move biomass between polygons
    	for(i=0; i < (*num_poly -1); i++){
            for(j=0; j < num_neighbours[i]; j++){
                N[i][j][0] = M[nk[start_nk[i] + j] -1][i];
                M[nk[start_nk[i] + j] -1][i] = 0;
            }
        }
        for(i=0; i < (*num_poly -1); i++){
            outside += M[i][*num_poly -1];
            M[i][*num_poly -1] = 0;
        }
        
        delta_time = delta_time + *delta;
    	counter = counter + 1;
    }
}


/**********************************************************************
 *
 * biomass_tracking_start_setup
 *
 * This function works exactly as "biomass_tracking", with the only
 * difference being that the here no inflow into a polygon is used 
 * to introduce biomass, but a starting biomass distribution is passed
 * as a parameter.
 *
 **********************************************************************/

void biomass_tracking_start_setup(double *N_delta, double *max_times,
                       double *P, int *num_neighbours, int *num_poly,
                       int *nk, int *num_age_classes, double *delta,
                       double* start_setup){
                       
    int i,j,k;
    int age, oldest_moving_age;
   
    double dbl_sum = 0;
    int counter = 0;    
    
    // Calculate the starting points of the polygons in P and nk
    int *start_P = (int*) R_alloc(*num_poly,sizeof(int));
    int *start_nk = (int*) R_alloc(*num_poly,sizeof(int));
    start_P[0] = 0;
    start_nk[0] = 0;
    for(i = 1; i < *num_poly; i++){
        start_P[i] = start_P[i-1] 
            + num_neighbours[i-1] * num_neighbours[i-1] * *num_age_classes;
        start_nk[i] = start_nk[i-1] + num_neighbours[i-1];
    }           
    
    // N holds the current biomass distribution. Its first index is the
    // polygon, second index is origin and third index age class.
    // It does not include the outside!
    double ***N = (double***) R_alloc(*num_poly,sizeof(double**));
    for(i = 0; i < (*num_poly -1); i++){
        N[i] = (double**) R_alloc(num_neighbours[i],sizeof(double*));
        for(j = 0; j < num_neighbours[i]; j++){
            N[i][j] = (double*) R_alloc(*num_age_classes,sizeof(double));
            for(age = 0; age < *num_age_classes; age++){
                N[i][j][age] = 0;
            }
        }
    }
    double outside = 0;

    // Distribute the starting biomass as "native" biomass
    for(i=0; i < (*num_poly-1); i++){
        // See what the latest age-class is in which native biomass moves
        oldest_moving_age = *num_age_classes -1;
        for(age=(*num_age_classes-1); age >= 0; age--){
            dbl_sum = 0;
            for(j = 0; j < num_neighbours[i]; j++){
                dbl_sum += N[i][j][age];
            }
            if(dbl_sum > 0){
                oldest_moving_age = age;
                break;
            }
        }
        //for(age=0; age <= oldest_moving_age; age++){
        //    N[i][i][age] = start_setup[i] / (oldest_moving_age+1);
        //}
        N[i][local_neighbour(i,i,nk,start_nk,num_neighbours[i])][0] = start_setup[i];
    }
    outside = start_setup[*num_poly -1];
    
    // M[i][j] holds the biomass to be moved from i to j
    // Put an additional row into M containing the inflow from the outside
   	double **M = (double**) R_alloc(*num_poly,sizeof(double*));
   	for(i=0; i < *num_poly; i++){
        M[i] = (double*) R_alloc(*num_poly,sizeof(double));
        for(j=0; j < *num_poly; j++){
            M[i][j] = 0;
        }
    }
    
     
    double delta_time = 0;
    counter = 0;

    while(delta_time < (*max_times + *delta)){
    	// Calculate the distribution
    	for(i = 0; i < (*num_poly-1); i++){
            dbl_sum = 0;
            for(j = 0; j < num_neighbours[i]; j++){
                for(age = 0; age < *num_age_classes; age++){          
                    dbl_sum += N[i][j][age];
                }
            }
    		N_delta[counter * *num_poly + i] = dbl_sum;
    	}
    	N_delta[counter * *num_poly + *num_poly -1] = outside;


    	// Determine biomass to move
    	for(i=0; i < (*num_poly -1); i++){
    		for(j=0; j < num_neighbours[i]; j++){
                dbl_sum = 0;
                for(k=0; k < num_neighbours[i]; k++){
                    for(age=0; age < *num_age_classes; age++){
                        dbl_sum += N[i][k][age] 
                          * P[start_P[i] + k +  (j * num_neighbours[i]) 
                          + age*(num_neighbours[i] * num_neighbours[i])];
                    }
                }
                M[i][nk[start_nk[i] + j] -1] = dbl_sum;
    		}
    	}	

	
    	// Move biomass between the age classes:
    	for(i=0; i < (*num_poly -1); i++){
            for(j=0; j < num_neighbours[i]; j++){
                for(age=0; age < *num_age_classes; age++){
                    dbl_sum = 0;
                    for(k=0; k < num_neighbours[i]; k++){    
                        dbl_sum += P[start_P[i] + j +  (k * num_neighbours[i])
                            + age*(num_neighbours[i] * num_neighbours[i])];
                    }
                    N[i][j][age] -= N[i][j][age] * dbl_sum;
                }
            }
    		for(age = (*num_age_classes-1); age > 0; age--){
                for(j=0; j < num_neighbours[i]; j++){
                    N[i][j][age] += N[i][j][age-1];
    			    N[i][j][age-1] = 0;
                }
    		}
    	}    
    	
    	// Move biomass between polygons
    	for(i=0; i < (*num_poly -1); i++){
            for(j=0; j < num_neighbours[i]; j++){
                N[i][j][0] = M[nk[start_nk[i] + j] -1][i];
                M[nk[start_nk[i] + j] -1][i] = 0;
            }
        }
        for(i=0; i < (*num_poly -1); i++){
            outside += M[i][*num_poly -1];
            M[i][*num_poly -1] = 0;
        }
        
        delta_time = delta_time + *delta;
    	counter = counter + 1;
    }
}




/* Here come the local helper functions which are accessed by the R scripts */

/* Get position of a single particle */
int get_single_position(double x, double y, int *S, int num_poly, int *m, int *n){
    int a,b;
    
    // Give a particle position -10,-10 if it has left the domain
    if((x >= (*n + 0.5)) || (x <= 0.5) 
    || (y >= (*m + 0.5)) || (y <= 0.5)){
    	x = -10;
    	y = -10;
    	return(num_poly);
    }
    else{
        a = (int) fround(x,0);
        b = (int) fround(y,0);
        return(S[(a-1)*(*m) + b-1]);
    }
}


/* Get the polygon position of all particles */
void get_position(int *pos, double *x, double *y, int *num_part, int *S, int num_poly, int *m, int *n){
    int i;
    
    for(i=0;i < *num_part;i++){
        pos[i] = get_single_position(x[i],y[i],S,num_poly,m,n);
    }
}

/* Calculate new coordinates of all particles.
   re_insert determines whether particles should be re_inserted and
   subdiv gives the number of calculation steps per time step.
   Note that re_insert is not used at the moment, but it is kept in 
   this function since it might be handy in the future. */
void move_particles(double *x,double *y,double *U,double *V,int *m,int *n,int *num_part,
                     int *S,int num_poly,double *t_step,int subdiv, int re_insert, double *infl_prob,
                     double mean_abs, double diff_sd){
                            
       GetRNGstate();
       
       int i,counter;
       int k,b,left_end,right_end;
       double rand_no;
       double little_t_step = (*t_step) / subdiv;

       for(counter=0; counter < subdiv; counter++){
       for(i=0;i < *num_part; i++){   
            // Do not move particle which have left anyway
            if((re_insert == 0) 
              && ((x[i] >= (*n + 0.5)) || (x[i] <= 0.5) || (y[i] >= (*m + 0.5)) || (y[i] <= 0.5))){
              // Put away properly
              x[i] = -10;
              y[i] = -10;
            }
            else{            
                // Advection
                runge_kutta(&(x[i]),&(y[i]),U,V,m,n,&little_t_step);
    		    // Diffusion               
    			x[i] = x[i] + (mean_abs * diff_sd * norm_rand() * little_t_step);
	    		y[i] = y[i] + (mean_abs * diff_sd * norm_rand() * little_t_step);
    			
		        // Has the particle hit land?
                if(get_single_position(x[i],y[i],S,num_poly,m,n) == 0){
                    x[i] = -10;
                    y[i] = -10;
                }
    			
    			// Reinsert if needed and wanted
    			if((re_insert == 1) 
                  && ((x[i] >= (*n + 0.5)) || (x[i] <= 0.5) || (y[i] >= (*m + 0.5)) || (y[i] <= 0.5))){
                    // Insert the particles again which have left the domain
     
                    rand_no = unif_rand();
                    // Search for entry point using divide-and-conquer strategy
                    b = 1;
                    left_end = 0;
                    right_end = (2*(*m) + 2*(*n))-1;
                    while(b==1){
                        k = (right_end + left_end)/2;
                        if ((k > 0) && (infl_prob[k-1] >= rand_no)){
                            right_end = k;
                        }
                        else if((infl_prob[k] <= rand_no) && (k < ((2*(*m) * 2*(*n))-1))){
                            left_end = k+1;
                        }
                        else b = 0;
                        if(k == ((2*(*m) * 2*(*n))-1)) b = 0;
                        if(k==0) b = 0;
                        if(left_end == right_end) b = 0;
                    }
                    
                    // Get particle position at entry point
                    if(k < (*m)){
                         x[i] = 0.5001;
                         y[i] = k + 0.999*unif_rand() - 0.4999;
                    }
                    else if(k < 2*(*m)){
                         x[i] = *n + 0.4999;
                         y[i] = k - (*m) +1 + 0.999*unif_rand() - 0.4999;
                    }
                    else if(k < (2*(*m) + (*n))){
                         x[i] = k - 2*(*m) +1 + 0.999*unif_rand() - 0.4999;
                         y[i] = 0.5001;
                    }
                    else{
                         x[i] = k - 2*(*m) - (*n) +1 + 0.999*unif_rand() - 0.4999;
                         y[i] = *m + 0.4999;
                    }
              }
          }             
      }
      }
      
      PutRNGstate();
}



// Return the local neighbourhood number of a polygon relative to another polygon
// Poly and neighbour have to be given in the C notation, i.e. 0...num_poly-1
int local_neighbour(int poly, int neighbour, int* nk, int* start_nk, int num_neighbours){
    int i;
    neighbour += 1;
    for(i=0; i < num_neighbours; i++){
        if(nk[start_nk[poly] + i] == neighbour) return(i);
    }     
    error("The local neighbour could not be found!\n");
    return(0);
}


// Advance a single particle by one step using a 4th order runge kutta algorithm
void runge_kutta(double *x, double *y, double *U, double *V, int *m, int *n, double *t_step){
    int lower_x_coord,lower_y_coord;
	double dx, dy;
	double u_velo1, u_velo2, u_velo3, u_velo4;
	double v_velo1, v_velo2, v_velo3, v_velo4;
	double u_velo, v_velo;
	double tmp_x, tmp_y;
	int k;
     
    // First coefficient
    // Get coordinates      
    lower_x_coord = (int) ftrunc(*x);
    lower_y_coord = (int) ftrunc(*y);
		
	// Get distances and therefore weights
	dx = *x - lower_x_coord;
	dy = *y - lower_y_coord;
		
	// Take care of exceptions!
	lower_x_coord = imin2(imax2(lower_x_coord,1),(*n -1));
	lower_y_coord = imin2(imax2(lower_y_coord,1),(*m -1));
	
	// Get velocities
	k = (lower_x_coord -1)*(*m) + lower_y_coord -1; 
	// -1 because array in numbered 0..n-1
	u_velo1 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(*m)] 
	           + dy*(1-dx)*U[k+1] + dx*dy*U[k+(*m)+1];
	v_velo1 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(*m)] 
	           + dy*(1-dx)*V[k+1] + dx*dy*V[k+(*m)+1];

	// Second coefficient
    // Get coordinates      
    tmp_x = *x + (*t_step * u_velo1 / 2);
    tmp_y = *y + (*t_step * v_velo1 / 2);
	lower_x_coord = (int) ftrunc(tmp_x);
	lower_y_coord = (int) ftrunc(tmp_y);
		
	// Get distances and therefore weights
	dx = tmp_x - lower_x_coord;
	dy = tmp_y - lower_y_coord;
			
	// Take care of exceptions!
	lower_x_coord = imin2(imax2(lower_x_coord,1),(*n -1));
	lower_y_coord = imin2(imax2(lower_y_coord,1),(*m -1));
			
	// Get velocities
	k = (lower_x_coord -1)*(*m) + lower_y_coord -1; 
	// -1 because array in numbered 0..n-1
	u_velo2 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(*m)] 
	           + dy*(1-dx)*U[k+1] + dx*dy*U[k+(*m)+1];
	v_velo2 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(*m)] 
	           + dy*(1-dx)*V[k+1] + dx*dy*V[k+(*m)+1];
	          
	// Third coefficient
    // Get coordinates   
    tmp_x = *x + (*t_step * u_velo2 / 2);
    tmp_y = *y + (*t_step * v_velo2 / 2);   
	lower_x_coord = (int) ftrunc(tmp_x);
	lower_y_coord = (int) ftrunc(tmp_y);
	
	// Get distances and therefore weights
	dx = tmp_x - lower_x_coord;
	dy = tmp_y - lower_y_coord;
			
	// Take care of exceptions!
	lower_x_coord = imin2(imax2(lower_x_coord,1),(*n -1));
	lower_y_coord = imin2(imax2(lower_y_coord,1),(*m -1));
		
	// Get velocities
	k = (lower_x_coord -1)*(*m) + lower_y_coord -1; 
	// -1 because array in numbered 0..n-1
	u_velo3 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(*m)] 
	           + dy*(1-dx)*U[k+1] + dx*dy*U[k+(*m)+1];
    v_velo3 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(*m)] 
	           + dy*(1-dx)*V[k+1] + dx*dy*V[k+(*m)+1];
	
    // Fourth coefficient
    // Get coordinates  
    tmp_x = *x + (*t_step * u_velo3);
    tmp_y = *y + (*t_step * v_velo3);
	lower_x_coord = (int) ftrunc(tmp_x);
	lower_y_coord = (int) ftrunc(tmp_y);
		
	// Get distances and therefore weights
	dx = tmp_x - lower_x_coord;
	dy = tmp_y - lower_y_coord;
			
	// Take care of exceptions!
	lower_x_coord = imin2(imax2(lower_x_coord,1),(*n -1));
	lower_y_coord = imin2(imax2(lower_y_coord,1),(*m -1));
			
	// Get velocities
	k = (lower_x_coord -1)*(*m) + lower_y_coord -1; 
	// -1 because array in numbered 0..n-1
	u_velo4 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(*m)] 
	           + dy*(1-dx)*U[k+1] + dx*dy*U[k+(*m)+1];
	v_velo4 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(*m)] 
	           + dy*(1-dx)*V[k+1] + dx*dy*V[k+(*m)+1];
		
	// Combine coefficients
	u_velo = u_velo1/6 + u_velo2/3
	          + u_velo3/3 + u_velo4/6;
	v_velo = v_velo1/6 + v_velo2/3
	          + v_velo3/3 + v_velo4/6;
                            
	// Move particles, add some random diffusion
    *x = *x + (u_velo * *t_step);
	*y = *y + (v_velo * *t_step);
}




// Get inflow probabilities from outside the boundary
// First the probabilities for inflow from the left, than right,
// than below and then above
double* get_infl_prob(double *U, double *V, int *S, int *m, int *n){
    int i;
    double dbl_tmp;

	double* infl_prob = (double*) R_alloc(2*(*m) + 2*(*n), sizeof(double));
	for(i=0; i < (2*(*m) + 2*(*n)); i++){
        infl_prob[i] = 0;
    }
	for (i=0;i < (2*(*m) + 2*(*n));i++){
        infl_prob[i] = 0;
    }
	dbl_tmp=0;
	for (i=0;i < *m;i++){
        if((U[i] > 0) && (S[i] != 0)) infl_prob[i] = U[i];
        if((U[(*n-1)*(*m) + i] < 0) && (S[(*n-1)*(*m) +i] != 0)) infl_prob[i+(*m)] = -U[(*n-1)*(*m) + i];
        dbl_tmp += (infl_prob[i] + infl_prob[i+(*m)]);
    }
    for (i=0;i < *n;i++){
        if((V[i*(*m)] > 0) && (S[i*(*m)] != 0)) infl_prob[2*(*m) + i] = V[i*(*m)];
        if((V[i*(*m) + (*m) -1] < 0) && (S[i*(*m) + (*m) -1] != 0)) infl_prob[2*(*m) + (*n) +i] = -V[i*(*m) + (*m) -1];
        dbl_tmp += (infl_prob[2*(*m) +i] + infl_prob[2*(*m)+(*n)+i]);
    }
    // Check whether there is any inflow at all
    if(dbl_tmp == 0){
        for (i=0;i < (2*(*m) + 2*(*n));i++){
            infl_prob[i] = 1 / (2*(*m) + 2*(*n));
            warning("There is no inflow into the domain. Is that sensible?\n");
        }
    }
    for (i=0;i < *m;i++){
        infl_prob[i] /= dbl_tmp;
        infl_prob[i+(*m)] /= dbl_tmp;
    }
    for (i=0;i < *n;i++){
        infl_prob[2*(*m) + i] /= dbl_tmp;
        infl_prob[2*(*m) + (*n) +i] /= dbl_tmp;
    }
    
    // Build cumulative distribution function
    dbl_tmp = 0;
    for (i=0;i < (2*(*m) + 2*(*n)); i++){
        dbl_tmp += infl_prob[i];
        infl_prob[i] = dbl_tmp;             
    }   
    
    return(infl_prob);
}
