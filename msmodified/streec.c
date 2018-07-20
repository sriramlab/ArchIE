/**********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure. 
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ms.h"
#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

#define SEGINC 80 

extern int flag;

int nchrom, begs, nsegs;
long nlinks ;
// Number of nodes and hence, the pointer to the root, in each segment.
// Initially it is nsam - 1
static int *nnodes = NULL ;  
double t, cleft , pc, lnpc ;

static unsigned seglimit = SEGINC ;
static unsigned maxchr ;

struct seg{
	int beg;
	int end;
	int desc;
	};

//ramcode
//Each chromosome is represented as a struct
// pop = which population
// nseg - number of segments
// pseg - pointer to segments (beg, end, individual number that inherits this chr or more generally the label of the root of the tree for this segment) 
// Important Invariant: For each chr and each segl in seglst, it is either absent in pseg or is contained within a single segment of pseg.
//ramcode

struct chromo{
	int nseg;
	int pop;
	struct seg  *pseg;
	} ;

static struct chromo *chrom = NULL ;

struct node{
	int abv;
	int ndes;
	float time;
	} *ptree1, *ptree2;

struct segl {
	int beg;
	struct node *ptree;
	int next;
    // ramcode
	int *desc;
    // ramcode
	}  ;
// Stores all the segments (crossovers)
// nseg -- number of segments
// Number of segments increases only when there is a xover
static struct segl *seglst = NULL ;

//ramcode
void printptree (struct node *ptree, int nsam){
	int i ;
	for (i  = 0 ; i <2*nsam-1;i++){
		printf ("(%d:%d:%2.3f),",i,(ptree+i)->abv, (ptree+i)->time);
	}
}

// Prints all segments
void printseglist () {
	int i = 0;
	int s  = 0;
	for (; i < nsegs ;  i++, s = seglst[s].next){
		printf ("%d,",seglst[s].beg);
	}

	printf ("\n");
}

void printseg ( int c) { 
	int i =  0;
	struct seg *s = chrom[c].pseg;
	for (; i < chrom[c].nseg; i++){ 
		printf ("(%d,%d,%d),", s[i].beg, s[i].end, s[i].desc);
	}
	printf ("\n");
}


double nhtime = 0;
//ramcode

	struct segl *
segtre_mig(struct c_params *cp, int *pnsegs ) 
{
	int i, j, k, seg, dec, pop, pop2, c1, c2, ind, rchrom, intn  ;
	int migrant, source_pop, *config, flagint ;
	double  ran1(), sum, x, tcoal, ttemp, rft, clefta,  tmin, p  ;
	double prec, cin,  prect, nnm1, nnm0, mig, ran, coal_prob, prob, rdum , arg ;
	char c, event ;
	int re(), xover(), cinr(), cleftr(), eflag, cpop, ic  ;
	int nsam, npop, nsites, nintn, *inconfig ;
	double r,  f, rf,  track_len, *nrec, *npast, *tpast, **migm ;
	double *size, *alphag, *tlast ;
	struct devent *nextevent ;

	nsam = cp->nsam;
    // ramcode
	int *done = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
    // ramcode

	npop = cp->npop;
	nsites = cp->nsites;
	inconfig = cp->config;
	r = cp->r ;
	f = cp->f ;
	track_len = cp->track_len ;
	migm = (double **)malloc( (unsigned)npop*sizeof(double *) ) ;
	for( i=0; i<npop; i++) {
	  migm[i] = (double *)malloc( (unsigned)npop*sizeof( double) ) ;
	  for( j=0; j<npop; j++) migm[i][j] = (cp->mig_mat)[i][j] ;
	  }
	// sc: List of events as given in the commandline
	nextevent = cp->deventlist ;
	
/* Initialization */
	if( chrom == NULL ) {
	   maxchr = nsam + 20 ;
	   chrom = (struct chromo *)malloc( (unsigned)( maxchr*sizeof( struct chromo) )) ;
	  if( chrom == NULL ) perror( "malloc error. segtre");
	  }
	if( nnodes == NULL ){
		nnodes = (int*) malloc((unsigned)(seglimit*sizeof(int)))  ;
		if( nnodes == NULL ) perror("malloc error. segtre_mig");
		}
	if( seglst == NULL ) {
		seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)) ) ;
		if( seglst == NULL ) perror("malloc error. segtre_mig.c 2");
		}

	config = (int *)malloc( (unsigned) ((npop+1)*sizeof(int) )) ;
	if( config == NULL ) perror("malloc error. segtre.");
	size = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( size == NULL ) perror("malloc error. segtre.");
	alphag = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( alphag == NULL ) perror("malloc error. segtre.");
	tlast = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( alphag == NULL ) perror("malloc error. segtre.");
	for(pop=0;pop<npop;pop++) {
	   config[pop] = inconfig[pop] ;
	   size[pop] = (cp->size)[pop] ;
	   alphag[pop] = (cp->alphag)[pop] ;
	   tlast[pop] = 0.0 ;
	}
	for(pop=ind=0;pop<npop;pop++)
		for(j=0; j<inconfig[pop];j++,ind++) {
			
			chrom[ind].nseg = 1;
			if( !(chrom[ind].pseg = (struct seg*)malloc((unsigned)sizeof(struct seg)) ))
			  ERROR("calloc error. se1");

			(chrom[ind].pseg)->beg = 0;
			(chrom[ind].pseg)->end = nsites-1;
			(chrom[ind].pseg)->desc = ind ;
			chrom[ind].pop = pop ;
			}


	int debug = 0;
	// ramcode
	//printf ("nind =  %d, nchrom = %d", ind, nchrom);
	/*
	for (ic =  0 ; ic <  ind; ic++){
		 printf ("Individuals = %d, %d, %d\n",
				(chrom[ic].pseg)->beg,
				(chrom[ic].pseg)->end,
				(chrom[ic].pseg)->desc);

	}*/
	// ramcode

	seglst[0].beg = 0;
	if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
		 perror("calloc error. se2");



	nnodes[0] = nsam - 1 ;
	nchrom=nsam;
	nlinks = ((long)(nsam))*(nsites-1) ;
	nsegs=1;
	t = 0.;
	r /= (nsites-1);
	if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
	else pc = 1.0 ;
	lnpc = log( pc ) ;
	cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
	if( r > 0.0 ) rf = r*f ;
	else rf = f /(nsites-1) ;
	rft = rf*track_len ;
	flagint = 0 ;

    /* Main loop */
    // ramcode
	// eflag = 1 if there is an event : coalescence, migration, recombination etc.
	// tmin is the time to this event
	// config [] - lineages in each population
	// size[] - Ne
	// alphag[] - growth 
    // ramcode

	while( nchrom > 1 ) {
		prec = nlinks*r;
		cin = nlinks*rf ;
		clefta = cleft*rft ;
		prect = prec + cin + clefta ;
		mig = 0.0;
		for( i=0; i<npop; i++) mig += config[i]*migm[i][i] ;
		if( (npop > 1) && ( mig == 0.0) && ( nextevent == NULL)) {
		   i = 0;
		   for( j=0; j<npop; j++) 
			if( config[j] > 0 ) i++;
		   if( i > 1 ) {
			fprintf(stderr," Infinite coalescent time. No migration.\n");
			exit(1);
		   }
		}
		eflag = 0 ;

		if( prect > 0.0 ) {      /* cross-over or gene conversion */
		  while( (rdum = ran1() )  == 0.0 ) ;
		  ttemp = -log( rdum)/prect ;
		  if( (eflag == 0) || (ttemp < tmin ) ){
		    tmin = ttemp;
		    event = 'r' ;
		    eflag = 1;
		  }
	        }
		if(mig > 0.0 ) {         /* migration   */
		  while( (rdum = ran1() ) == 0.0 ) ;
		  ttemp = -log( rdum)/mig ;
		  if( (eflag == 0) || (ttemp < tmin ) ){
		    tmin = ttemp;
		    event = 'm' ;
		    eflag = 1 ;
		  }
	        }
		
	    for(pop=0; pop<npop ; pop++) {     /* coalescent */
		coal_prob = ((double)config[pop])*(config[pop]-1.) ;
	        if( coal_prob > 0.0 ) {
		   while( ( rdum = ran1() )  == .0 ) ;
	  	   if( alphag[pop] == 0 ){
			 ttemp = -log( rdum )*size[pop] /coal_prob ;
		         if( (eflag == 0) || (ttemp < tmin ) ){
		            tmin = ttemp;
		            event = 'c' ;
		            eflag = 1 ;
			    cpop = pop;
			 }
		    }
	  	   else {
		     arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob     ;
		     if( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */ 
		         ttemp = log( arg ) / alphag[pop]  ;
		         if( (eflag == 0) || (ttemp < tmin ) ){
		            tmin = ttemp;
		            event = 'c' ;
		            eflag = 1 ;
			    cpop = pop ;
			 }
		     }
		   }
	         }		
 	      }

	    if( (eflag == 0) && ( nextevent == NULL) ) {
	      fprintf(stderr,
               " infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
	      exit( 0);
	    }
	if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) ) {
		// ramcode
		/*
		printf ("Population event\n");
		for (ic  = 0 ;ic < nchrom; ic++){
				  printf ("Individuals  = %d, %d, %d\n",
				(chrom[ic].pseg)->beg,
				(chrom[ic].pseg)->end,
				(chrom[ic].pseg)->desc);
		}*/
		// ramcode

	    t = nextevent->time ;
	    switch(  nextevent->detype ) {
		case 'N' :
		   for(pop =0; pop <npop; pop++){
			 size[pop]= nextevent->paramv ;
			 alphag[pop] = 0.0 ;
			}
		   nextevent = nextevent->nextde ;
		   break;
		case 'n' :
		   size[nextevent->popi]= nextevent->paramv ;
		   alphag[nextevent->popi] = 0.0 ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'G' :
		   for(pop =0; pop <npop; pop++){
		     size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
		     alphag[pop]= nextevent->paramv ;
		     tlast[pop] = t ;
		   }
		   nextevent = nextevent->nextde ;
		   break;
		case 'g' :
		     pop = nextevent->popi ;
		     size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
		     alphag[pop]= nextevent->paramv ;
		     tlast[pop] = t ;
		     nextevent = nextevent->nextde ;
		     break;
		case 'M' :
		   for(pop =0; pop <npop; pop++)
		     for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
		   for( pop = 0; pop <npop; pop++)
		     migm[pop][pop]= nextevent->paramv ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'a' :
		   for(pop =0; pop <npop; pop++)
		     for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'm' :
		  i = nextevent->popi ;
		  j = nextevent->popj ;
		  migm[i][i] += nextevent->paramv - migm[i][j];
		   migm[i][j]= nextevent->paramv ;
		   nextevent = nextevent->nextde ;
		   break;
	    case 'j' :         /* merge pop i into pop j  (join) */
		  i = nextevent->popi ;
		  j = nextevent->popj ;
          //
          //ramcode
          if (i==2 && j==0) {
              nhtime = nextevent->time;
          }
		  if (debug) {
              printf (" Merging population %d (%d) to population %d (%d) at time %2.3f\n",i,config[i], j, config[j],nhtime);
          }
          //ramcode
          
		  config[j] += config[i] ;
		  config[i] = 0 ;
		  for( ic = 0; ic<nchrom; ic++) if( chrom[ic].pop == i ) chrom[ic].pop = j ;
		/*  the following was added 19 May 2007 */
		  for( k=0; k < npop; k++){
		     if( k != i) {
		        migm[k][k] -= migm[k][i] ;
		        migm[k][i] = 0. ;
			 }
		   }
		/* end addition */
		   nextevent = nextevent->nextde ;
		   break;
	    case 's' :         /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
		  i = nextevent->popi ;
		  p = nextevent->paramv ;
		  npop++;
		  config = (int *)realloc( config, (unsigned)(npop*sizeof( int) )); 
		  size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ));
		  alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ));
		  tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ;
		  tlast[npop-1] = t ;
		  size[npop-1] = 1.0 ;
		  alphag[npop-1] = 0.0 ;
		  migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
		  for( j=0; j< npop-1; j++)
			 migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
		  migm[npop-1] = (double *)malloc( (unsigned)(npop*sizeof( double) ) ) ;
		  for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
		  config[npop-1] = 0 ;
		  config[i] = 0 ;
          // Assign chromosomes to pop i or npop-1 with probability (p,1-p).
		  for( ic = 0; ic<nchrom; ic++){
		    if( chrom[ic].pop == i ) {
		      if( ran1() < p ) config[i]++;
		      else {
    			 chrom[ic].pop = npop-1 ;
	    		 config[npop-1]++;
		      }
		    }
		  }
		   nextevent = nextevent->nextde ;
		  //ramcode
		  if (debug) 
			 printf ("Split\n");

     	  for (ic = 0 ; ic < nsegs; ic++) {
			 seglst[ic].desc =  (int *)malloc( (unsigned)(nsam)*sizeof( int) );
			 int j  ;
  			 for (j = 0 ; j < nsam; j++) 
				 seglst[ic].desc[j] = 0;
			 
		  }
		 
		  if (debug)
			  printf ("nchrom = %d, npop = %d, nsam = %d\n", nchrom, npop, nsam);

		  for (ic  = 0 ;ic < nchrom; ic++){
			  if (chrom[ic].pop == npop-2){
				  if (debug)
					  printf ("ic=%d\t",ic);

				  struct seg *s = chrom[ic].pseg;
				  for (i = 0 ; i < 2*nsam - 1; i++)
					done[i] = 0;

				  if (debug) {
					  int j  = 0 ;
					  for (  j = 0 ;  j < chrom[ic].nseg; j++) {
						  printf ("(%d,%d),",s[j].beg, s[j].end);
					  }
					  printf ("\n");
				  }

                  for (  j = 0 ;  j < chrom[ic].nseg; j++) {
                      int l = 0;
                      int k = 0;
                      for( l=0,k=0; (k<nsegs)&&(s[j].beg > seglst[l].beg);
                              l=seglst[l].next, k++) ;
                      getdesc (seglst[l].ptree, nsam, s[j].desc, done);

                      if (debug){
                          for (i  = 0; i < 2*nsam-1; i++)
                              printf ("%d,",done[i]);
                          printf ("\n");
                      }
                  }
			  }
		  }

		  for (ic  = 0 ;ic < nchrom; ic++){
              // for all chromosomes that are ascending the Neandertal branch (picked by their population number)
              // for each segment, find all individuals that descend from that segment
              //
			  if (chrom[ic].pop== npop-1 || chrom[ic].pop == 2){
				  
				  	struct seg *s = chrom[ic].pseg;
					int i  = 0 ;
					int j  = 0;
					for (i = 0 ; i < 2*nsam - 1; i++)
						done[i] = 0;
					
					int k = 0 ;
					for (i = 0; i < chrom[ic].nseg; i++) {

                        // Discard all segments that are not part of this chromosome
						for( j=0,k=0; (k<nsegs)&&(s[i].beg > seglst[j].beg);
					   		j=seglst[j].next, k++) ;

                        // for all segments that are in this chromosome
                        // find all descendants (done [l] == 2)
                        // and store in seglst[j].desc
						for (;(k<nsegs) && (s[i].end >= seglst[j].beg);
							j=seglst[j].next, k++){
							getdesc (seglst[j].ptree, nsam, s[i].desc, done);				


							int l = 0 ;
							for (l = 0 ; l < nsam; l++)
								if (done[l]==2)
									seglst[j].desc[l]  = 1;
						
						}

					}
			  }
		  }

		  if (debug){
			  printf ("More desc\n");
			  int k =  0;
			  for (ic = 0 ; ic < nsegs; k=seglst[k].next, ic++){
				  int *desc  = seglst[k].desc;
				  int j ;
				  printf ("%d\t",seglst[k].beg);
				  for ( j = 0  ; j < nsam; j++){
					  printf ("%d,",desc[j]);
				  }
				  printf ("\n");
			  }
		  }
          //ramcode
          //
		  break;
		}
 	   } 
    else {
        t += tmin ;	
        // ramcode
        if (debug) {
            printf ("event = %c\n", event);
        }
        // ramcode
        if( event == 'r' ) {   
            if( (ran = ran1()) < ( prec / prect ) ){ /*recombination*/
                rchrom = re(nsam);
                config[ chrom[rchrom].pop ] += 1 ;
            }
            else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
                rchrom = cleftr(nsam);
                config[ chrom[rchrom].pop ] += 1 ;
            }
            else  {         /* cin event */
                rchrom = cinr(nsam,nsites);
                if( rchrom >= 0 ) config[ chrom[rchrom].pop ] += 1 ;
            }
            //ramcode
            if (debug){
                printf ("Recombination event:");
                printf ("Recombined %d to give (%d,%d,%d),(%d,%d,%d)\n", rchrom, (chrom[rchrom].pseg)->beg,(chrom[rchrom].pseg)->end,(chrom[rchrom].pseg)->desc,(chrom[nchrom-1].pseg)->beg,(chrom[nchrom-1].pseg)->end,(chrom[nchrom-1].pseg)->desc);
                printf ("New chromosomes\n");
                printseg (rchrom);
                printseg (nchrom-1);
                printseglist (nsegs);
            }
            //ramcode

        }
        else if ( event == 'm' ) {  /* migration event */
            x = mig*ran1();
            sum = 0.0 ;
            for( i=0; i<nchrom; i++) {
                sum += migm[chrom[i].pop][chrom[i].pop] ;
                if( x <sum ) break;
            }
            migrant = i ;
            x = ran1()*migm[chrom[i].pop][chrom[i].pop];
            sum = 0.0;
            for(i=0; i<npop; i++){
                if( i != chrom[migrant].pop ){
                    sum += migm[chrom[migrant].pop][i];
                    if( x < sum ) break;
                }
            }
            source_pop = i;
            config[chrom[migrant].pop] -= 1;
            config[source_pop] += 1;
            chrom[migrant].pop = source_pop ;
        }
        else { 								 /* coalescent event */
            /* pick the two, c1, c2  */
            pick2_chrom( cpop, config, &c1,&c2);  /* c1 and c2 are chrom's to coalesce */

            //ramcode
            if (debug){
                printf ("Coalescence event:");
                printf ("Coalescing %d(%d) and %d(%d)\n", c1, chrom[c1].nseg, c2, chrom[c2].nseg);
                printseg (c1);
                printseg (c2);
                printseglist ();
                printf ("Coalescing at time %f\n", t);
            }
            //ramcode

            dec = ca(nsam,nsites,c1,c2 );
            config[cpop] -= dec ;
            //ramcode
            // At the end of this process
            // c1 contains the new chromosome, c2 contains the last chromosome , nchrom = nchrom - 1
            if (debug){
                printf ("Post-coalescence event:");
                printseg (c1);
                printseg (c2);
                printseglist ();
            }
            /*
               printf ("forming %d(%d) and %d(%d)\n", c1, chrom[c1].nseg, c2, chrom[c2].nseg);
               printseg (c1);
               printseg (c2);*/

            /*
               int i = 0;
               struct seg *s = chrom[c1].pseg;
               printf ("Trees\n");
               for (; i < chrom[c1].nseg; i++) {
               printf ("Beg = %d, Desc = %d\t",s[i].beg,s[i].desc);
               prforest (seglst[i].ptree, nsam, s[i].desc);
               printptree (seglst[i].ptree, nsam);
               int j  ;
               for (j = 0 ; j < 2*nsam - 1; j++)
               done[j] = 0;
               getdesc (seglst[i].ptree, nsam, s[i].desc, done);				
               printf ("\nDescendants\t");
               for (j = 0 ; j < nsam; j++){
               if (done[j]==2)
               printf ("%d,",j);
               }
               printf ("\n");
               }*/
            //			printf ("Forming %d(%d,%d,%d),%d(%d,%d,%d)\n", c1,(chrom[c1].pseg)->beg,(chrom[c1].pseg)->end,(chrom[c1].pseg)->desc,c2,(chrom[c2].pseg)->beg,(chrom[c2].pseg)->end,(chrom[c2].pseg)->desc);
            //ramcode
        }
    }
	}
    //ramcode
    if (debug){
        printf ("desc\n");
        for (i = 0 ; i < nsegs; i++){
            int *desc = seglst[i].desc;
            int j =  0;
            for (j = 0 ; j < nsam;j++){
                printf ("%d,",desc[j]);
            }
            printf ("\n");
        }
	}
    //ramcode

	*pnsegs = nsegs ;
	free(config); 
	free( size ) ;
	free( alphag );
	free( tlast );
	for( i=0; i<npop; i++) free ( migm[i] ) ;
	free( migm ) ;

    //ramcode
	free (done);
    //ramcode
	return( seglst );
}

/******  recombination subroutine ***************************
   Picks a chromosome and splits it in two parts. If the x-over point
   is in a new spot, a new segment is added to seglst and a tree set up
   for it.   ****/


	int
re(nsam)
	int nsam;
{
	struct seg *pseg ;
	int  el, lsg, lsgm1,  ic,  is, in, spot;
	double ran1();


/* First generate a random x-over spot, then locate it as to chrom and seg. */

	spot = nlinks*ran1() + 1.;

    /* get chromosome # (ic)  */

	for( ic=0; ic<nchrom ; ic++) {
		lsg = chrom[ic].nseg ;
		lsgm1 = lsg - 1;
		pseg = chrom[ic].pseg;
		el = ( (pseg+lsgm1)->end ) - (pseg->beg);
		if( spot <= el ) break;
		spot -= el ;
		}
	is = pseg->beg + spot -1;
	xover(nsam, ic, is);
	return(ic);	
}

	int
cleftr( int nsam)
{
	struct seg *pseg ;
	int   lsg, lsgm1,  ic,  is, in, spot;
	double ran1(), x, sum, len  ;
	int xover(int, int, int);

while( (x = cleft*ran1() )== 0.0 )  ;
	sum = 0.0 ;
	ic = -1 ;
	while ( sum < x ) {
		sum +=  1.0 - pow( pc, links(++ic) )  ;
		}
	pseg = chrom[ic].pseg;
	len = links(ic) ;
	is = pseg->beg + floor( 1.0 + log( 1.0 - (1.0- pow( pc, len))*ran1() )/lnpc  ) -1  ;
	xover( nsam, ic, is);
	return( ic) ;
}

int
cinr( int nsam, int nsites)
{
	struct seg *pseg ;
	int len,  el, lsg, lsgm1,  ic,  is, in, spot, endic ;
	double ran1();
	int xover(), ca() ;


/* First generate a random x-over spot, then locate it as to chrom and seg. */

	spot = nlinks*ran1() + 1.;

    /* get chromosome # (ic)  */

	for( ic=0; ic<nchrom ; ic++) {
		lsg = chrom[ic].nseg ;
		lsgm1 = lsg - 1;
		pseg = chrom[ic].pseg;
		el = ( (pseg+lsgm1)->end ) - (pseg->beg);
		if( spot <= el ) break;
		spot -= el ;
		}
	is = pseg->beg + spot -1;
	endic = (pseg+lsgm1)->end ;
	xover(nsam, ic, is);

	len = floor( 1.0 + log( ran1() )/lnpc ) ;
	if( is+len >= endic ) return(ic) ;  
	if( is+len < (chrom[nchrom-1].pseg)->beg ){
	   ca( nsam, nsites, ic, nchrom-1);
	    return(-1) ;
	    }
	xover( nsam, nchrom-1, is+len ) ;
	ca( nsam,nsites, ic,  nchrom-1);
	return(ic);	

}

// Make changes to the data for a xover event.
// ic - chromosome, is - xover location
	int
xover(int nsam,int ic, int is)
{
	struct seg *pseg, *pseg2;
	int i,  lsg, lsgm1, newsg,  jseg, k,  in, spot;
	double ran1(), len ;


	pseg = chrom[ic].pseg ;
	lsg = chrom[ic].nseg ;
	len = (pseg + lsg -1)->end - pseg->beg ;
	cleft -= 1 - pow(pc,len) ;
   /* get seg # (jseg)  */

    // jseg -- number of segments retained on current chr (left part)
    // newsg -- number of segments moved to new chr (right part) 
	for( jseg=0; is >= (pseg+jseg)->end ; jseg++) ;
	if( is >= (pseg+jseg)->beg ) in=1;
	else in=0;
	newsg = lsg - jseg ;

   /* copy last part of chrom to nchrom  */

	nchrom++;
	if( nchrom >= maxchr ) {
	    maxchr += 20 ;
	    chrom = (struct chromo *)realloc( chrom, (unsigned)(maxchr*sizeof(struct chromo))) ;
	    if( chrom == NULL ) perror( "malloc error. segtre2");
	    }
	if( !( pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg)) ) )
		ERROR(" alloc error. re1");
	chrom[nchrom-1].nseg = newsg;
	chrom[nchrom-1].pop = chrom[ic].pop ;
    // Pointer to segment for new chr
	pseg2->end = (pseg+jseg)->end ;
	if( in ) {
		pseg2->beg = is + 1 ;
		(pseg+jseg)->end = is;
		}
	else pseg2->beg = (pseg+jseg)->beg ;
	pseg2->desc = (pseg+jseg)->desc ;
	for( k=1; k < newsg; k++ ) {
		(pseg2+k)->beg = (pseg+jseg+k)->beg;
		(pseg2+k)->end = (pseg+jseg+k)->end;
		(pseg2+k)->desc = (pseg+jseg+k)->desc;
		}

	lsg = chrom[ic].nseg = lsg-newsg + in ;
	lsgm1 = lsg - 1 ;
	nlinks -= pseg2->beg - (pseg+lsgm1)->end ;
	len = (pseg+lsgm1)->end - (pseg->beg) ;
	cleft += 1.0 - pow( pc, len) ;
	len = (pseg2 + newsg-1)->end - pseg2->beg ;
	cleft += 1.0 - pow(pc, len) ;
if( !(chrom[ic].pseg = 
     (struct seg *)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)) )) )
		perror( " realloc error. re2");
	if( in ) {
		begs = pseg2->beg;
		for( i=0,k=0; (k<nsegs-1)&&(begs > seglst[seglst[i].next].beg-1);
		   i=seglst[i].next, k++) ;
		if( begs != seglst[i].beg ) {
						/* new tree  */

	   	   if( nsegs >= seglimit ) {  
	   	   	  seglimit += SEGINC ;
	   	      nnodes = (int *)realloc( nnodes,(unsigned)(sizeof(int)*seglimit)) ; 
	   	      if( nnodes == NULL) perror("realloc error. 1. segtre_mig.c");
	   	      seglst =
	   	      	 (struct segl *)realloc( seglst,(unsigned)(sizeof(struct segl)*seglimit));
	   	      if(seglst == NULL ) perror("realloc error. 2. segtre_mig.c");
	   	      /*  printf("seglimit: %d\n",seglimit);  */
	   	      } 
	   	   seglst[nsegs].next = seglst[i].next;
	   	   seglst[i].next = nsegs;
	   	   seglst[nsegs].beg = begs ;

           //ramcode
		   seglst[nsegs].desc = seglst[i].desc;
           //ramcode
		   if( !(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam), sizeof(struct
			 node)) )) perror("calloc error. re3.");
		   nnodes[nsegs] = nnodes[i];
		   ptree1 = seglst[i].ptree;
		   ptree2 = seglst[nsegs].ptree;
		   nsegs++ ;
		   for( k=0; k<=nnodes[i]; k++) {
		      (ptree2+k)->abv = (ptree1+k)->abv ;
		      (ptree2+k)->time = (ptree1+k)->time;
		      }
		   }
	}
	return(ic) ;
}

/***** common ancestor subroutine **********************
   Pick two chromosomes and merge them. Update trees if necessary. **/

	int
ca(nsam, nsites,c1,c2)
	int nsam,c1,c2;
	int  nsites;
{
	int yes1, yes2, seg1, seg2, seg ;
	int tseg, start, end, desc, k;
	struct seg *pseg;
	struct node *ptree;

	seg1=0;
	seg2=0;

	if( !(pseg = (struct seg *)calloc((unsigned)nsegs,sizeof(struct seg) ))) 
		perror("alloc error.ca1");

	tseg = -1 ;

	for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		start = seglst[seg].beg;
		yes1 = isseg(start, c1, &seg1);
		yes2 = isseg(start, c2, &seg2);
		if( yes1 || yes2 ) {
			tseg++;
			(pseg+tseg)->beg=seglst[seg].beg;
			end = ( k< nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1 ) ;
			(pseg+tseg)->end = end ;

			if( yes1 && yes2 ) {
				nnodes[seg]++;
				if( nnodes[seg] >= (2*nsam-2) ) tseg--;
				else
					(pseg+tseg)->desc = nnodes[seg];
				ptree=seglst[seg].ptree;
				desc = (chrom[c1].pseg + seg1) ->desc;
				(ptree+desc)->abv = nnodes[seg];
				desc = (chrom[c2].pseg + seg2) -> desc;
				(ptree+desc)->abv = nnodes[seg];
				(ptree+nnodes[seg])->time = t;

			}
			else {
				(pseg+tseg)->desc = ( yes1 ?
				   (chrom[c1].pseg + seg1)->desc :
				  (chrom[c2].pseg + seg2)->desc);
				}
			}
		}
	nlinks -= links(c1);
	cleft -= 1.0 - pow(pc, (double)links(c1));
	free(chrom[c1].pseg) ;
	if( tseg < 0 ) {
		free(pseg) ;
		chrom[c1].pseg = chrom[nchrom-1].pseg;
		chrom[c1].nseg = chrom[nchrom-1].nseg;
		chrom[c1].pop = chrom[nchrom-1].pop ;
		if( c2 == nchrom-1 ) c2 = c1;
		nchrom--;
		}
	else {
		if( !(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
			perror(" realloc error. ca1");
		chrom[c1].pseg = pseg;
		chrom[c1].nseg = tseg + 1 ;
		nlinks += links(c1);

	   	cleft += 1.0 - pow(pc, (double)links(c1));
		}
	nlinks -= links(c2);
	cleft -= 1.0 - pow(pc, (double)links(c2));
	free(chrom[c2].pseg) ;
	chrom[c2].pseg = chrom[nchrom-1].pseg;
	chrom[c2].nseg = chrom[nchrom-1].nseg;
	chrom[c2].pop = chrom[nchrom-1].pop ;
	nchrom--;
	if(tseg<0) return( 2 );  /* decrease of nchrom is two */
	else return( 1 ) ;
}

/*** Isseg: Does chromosome c contain the segment on seglst which starts at
	    start? *psg is the segment of chrom[c] at which one is to begin 
	    looking.  **/

	int
isseg(start, c, psg)
	int start, c, *psg;
{
	int ns;
	struct seg *pseg;

	ns = chrom[c].nseg;
	pseg = chrom[c].pseg;

/*  changed order of test conditions in following line on 6 Dec 2004 */
	for(  ; ((*psg) < ns ) && ( (pseg+(*psg))->beg <= start ) ; ++(*psg) )
		if( (pseg+(*psg))->end >= start ) return(1);
	return(0);
}
	


	int
pick2_chrom(pop,config,pc1,pc2)
	int pop, *pc1, *pc2, config[];
{
	int c1, c2, cs,cb,i, count;
	
	pick2(config[pop],&c1,&c2);
	cs = (c1>c2) ? c2 : c1;
	cb = (c1>c2) ? c1 : c2 ;
	i=count=0;
	for(;;){
		while( chrom[i].pop != pop ) i++;
		if( count == cs ) break;
		count++;
		i++;
		}
	*pc1 = i;
	i++;
	count++;
	for(;;){
		while( chrom[i].pop != pop ) i++;
		if( count == cb ) break;
		count++;
		i++;
		}
	*pc2 = i ;
}
	
	

/****  links(c): returns the number of links between beginning and end of chrom **/

	int
links(c)
	int c;
{
	int ns;

	ns = chrom[c].nseg - 1 ;

	return( (chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}

