/**************************************************************************************
20 Sept, 1995
burst.c, modified from pois_bst.c

Modified the code to run as an .exe program. kgt

The user will be prompted to provide a text (ASCII) file with
spike times sorted in ascending order in a single column
(a <return> between each time).

A sample file called spike.txt is provided as an example.

The user will be prompted to give a starting time and an ending time of
the spike train to be analyzed.

Also, the significance level needs to be specified.  All 'bursts' with a
surprise index corresponding to a significance level above the level given by
the user will not be listed.

Activation times will be listed at the users request.




11 May, 1995
pois_bst.c
The code for the Poisson spike train analysis that you requested is below.
The rationale of our algorithm is presented in a manuscript (Hanes, Thompson and Schall
(1995) Relationship of presaccadic activity in frontal eye field and supplementary
eye field to saccade initiation in macaque: Poisson spike train analysis.
Eperimental Brain Research, 103:85-96).
	The following algorithm is an extension of the algorithm utilized in Hanes
et al, 1995 that can detect multiple bursts of activity within a spike train.
The code for the multiple burst detector is below.  We would be very interested
in your comments about the performance of the algorithm on your data.


Doug Hanes, Kirk Thompson, & Jeff Schall
Dept of Psychology
Vanderbilt University

email -- kirk.g.thompson@vanderbilt.edu

******/

/*  pois_bst.c

	Poisson Burst Detector

	Algorithm adapted from Legendy and Salcman, 1985

	This algorithm searches through a spike train and determines periods of
	"bursts" that maximally deviate from the expected number of spikes given
	the random Poisson distribution calculated from the average rate of
	spike train.

	It also defines times around these "bursts" of activity that significantly
	differ from the expected number of spikes.

	These times of activity are stored in struct activtn[] as four times:
		.boa = beginning of significant activation for the current "burst"
		.bob = beginning of the "burst"
		.eob = end of the "burst"
		.eoa = end of significant activation for the current "burst"
	Also,
		.surprise = surprise index of the current burst

	Three variables are required by the user:

	int sp_time[]           an array of spike times
	int burst_spike         the number of spikes in the time period that you are looking for bursts
	int burst_time          the duration of the time period in which you are looking for bursts

	For the purpose of demonstration, values are given to these three variables below.
*/

#include <math.h>       /* math functions */
#include <stdlib.h>		/* standard functions. max(), min() */
#include <stdio.h>
#include <curses.h>

#define MAXSPIKES      1000
#define MAXACTIVATION  15
#define MAXPROB        0.99999999999999
#define NULL		   0
#define NO			   0
#define YES			   1

/* Variables the multiple burst detector requires from user */
char spikefile[15];
FILE *in;
int starttime = 0;
int endtime = 0;


float sp_time[MAXSPIKES] = {NULL};/*{59,83,115,158,347,460,464,472,519,535,551,559,566,
							571,575,584,587,593,596,604,607,609,612,614,619,
							620,622,624,625,628,631,632,633,638,639,650,657,
							681,715,855,896,924,988};    array holding spike times */

int burst_spike = 0;   /* 43;                the number of spikes in the time period that you are looking for bursts */
int burst_time = 0;    /*1000;				 the duration of the time period in which you are looking for bursts */


/* BURST DETECTOR VARIABLES  */
double isi = 0.;					/* inter-spike-interval for spikes being tested	*/
double isi_rate = 0.;       		/* rate for spikes being tested			*/
double sup_max = 0.;        		/* maximum surprise value			*/
float prob = 0.;					/* prob level (suprise index) the burst must be above to be determined a burst. derived from signif_level*/
float avg_rate = 0;        			/* avg rate during baseline period 100 ms after fixation to 100 ms after saccade */

float sup_new = 0.;        			/* surprise value for current spikes being testd*/
float signif_level = .0;			/* probability for burst detection. input by the user 		*/
float spike_time = 0;					/* spike times as they are read from the input file */
float previous_spike_time = 0;		/* the previous spike time as they are read from the input file */
int MinBurstSpike = 2;				/* minimum spikes allowable in a burst */
int min_isi_between_burst = 8;		/* flag for the min isi b/w/ bursts where they are considered as indistinguishable */
int num_isi = 0;            		/* # of isi being tested			*/
int num_spike = 0;          		/* # of spikes being tested num_isi + 1		*/
int t = 0;                  		/* time difference of spike interval being tested*/
int burst = NO;             		/* was a burst detected				*/
int on = NO;                   		/* flag set when the surprise value first goes above signif prob level when finding boa and eoa */
int first_spike_before_burst = 0;   /* used in determining eoa */
int first_spike_after_burst = 0;    /* the spike at which to begin looking for bursts */
int number_of_activations = 0;		/* the number of activations found */
int act_num = 0;                    /* the current activation number */
int MaxNumberExtraSpikes = 10; 		/* Number of spikes after eob that burst detector searches for increased SI */
int MaxExtraTime = 30;				/* time after eob that burst detector searches for increased SI */

/* structure for activation and burst variables */
struct {
	int boa;
	int eoa;
	int bob;
	int eob;
	double surprise;
} activtn[MAXACTIVATION];

/* FUNCTION DECLARATIONS */
// int multiple_burst(void);
double surprise(int,int,double);
void find_activation(void);

int calcSurprise(void)
{
    double S;
    S = surprise (100, 10, 0.01);
    printf ("Surprise index: %2.7f\n",S);

    S = surprise (100, 10, 0.04);
    printf ("Surprise index: %2.7f\n",S);

    S = surprise (100, 50, 0.04);
    printf ("Surprise index: %2.7f\n",S);

    S = surprise (100, 7, 0.1);
    printf ("Surprise index: %2.7f\n",S);

    S = surprise (5, 1, 1);
    printf ("Surprise index: %2.7f\n",S);

}


// void multiple_burst() /*Poisson burst detector  */
int main(void)
{
	register i = 0;
	register m = 0;
	int CurrentNumberExtraSpikes = 0; /* Current number of spikes after eob */
	int CurrentExtraTime = 0;         /* Current time after eob */
	int sum_burst_time = 0;     /* total amount of time in each trial within a burst period */
	int sum_burst_spike = 0;    /* total amount of spikes in each trial within a burst period */
	int bob_anchor = 0;			/* number of spikes for anchor point for determining bob */
	int eob_anchor = 0;			/* number of spikes for anchor point for determining eob */
	int anchor_time = 0;		/* variable time for bob_anchor and eob_anchor determination */
	int old_bob = -1;			/* used in while loop for determining if bob changes with next iteration of burst detector */
	int old_eob = -1;			/* used in while loop for determining if eob changes with next iteration of burst detector */
	int iterate = 0;			/* flag that determines whether burst detector needs to go into the next iteration */
	int iterations = 0;			/* number of iterations on current burst */
	int max_eob = 0;			/* max eob value on current burst with multiple iterations */
	int min_bob = 5000;			/* (5000 is a bogus value) min bob value on current burst with multiple iterations */
	int current_act_num = 0;	/* current act_num value */

    
    calcSurprise();


	/* inputs file of spike times */
	printf("\nSpike file (with spikes times sorted in ascending order):");
	gets(spikefile);
	printf("opening spike file: %s", spikefile);
	if ((in = fopen(spikefile, "rt")) == NULL){
		printf("\nerror: Cannot find spike file:  %s", spikefile);
		return 1;
	}

	/* inputs starting and ending times of spike train */
	printf("\n\nSpike train start time:");
	scanf ("%i",&starttime);
	printf("Spike train end time:");
	scanf ("%i",&endtime);


	/* inputs spike times from the spike file to sp_time[] */
	for (i=0; i <= MAXSPIKES; i++){
		fscanf(in, "%f\n", &spike_time);

		if (feof(in))
			break;

		if ((spike_time > starttime) && (spike_time < endtime)){

			/* if you are collapsing across trials it is possible that spikes
			** will occur simultaneously so to prevent an error .001 is added
			** to spikes until they are sorted in ascending order with different
			** times.
			*/
			while (spike_time == sp_time[burst_spike])
				spike_time += .001;

			sp_time[burst_spike] = spike_time;
			burst_spike++;
		}

        printf("spike_time: %2.7f\n",spike_time);

	}
	/* closes spike file */
	fclose(in);

	burst_time = endtime - starttime;
	if (burst_time<=0){
		printf ("\nerror: Your end time must be greater than your start time.");
		return 1;
	}

	/*
	** To treat each burst as equally as possible a consistent time before bob
	** and after eob is used as the anchor point of the interval 't' to find
	** eob and bob during multiple iterations.
	*/
	anchor_time = 50;

	/* inputs significance level for burst to be counted as a burst */
	printf("\nsignificance level: p < ");
	scanf ("%f",&signif_level);

	/*
	** Use natural ln() of probability, e.g.
	** ln(0.01) = 4.61, i.e., 0.01 probability
	*/
	prob = -log ((double)signif_level);
    printf("Probability: %2.2f\n",prob);
	/*
	**	determine average rate
	*/
	avg_rate = (float)(burst_spike)/(double)(burst_time);
    printf("Burst Spike: %2.7f, Burst Time: %2.7f, Average Rate: %2.7f",(float)burst_spike, (double)burst_time,avg_rate);
	/*
	** the starting point in the spike array in which to look for bursts
	*/
	first_spike_after_burst = 0;
	i = 0;

	/*
	** if there are less than MinBurstSpike spikes in interval,
	** then no need to run burst detector
	*/
	if (burst_spike <= MinBurstSpike)
		return(NO);

	while((first_spike_after_burst <= burst_spike-2) && (i < burst_spike-2)) {

		/*
		** Takes first 2 ISIs and finds average. If the average rate for
		** these 3 spikes is greater than the overall rate (avg_rate), the burst
		** is said to begin here.  If not, then the rate for the next 2 isi's
		** is calculated.
		*/

		for (i = first_spike_after_burst; i < (burst_spike-2); i++) {

			isi = 0;
			isi_rate = 0;
			isi = (sp_time[i + 2] - sp_time[i]);
			isi /= 3;
			isi_rate = 1 / (double)isi;
			if (isi_rate >= avg_rate) {

				/* initialize variables */
				activtn[act_num].bob = 0;
				activtn[act_num].eob = 0;
				burst = NO;
				sup_new = 0;
				sup_max = 0;
				on = NO;
				CurrentNumberExtraSpikes = 0;
				CurrentExtraTime = 0;
				eob_anchor = 0;

				/*
				** The beginning and end of each putative burst is calculated
				** multiple times until their time values are consistent for
				** two iterations
				*/
				/* do-while loop for multiple iterations of burst detector */
				do {
					/* initialize variables */
					sup_max = 0;
					sup_new = 0;
					old_bob = activtn[act_num].bob;
					old_eob = activtn[act_num].eob;

					/* finds the number of spikes in anchor interval before current starting point */
					if (!activtn[act_num].bob){
						for (m = i-1; m >= 0; m--){
							if (sp_time[m] >= (sp_time[i] - anchor_time))
								eob_anchor++;
							else
								break;
						}
					}else{
						/* determines number of spikes after first iteration when there is a current bob */
						for (m = activtn[act_num].bob-1; m >= 0; m--){
							if (sp_time[m] >= (sp_time[activtn[act_num].bob] - anchor_time))
								eob_anchor++;
							else
								break;
						}

					}

					/*
					** This loop adds spikes for until the surprise index
					** decreases for a user defined time or number of spikes.
					** Finds End of Burst (eob)
					*/
					for (num_isi = 2; num_isi < burst_spike - i; num_isi++) {

						if (!activtn[act_num].bob){
							/* first iteration */
							t = sp_time[i + num_isi] - (sp_time[i] - anchor_time);
						}else{
							if (activtn[act_num].bob + num_isi < burst_spike)
							/* after first iteration */
								t = sp_time[activtn[act_num].bob + num_isi] - (sp_time[activtn[act_num].bob] - anchor_time);
							else
								break;
						}

						num_spike = eob_anchor + num_isi;
						sup_new = surprise(t, num_spike, avg_rate);

						if (sup_new >= sup_max) {

							/* first iteration */
							if (!activtn[act_num].bob)
								activtn[act_num].eob = i + num_isi;
							/* after first iteration */
							else
								activtn[act_num].eob = activtn[act_num].bob + num_isi;

							sup_max = sup_new;
							CurrentNumberExtraSpikes=0;

						} else {
							/*
							** MaxNumberExtraSpikes is the number of spikes after
							** eob in which the surprise index can increase above the
							** current sup_max.  Currently set to 10.
							** MaxExtraTime is the time after eob in which the
							** surprise index can increase above current sup_max.
							** Currently set to 30ms; approximates EPSP decay.
							** These numbers need to be determined by the user.
							*/
							CurrentNumberExtraSpikes++;

							/* first iteration */
							if (!activtn[act_num].bob)
								CurrentExtraTime = sp_time[i+num_isi] - sp_time[activtn[act_num].eob];
							/* after first iteration */
							else
								CurrentExtraTime = sp_time[activtn[act_num].bob + num_isi] - sp_time[activtn[act_num].eob];

							if (CurrentNumberExtraSpikes >= MaxNumberExtraSpikes ||
							  CurrentExtraTime >= MaxExtraTime)
								break;
						}

					}


					/* initialize variables */
					CurrentNumberExtraSpikes = 0;
					sup_new = 0;
					sup_max = 0;
					on = NO;
					bob_anchor = 0;
					eob_anchor = 0;

					/* resets min_bob at the 10th iteration */
					if (iterate && iterations == 10)
						min_bob = activtn[act_num].bob;

					/* counts number of spikes in anchor_time after eob */
					for (m = activtn[act_num].eob+1; m < burst_spike; m++){
						if (sp_time[m] <= (sp_time[activtn[act_num].eob] + anchor_time))
							bob_anchor++;
						else
							break;
					}


					/*
					** This loop indexes backward in time to the end of the previous
					** burst and finds begin spike when surprise is max.
					** Finds Beginning Of Burst (bob).
					*/
					for (m = activtn[act_num].eob-1; m >= first_spike_after_burst; m--){
						t = anchor_time + sp_time[activtn[act_num].eob] - sp_time[m];
						num_spike = bob_anchor + activtn[act_num].eob - m;
						sup_new = surprise(t,num_spike,avg_rate);


						if (sup_new >= sup_max) {
							sup_max = sup_new;
							activtn[act_num].bob = m;
						}
					}

					if ((activtn[act_num].bob == old_bob) && (activtn[act_num].eob == old_eob))
						iterate = NO;
					else
						iterate = YES;

					/*
					** if caught in infinite loop on the same burst then only 20
					** iterations will be allowed and the max_bob and the min_eob
					** determined during these iterations will become bob and eob
					** for this burst.
					*/
					current_act_num = max(act_num, current_act_num);
					if (current_act_num == act_num){
						max_eob = max(activtn[act_num].eob,max_eob);
						min_bob = min(activtn[act_num].bob,min_bob);
						iterations++;
					}
					if (iterations > 20){
						activtn[act_num].eob = max_eob;
						activtn[act_num].bob = min_bob;
						iterate = NO;
						iterations = 0;
					}

				} while (iterate);

				/* finds surprise index of the current burst */
				sup_max = surprise(sp_time[activtn[act_num].eob]-sp_time[activtn[act_num].bob],(activtn[act_num].eob-activtn[act_num].bob+1),avg_rate);

				/* resets first_spike_after_burst to look for the next burst */
				first_spike_after_burst = activtn[act_num].eob + 1;

                printf("Through Loop: (%f,%f) Surprise: %2.7f\n",sp_time[activtn[act_num].bob], sp_time[activtn[act_num].eob], sup_max);

				/* test to check if burst meets minimum spike & probability criteria */
				if ((sup_max > prob) && ((activtn[act_num].eob - activtn[act_num].bob) >= MinBurstSpike)){
					number_of_activations++;
					burst = YES;
					activtn[act_num].surprise = sup_max;
					act_num++;
				}

				/* if no burst was found, reset variables */
				if (!burst){
					activtn[act_num].bob = activtn[act_num].eob = 0;
					max_eob = 0;
					min_bob = 5000;				/* a bogus value */
				}

				/*
				** resets iteration variables when the next activation
				** period has been found
				*/
				if (current_act_num < act_num){
					iterations = 0;
					max_eob = 0;
					min_bob = 5000;				/* a bogus value */
					current_act_num = 0;
				}

				break;
			}
		}
	}


	/*
	** To remove the influence of bursts in determining significant activation
	** times, we remove bursts when calculating average rate for activation times
	*/
	for (act_num = 0; act_num < number_of_activations; act_num++){
		sum_burst_time += sp_time[activtn[act_num].eob] - sp_time[activtn[act_num].bob];
		sum_burst_spike += activtn[act_num].eob - activtn[act_num].bob;
	}
	/* recalculate average rate without bursts */
	avg_rate = (float)(burst_spike-sum_burst_spike)/(double)(burst_time-sum_burst_time);

	/* if bursts are found then find beginning of activation (boa) and end of activation (boa) */
	if (act_num) find_activation();

	/*
	** if times of activation overlap then the furthest activation times are used
	*/
	for (i = 0; i < number_of_activations - 1; i++) {
		for (act_num=0; act_num<number_of_activations-1; act_num++) {
			if(activtn[act_num].eoa >= activtn[act_num+1].boa) {
				activtn[act_num].eoa = activtn[act_num+1].eoa;
				activtn[act_num+1].boa = activtn[act_num].boa;
			}
		}
	}
	printf ("\n");
	printf ("\n");

	if (number_of_activations){
		for (i = 0; i < number_of_activations; i++){
			printf ("Burst begin: %2.1f ms  ",sp_time[activtn[i].bob]);
			printf ("Burst end: %2.1f ms  ",sp_time[activtn[i].eob]);
			printf ("Surprise index: %2.2f",activtn[i].surprise);
			printf ("\n");
		}
		printf ("\n");
		printf ("List activation times? y/n  ");
		if (getch() == 'y'){
			printf ("\n");
			printf ("\n");
			for (i = 0; i < number_of_activations; i++){
				if (activtn[i].boa!= activtn[i-1].boa){
					printf ("Activation begin: %2.1f ms  ",sp_time[activtn[i].boa]);
					printf ("Activation end: %2.1f ms  ",sp_time[activtn[i].eoa]);
					printf ("\n");
				}
			}
		}else{
			printf ("\n");
		}
	}else{
		printf ("\n No bursts with a significance level less than %1.5f.",signif_level);
	}

	printf("\n\nPoisson spike train analysis algorithm by Kirk Thompson & Doug Hanes, 1995");

	return(NO);
}




inline int max(int a, int b)
{
        return (a > b) ? a : b;
}



inline int min(int a, int b)
{
        return (a < b) ? a : b;
}



/***************************************************************************
**                        Determine's boa and eoa                         **
***************************************************************************/
void find_activation()
{
	register i,m;

	first_spike_before_burst = 0;

	for (act_num = 0; act_num < number_of_activations; act_num++){

		/*
		** finds spike to begin process of EOA determination
		*/
		if(act_num > 0 && number_of_activations)
			first_spike_after_burst = activtn[act_num-1].eob;
		else
			first_spike_after_burst = 0;

		/*
		** This loop goes from the beginning of burst to the end
		** of the previous burst and decrements spikes from the beginning
		** of this interval and finds the first spike where the surprise
		** value goes above 'prob', this defines the
		** Beginning Of Activation (boa) for this burst.
		*/
		on = NO;
		for (m = first_spike_after_burst; m < activtn[act_num].bob; m++){
			t = abs(sp_time[activtn[act_num].bob] - sp_time[m]);
			num_spike = activtn[act_num].bob - m + 1;
			sup_new = surprise(t,num_spike,avg_rate);

			if (sup_new >= prob) {
				activtn[act_num].boa = m;
				on = YES;
				break;
			}
		}

		/* if no boa is found then boa is made equal to bob */
		if (!on) activtn[act_num].boa = activtn[act_num].bob;
		on = NO;

		/*
		** finds last spike for process of determining EOA
		*/
		if(activtn[act_num+1].bob){
			first_spike_before_burst = activtn[act_num+1].bob;
		}else{
			first_spike_before_burst = burst_spike-1;
		}

		/*
		** This loop goes from the end of burst to the beginning
		** of the next burst and decrements spikes from the end
		** of this interval and finds the first spike where the surprise
		** value goes above 'prob', this defines the
		** End Of Activation (eoa) for this burst.
		*/
		for (m = first_spike_before_burst; m > activtn[act_num].eob; m--) {
			t = abs(sp_time[m] - sp_time[activtn[act_num].eob]);
			num_spike = m - activtn[act_num].eob + 1;
			sup_new = surprise(t,num_spike,avg_rate);

			if (sup_new >= prob) {
				activtn[act_num].eoa = m;
				on = YES;
				break;
			}
		}
		if (!on)
			activtn[act_num].eoa = activtn[act_num].eob;
		on = NO;
	}
}

/****************************************************************************/

double surprise (t, num_spike, avg_rate)
int num_spike;
int t;
double avg_rate;

{
	register j;
	double rt = 0.;
	double sum = 0.;
	double p = 0.;
	double new_prob = 0.;
	double old_prob = 0.;
	double sup_new;



	/*
	** does prob summation
	** for 0 to n spikes
	*/
	rt = (double)avg_rate * t;
	sum = old_prob = exp(-rt);

	for (j = 1; j < num_spike; j++) {
		new_prob = (double)(rt * old_prob) / j;
		sum += new_prob;
		/*
		** If sum gets sufficiently close to 1.000...,
		** break out of the loop
		*/
		if (sum >= 1){
			sum = MAXPROB;
			break;
						}
		old_prob = new_prob;
	}

	/*
	** Subtract obtained probability from 1.0 to determine tail
	** of distribution
	*/
	p = 1.0 - sum;

	/* the next line is for when p = 1 (log(1) = 0); the algorithm will sometimes stop searching for bob and eob */
	if (p>=0.9999999999999)
		p = 0.9999999999999;

	sup_new = -log((double)p);
	return (double)sup_new;
}
