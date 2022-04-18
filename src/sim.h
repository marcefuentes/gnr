#include <stdbool.h>

#define CONTINUOUS_V 5
#define BINS 64
#define BOOLEAN_V 2

// Structures

struct ptype
{
	int		time;
	double		sumc[CONTINUOUS_V][BINS], sumc2[CONTINUOUS_V][BINS];
	double		sumBD[CONTINUOUS_V], sumBD2[CONTINUOUS_V];
	double		summedian[CONTINUOUS_V], summedian2[CONTINUOUS_V];
	double		sumb[BOOLEAN_V], sumb2[BOOLEAN_V];
};

struct pruntype
{
	int		time;
	double		frc[CONTINUOUS_V][BINS];
	double		median[CONTINUOUS_V];
	double		frb[BOOLEAN_V];
};

struct itype
{
	double		a2Default; 	
	double		a2Decided;		// a2 for next round
	double		a2Seen;			// a2 in present round
	double		wCumulative; 	
	double		ChooseGrain;
	double		MimicGrain;
	double		cost;			// Information costs
	bool		isRecruit;		// It can't be killed. It isn't known to, and doesn't know, group mates
	bool		chose_partner;		// It has chosen a new partner
	bool		changed_a2;		// It has changed a2Decided
	struct itype	*partner;
	struct itype	*newpartner;
};

struct rtype
{
	double		randomwc;
	double		a2Default;
	double		ChooseGrain;
	double		MimicGrain;
	struct rtype	*next;
};

// Functions

void		decide_a2		(struct itype *i, struct itype *i_last, double amax);
void		decide_a2_ir		(struct itype *i, struct itype *i_last, double amax);
void		shuffle_partners	(struct itype *i, struct itype *i_last, int groupsize);
void		choose_partner		(struct itype *i, struct itype *i_last, int groupsize);
struct rtype	*create_recruits	(int deaths, double wc);
void		free_recruits		(struct rtype *recruit);
void		kill			(struct rtype *recruit, struct itype *i_first, int n);
void		stats_period		(struct itype *i_first, struct itype *i_last, struct pruntype *prun, int n, double amin, double amax, double wmax);
void		stats_end		(struct pruntype *prun_first, struct pruntype *prun_last, struct ptype *p_first);
void		stats_runs		(struct ptype *p, struct ptype *p_last, int runs);
void		write_headers		(char *filename, char *header1, char *header2);
void		write_stats		(char *filename, float factor1, float factor2, int groupsize, struct ptype *p, struct ptype *p_last);
void		write_time_elapsed	(char *filename, float time_elapsed);

// Functions in main.c needed by functions in other files

double		calculate_cost		(double choose, double mimic);		// gChooseCost, gMimicCost
