#include <stdbool.h>

#define CONTINUOUS_V 5
#define BINS 64
#define BOOLEAN_V 2

// Structures

struct ptype
{
	int		time;
	double		alpha, logES, Given, wmax;
	double		summean[CONTINUOUS_V], summean2[CONTINUOUS_V];
	double		sumsd[CONTINUOUS_V], sumsd2[CONTINUOUS_V];
	double		sumc[CONTINUOUS_V][BINS], sumc2[CONTINUOUS_V][BINS];
	double		sumBD[CONTINUOUS_V], sumBD2[CONTINUOUS_V];
	double		summedian[CONTINUOUS_V], summedian2[CONTINUOUS_V];
	double		sumiqr[CONTINUOUS_V], sumiqr2[CONTINUOUS_V];
	double		sumb[BOOLEAN_V], sumb2[BOOLEAN_V];
	double		summeanTS, summean2TS;
	double		sumsdTS, sumsd2TS;
};

struct pruntype
{
	int		time;
	double		alpha, logES, Given, wmax;
	double		mean[CONTINUOUS_V];
	double		sd[CONTINUOUS_V];
	double		frc[CONTINUOUS_V][BINS];
	double		median[CONTINUOUS_V];
	double		iqr[CONTINUOUS_V];
	double		frb[BOOLEAN_V];
	double		meanTS;
	double		sdTS;
};

struct itype
{
	double		a2Default; 	
	double		a2Decided;		// a2 for next round
	double		a2Seen;			// a2 in present round
	double		a2SeenSum;		// sum of a2Seen since birth
	double		wCumulative; 	
	double		ChooseGrain;
	double		MimicGrain;
	double		cost;			// Information costs
	int		age;			// It can't be killed. It isn't known to, and doesn't know, group mates
	bool		chose_partner;		// It has chosen a new partner
	bool		changed_a2;		// It has changed a2Decided
	struct itype	*oldpartner;
	struct itype	*partner;
};

struct rtype
{
	double		randomwc;
	double		a2Default;
	double		ChooseGrain;
	double		MimicGrain;
	double		cost;
	struct rtype	*next;
};

// Functions

void		decide_a2		(struct itype *i, struct itype *i_last, double amax, int indirectr, int discrete);
void		shuffle_partners	(struct itype *i, struct itype *i_last, int groupsize);
void		choose_partner		(struct itype *i, struct itype *i_last, int groupsize);
struct rtype	*create_recruits	(int deaths, double wc);
void		free_recruits		(struct rtype *recruit);
void		kill			(struct rtype *recruit, struct itype *i_first, int n, int discrete);
void		stats_period		(struct itype *i, struct itype *i_last, struct pruntype *prun, int n, double amin, double amax);
void		stats_end		(struct pruntype *prun, struct pruntype *prun_last, struct ptype *p);
void		stats_runs		(struct ptype *p, struct ptype *p_last, int runs);
void		write_headers_csv	(char *filename);
void		write_headers_i		(char *filename);
void		write_headers_frq	(char *filename);
void		write_stats_csv		(char *filename, struct ptype *p, struct ptype *p_last);
void		write_stats_frq		(char *filename, struct ptype *p, struct ptype *p_last);
void		write_i			(char *filename, float alpha, float logES, float Given, struct itype *i, struct itype *i_last);
void		write_time_elapsed	(char *filename, float time_elapsed);
void		file_write_error	(char *filename);
double		stdev			(double sum, double sum2, int n);
