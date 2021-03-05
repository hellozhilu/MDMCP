#ifndef STRUCT_H_
#define STRUCT_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
using namespace std;

#define DEBUG_SWT 0
#define MIN(a,b) ((a)<(b)?(a):(b))
#define PARENT(idx) ((idx-1)/2)
#define LEFT(idx) (2*idx+1)
#define RIGHT(idx) (2*idx+2)
#define EMPTY_IDX 0

extern const int MAX_CHAR ;
extern const int MAX_VAL ;
extern const float CONST_E;

/***********************Definition of Graph*********************/
extern int nnode;
extern int **matrix;

/****************Declaration of all the parameters**************/
//configurations
extern char param_filename[1000];
extern int param_knownbest ;
extern int param_time;		     //the max time for memetic procedure, unit: seconds
extern int param_seed;
extern int param_max_generations;//not used
//5 parameters
extern int param_sizefactor;     //theta_size of SA
extern double param_tempfactor;  //theta_cool of SA
extern double param_minpercent;  //theta_minper of SA
extern double param_shrink;      //shrink percent in crossover
extern int param_pool_size;      //size of population

/***********************Partition data**************************/
extern FILE *fout;

typedef struct{
	int *ppos;    /*ppos[pid] is the position of partition "pid" in "pbkt"*/
	int *pbkt;    /*Elements from pbkt[0] to pbkt[pbkt_size-1] include all the partition ids*/
	int pbkt_size;/*The number of partitions*/
	int *pcnt;    /*The number of vertices in each partition*/
	int *pvertex; /*The partition of each vertex, for example, v is in partition pvertex[v]*/
}PartitionType;

extern void printPartition(PartitionType *ppt, int gnnode, FILE *fout);
extern PartitionType* allocatePartitionData(int gnnode);
extern void buildPartition(PartitionType *ppt, int *vpart, int gnnode);
extern void copyPartition(PartitionType *dest, PartitionType *source, int gnnode);
extern void disposePartition(PartitionType *ppt);
extern void updatePartition(PartitionType *ppt, int v, int target);
extern int calculateDistance(int *p1, int *p2, int gnnode);
extern void generateRandList(int *randlist, int len);
extern void swapAry(int *ary, int idx1, int idx2);
extern int size_inter_section(vector<int> *v1, vector<int> *v2);
extern int calculateMaxMatch(int *p1, int n_part1, int *p2, int n_part2);

/***************************Statistic***************************/
typedef struct ST_Stats{
	PartitionType *best_partition;

	int best_val; 			//largest objective value
	double best_foundtime;  //The first time best_obj was found
	int best_generation;    //The iteration when best_obj was found

}BestSolutionInfo;

extern void clearResult(BestSolutionInfo *sts);
//extern void printResult(BestSolutionInfo *sts, int index);

/****************Simulated annealing local search***************/
/*SA internal data structure*/
typedef struct SA_RT_Data{
	int **gmatrix;
	int gnnode;

	PartitionType *ppt;      /*current solution*/
	int fcurrent; 	         /*The current objective function*/

	PartitionType *ppt_best; /*best solution*/
	int fbest;	             /*The ever found best solution*/

	int itr;
	int totalIter;

	int **gammatbl;
}SA_RT_Data;

extern SA_RT_Data *lsdata;
extern int generateInitSolution(int **gmatrix, int gnnode, PartitionType *ppt);
extern void buildCurGamma(int *pvertex);
extern void updateCurGamma(int u, int src, int dest);
extern void setLSEnvironment(int **graph_mat, int gnnode);
extern void setLSStart(PartitionType *start_sol, int start_obj_value);
extern void disposeLSEnvironment();
//static int decideTarget(int dest);
//static int changeCurSolution(int u, int dest);
extern void pureDescent();
extern double calibrateTemp();
extern void annealingSearch(double initTemp, clock_t start_time, double &best_time);

/**************************memetic run**************************/
extern void memeticRun(BestSolutionInfo *frt, int *totoal_gen, int iter_limit, int max_seconds);

#endif /* STRUCT_H_ */
