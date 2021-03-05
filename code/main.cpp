/* COPYRIGHT
 * Zhi Lu, Yi Zhou, Jin-Kao Hao,
 * A Hybrid Evolutionary Algorithm for the Clique Partitioning Problem,
 * IEEE Transactions on Cybernetics, January 2021
 */

#include "defines.h"


int nnode;		//The number of the nodes
int nedge;		//The number of the edges
int **matrix;	//Adjacent matrix of the input graph

/*Default configuration*/
char param_filename[1000] = "./instance/rand500-100.txt";
int param_knownbest = 309125;
int param_time = 500;		       //the max time for memetic procedure, unit: second
int param_seed = 123456;
int param_max_generations = 100000;//not used
//five parameters
int param_sizefactor = 8;          //theta_size of SA
double param_tempfactor = 0.96;    //theta_cool of SA
double param_minpercent = 1.0;     //theta_minper of SA
double param_shrink = 0.6;		   //shrink ratio in crossover
int param_pool_size = 10;          //size of population

double totaltime;
int totalgen;
BestSolutionInfo finalBest;
FILE *fout = NULL;


/* print error information */
void showUsage(){
	cerr << "usage: [-f <file path>] [-t <run time>] [-g <seed>] [-v <best solution>] "
			"[-x <max generation>] [-b <theta size>] [-c <theta cool>] [-d <theta minper>] "
			"[-s <shrink ratio>] [-p <pool size>]\n" << endl;
}

/* You can input parameters in arbitrary order */
void readParameters(int argc, char **argv){
	for (int i = 1; i < argc; i+=2){
		if (argv[i][0] != '-' || argv[i][2] != 0){
			showUsage();
			exit(0);
		}else if(argv[i][1] == 'f'){//file path
			strncpy(param_filename, argv[i+1], 1000);
		}else if (argv[i][1] == 't'){//run time
			param_time = atoi(argv[i+1]);
		}else if(argv[i][1] == 'g'){//seed
			param_seed = atoi(argv[i+1]);
		}else if(argv[i][1] == 'v'){//best solution
			param_knownbest = atoi(argv[i+1]);
		}else if(argv[i][1] == 'x'){//max generation (not used)
			param_max_generations = atoi(argv[i+1]);
		}else if(argv[i][1] == 'b'){//theta size
			param_sizefactor = atoi(argv[i+1]);
		}else if(argv[i][1] == 'c'){//theta cool
			param_tempfactor = atof(argv[i+1]);
		}else if(argv[i][1] == 'd'){//theta minper
			param_minpercent = atof(argv[i+1]);
		}else if(argv[i][1] == 's'){//shrink ratio
			param_shrink = atof(argv[i+1]);
		}else if(argv[i][1] == 'p'){//pop size
			param_pool_size = atoi(argv[i+1]);
		}
	}

	/*check parameters*/
	if (strlen(param_filename) == 0){
		cerr << "No input data" << endl;
		exit(1);
	}
}

/* for CPP graph format
 * (for other graph format, you can translate it in CPP format)*/
void loadCompletedGraph(char *filename){
	ifstream fin;
	fin.open(filename);
	if (fin.fail()){
		cerr << "Can not open file " << filename << endl;
		exit(0);
	}
	fin >> nnode;
	if (fin.eof()){
		cerr << "Empty file" << filename << endl;
		exit(0);
	}
	matrix = new int*[nnode];
	for (int i = 0; i < nnode; i++)
		matrix[i] = new int[nnode];
	int val;
	for (int ni = 0; ni < nnode; ni++){
		for (int nj = ni; nj < nnode; nj++){
			fin >> val;
			matrix[ni][nj] = matrix[nj][ni] = -val;
		}
	}
	nedge = nnode * (nnode-1) /2;
}

void clearResult(BestSolutionInfo *sts){
	sts->best_generation = 0;
	sts->best_val = -MAX_VAL;
	sts->best_foundtime = 0.0;
}

FILE* setupRecordFile(){
	/* creat file in current direct */
	char path_cwd[FILENAME_MAX];
	char *graph_name = basename(param_filename);
	char file_name[FILENAME_MAX];

	getcwd(path_cwd, FILENAME_MAX);
	sprintf(file_name, "%s/rec/%s_%5d.rec", path_cwd, graph_name, param_seed);
	cout << file_name << endl;

	FILE *f = fopen(file_name, "w+");
	if (f == NULL){
		return 0;
	}
	return f;
}

void reportResult(){
	printf("\n");
	printf("$seed=%d\n", param_seed);
	printf("@solver=CPP\n");
	printf("@para_file=%s\n", param_filename);
	printf("@para_kbv=%d\n", param_knownbest);
	printf("@para_seconds=%d\n", param_time);
	printf("#vnum=%d\n", nnode);
	printf("#objv=%d\n", finalBest.best_val);
	printf("#bestiter=%d\n",finalBest.best_generation);
	printf("#besttime=%.3f\n", finalBest.best_foundtime);
	printf("#totoaltime=%.2f\n", totaltime);
	printf("#totoaliter=%d\n", totalgen);
	printf("#bestpartition=%d\n", finalBest.best_partition->pbkt_size-1);
	printf("Solution:");
	for (int i = 0; i < nnode; i++){
		printf("%d ", (finalBest.best_partition)->pvertex[i]);
	}
	printf("\n");
}

int main(int argc, char** argv){
	fout= NULL;

	/* Reading parameters from input */
	readParameters(argc, argv);

	/* initialize some random parameter */
	srand(param_seed);
	//srand((unsigned int)time(NULL));

	/* set logging device */
	fout = setupRecordFile();

	//load graph by jostle format
	loadCompletedGraph(param_filename);
	finalBest.best_partition = allocatePartitionData(nnode);

	for (int i = 0; i < argc; i++) {
		fprintf(fout, "%s ", argv[i]);
	}
	fprintf(fout, "\n");
	fprintf(fout, "idx\t best_v\t npat\t find_t\t find_i\t ttl_t\t  ttl_i\n");

	int run_cnt = 0;
	float sumtime = 0.0;
	int sumres = 0;
	int sumiter = 0;
	int bestInAll = -MAX_VAL;
	int *bestInAlllPartition = new int[nnode];

	/* To accelerate the test speed, we can allocate one run to each CPU node on cluster */
	while (run_cnt < 1){
		clearResult(&finalBest);
		clock_t starttime = clock();

		/* Start memetic algorithm */
		memeticRun(&finalBest, &totalgen, param_max_generations, param_time);

		totaltime = (double)(clock() - starttime)/ CLOCKS_PER_SEC;
		fprintf(fout, "%-d\t %-d\t %-d\t %-.2f\t %-d\t %-.2f\t %-d\n",run_cnt+1, finalBest.best_val, (finalBest.best_partition)->pbkt_size-1,
				finalBest.best_foundtime, finalBest.best_generation, totaltime, totalgen);

		if (finalBest.best_val > bestInAll){
			bestInAll = finalBest.best_val;
		}
		sumtime += finalBest.best_foundtime;
		sumres += finalBest.best_val;
		sumiter += finalBest.best_generation;

		reportResult();

		run_cnt++;
	}
	fprintf(fout, "best result: %d\n", bestInAll);
	fprintf(fout, "average time: %.2f\n", sumtime/run_cnt);
	fprintf(fout, "average result: %.2f\n", (float) sumres/run_cnt);
	fprintf(fout, "average best iteration: %d\n", sumiter/run_cnt);
	fclose(fout);
	delete[] bestInAlllPartition;

	return 0;
}
