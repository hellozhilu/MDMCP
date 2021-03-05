#include "defines.h"

SA_RT_Data *lsdata;


int generateInitSolution(int **gmatrix, int gnnode, PartitionType *ppt){
	int *initpart = new int[gnnode];
	int sum = 0;
	assert(ppt != NULL);
	for (int i = 0; i < gnnode; i++){
		initpart[i] = i + 1;  //=(1,2,...,gnnode)
		sum += gmatrix[i][i]; //=0
	}
	buildPartition(ppt, initpart, gnnode);
	delete[] initpart;
	return sum;
}

void buildCurGamma(int* pvertex){
	int lsnnode = lsdata->gnnode;
	for (int i = 0; i < lsnnode; i++){
		memset(lsdata->gammatbl[i], 0, sizeof(int) * (lsnnode+1));
		for (int j = 0; j < lsnnode; j++){
			lsdata->gammatbl[i][pvertex[j]] += lsdata->gmatrix[i][j];
		}
	}
}

void updateCurGamma(int u, int src, int dest){
	assert(src != EMPTY_IDX && dest != EMPTY_IDX);
	for (int i = 0; i < lsdata->gnnode; i++){
		lsdata->gammatbl[i][src] -= lsdata->gmatrix[u][i];
		lsdata->gammatbl[i][dest] += lsdata->gmatrix[u][i];
	}
}

/* allocate and setup all the environment for local search */
void setLSEnvironment(int **graph_mat, int gnnode){
	/* we just repoint the graph pointer,
	 * thus the graph_mat will be kept even we dispose the runtime data */
	if (lsdata == NULL)
		lsdata = new SA_RT_Data;
	assert(graph_mat != NULL);
	lsdata->gmatrix = graph_mat;
	lsdata->gnnode = gnnode;

	/* set current solution, assume the given obj value is compatible
	 *  also set the known best value for the stop condition */
	lsdata->ppt = allocatePartitionData(gnnode);
	lsdata->ppt_best = allocatePartitionData(gnnode);

	/* build pgamma data */
	lsdata->gammatbl = new int*[gnnode];
	for (int i = 0; i < gnnode; i++){
		lsdata->gammatbl[i] = new int[gnnode + 1];
	}
}

/* Reset a new start solution for local search */
void setLSStart(PartitionType *start_sol, int start_obj_value){
	assert(lsdata != NULL);

	/* reset current partition and current solution
	 * The fbest is kept */
	copyPartition(lsdata->ppt, start_sol, lsdata->gnnode);
	lsdata->fcurrent = start_obj_value;

	copyPartition(lsdata->ppt_best, start_sol, lsdata->gnnode);
	lsdata->fbest = start_obj_value;

	/* rebuild gamma table */
	buildCurGamma(start_sol->pvertex);
}

void disposeLSEnvironment(){
	disposePartition(lsdata->ppt);
	disposePartition(lsdata->ppt_best);

	for (int i = 0; i < lsdata->gnnode; i++)
		delete[] lsdata->gammatbl[i];
	delete[] lsdata->gammatbl;

	delete lsdata;
	lsdata = NULL;
}

static int decideTarget(int dest){
	PartitionType *cur_ppt = lsdata->ppt;
	int gnnode = lsdata->gnnode;
	/* Move v to a new cluster */
	if (cur_ppt->pbkt_size > gnnode+1)
		printf("ERROR\n");
	if (dest == EMPTY_IDX){
		assert(cur_ppt->pbkt_size > 0 && cur_ppt->pbkt_size <= gnnode+1);
		dest = cur_ppt->pbkt[cur_ppt->pbkt_size];
		cur_ppt->pbkt_size++;
	}
	return dest;
}

/* update partition */
static int changeCurSolution(int u, int dest){
	PartitionType *cur_ppt = lsdata->ppt;
	int target = decideTarget(dest);
	int src = cur_ppt->pvertex[u];
	if (src == target)
		fprintf(stderr, "warning: iter %d invalid move, move %d from %d to %d\n",
				lsdata->itr, u, src, dest);
	//update gamma table
	updateCurGamma(u, src, target);
	updatePartition(cur_ppt, u, target);
	return target;
}

void pureDescent(){
	int *randlst = new int[lsdata->gnnode];
	int improved = 1;

	PartitionType *ppt = lsdata->ppt;
	int **pgamma = lsdata->gammatbl;
	int nnode = lsdata->gnnode;

	while (improved){
		improved = 0;
		generateRandList(randlst, lsdata->gnnode);
		for (int i = 0; i < lsdata->gnnode; i++){
			int currentnode = randlst[i];
			for (int k = 0; k < ppt->pbkt_size; k++){
				int curpart = ppt->pbkt[k];
				int gain = pgamma[currentnode][curpart] - pgamma[currentnode][ppt->pvertex[currentnode]];;
//				printf("i=%d, k=%d, gain=%d\n", i, k, gain);
				if (gain > 0){
					improved = 1;
					changeCurSolution(currentnode, curpart);
					lsdata->fcurrent += gain;
//					cout << "f=" << lsdata->fcurrent << endl;
				}
			}
		}
	}
	if (lsdata->fcurrent > lsdata->fbest){
		copyPartition(lsdata->ppt_best, lsdata->ppt, nnode);
		lsdata->fbest = lsdata->fcurrent;
	}
	delete[] randlst;

}

/* Determine the initial temp by binary search in SA */
double calibrateTemp(){
	double lt = 1.0, ut = 2000.0;
	double tempTolerate = 0.05;
	double temp = (lt + ut) /2;
	int L;
	int find = 0;
	printf("Calibrate the initial temperature.\n");

	PartitionType *ppt = lsdata->ppt;
	int **gammatbl = lsdata->gammatbl;

	while (!find){
		int iter = 0;
		int accpCnt =  0;
		pureDescent();
		L = nnode * ppt->pbkt_size * param_sizefactor;
		while (1){
			int i = rand() % nnode;
			int bestpid = -1;
			int bestdelta = -MAX_VAL;
			for (int k = 0; k < ppt->pbkt_size; k++){
				int pid = ppt->pbkt[k];
				if (pid == ppt->pvertex[i] ||//source == dest
					(ppt->pcnt[ppt->pvertex[i]] == 1 && pid == EMPTY_IDX))//cannot change the fobj
					continue;
				int delta = gammatbl[i][pid] - gammatbl[i][ppt->pvertex[i]];
				if (delta > bestdelta){
					bestdelta = delta;
					bestpid = pid;
				}
			}
			if (bestdelta > 0){
				//accept the solution
				changeCurSolution(i, bestpid);
				lsdata->fcurrent += bestdelta;
				accpCnt++;
			}else{
				double prob = exp((double)bestdelta / temp);
				if (rand() % 1000 < prob * 1000){ //accept the solution
					changeCurSolution(i, bestpid);
					lsdata->fcurrent += bestdelta;
					accpCnt++;
				}
			}

			iter++;
			if (iter > L){
				goto evaluate;
			}
		}
evaluate:/* a binary search */
		double accpProb = (double)accpCnt /L;
		if (accpProb > 0.5 + tempTolerate){
			ut = temp;
			temp = (lt + ut) /2;
//			printf("Accept %.3f, Too hot, turn down to %.2f\n", accpProb, temp);
		}else if (accpProb < 0.5 - tempTolerate){
			lt = temp;
			temp = (lt + ut) /2;
//			printf("Accept %.3f, Too cold, turn up to %.2f\n", accpProb, temp);
		}else{
			printf("Accept probability %.3f, Calibrate temp %.2f\n\n", accpProb, temp);
			break;
		}
	}
	return temp;
}

void annealingSearch(double initTemp, clock_t start_time, double &best_time){
	PartitionType *ppt = lsdata->ppt;
	int **gammatbl = lsdata->gammatbl;
	int L;
	double T = initTemp;
	int frozenCounter = 0;
	int accpCnt = 0;
	int saitr = 0;
	int *vecTmpBest = new int[lsdata->gnnode];
	pureDescent();

	memcpy(vecTmpBest, ppt->pvertex, sizeof(int) * lsdata->gnnode);
	L = lsdata->gnnode * ppt->pbkt_size * param_sizefactor;
	while (1){
		int i = rand() % lsdata->gnnode;
		/*Find the best partition for node i*/
		int bestpid = -1;
		int bestdelta = -MAX_VAL;
		for (int k = 0; k < ppt->pbkt_size; k++){
			int pid = ppt->pbkt[k];
			if (pid == ppt->pvertex[i] || (ppt->pcnt[ppt->pvertex[i]] == 1 && pid == EMPTY_IDX))
				continue;
			int delta = gammatbl[i][pid] - gammatbl[i][ppt->pvertex[i]];
			if (delta > bestdelta){
				bestdelta = delta;
				bestpid = pid;
			}
		}

		if (bestdelta >= 0){
			//accept the solution
			changeCurSolution(i, bestpid);
			lsdata->fcurrent += bestdelta;
			accpCnt++;
		}else{//delta < 0
			double prob = exp((double)bestdelta /T);
			if (rand() % 1000 < prob * 1000){
				//accept
				changeCurSolution(i, bestpid);
				lsdata->fcurrent += bestdelta;
				accpCnt++;
			}
		}
		//record local
		lsdata->itr++;
		if (lsdata->fcurrent > lsdata->fbest){
			memcpy(vecTmpBest, lsdata->ppt->pvertex, sizeof(int) * lsdata->gnnode);
			lsdata->fbest = lsdata->fcurrent;
		}

		saitr++;
//		printf("%d-%d\n", saitr, lsdata->fcurrent);
//		printf("%d-%d-fbest=%d\n", frozenCounter, saitr, lsdata->fbest);
		if (saitr > L){
			saitr = 0;
			T = T * param_tempfactor; //cool down
			if (accpCnt * 100 < (int)param_minpercent * L){//=(double)accpCnt/ L
				frozenCounter++;
			}
			if (frozenCounter > 5){
				goto anneal_stop;
			}
			accpCnt = 0;
			saitr = 0;
		}
	}
anneal_stop:
	buildPartition(lsdata->ppt_best, vecTmpBest, lsdata->gnnode);
	delete[] vecTmpBest;
	return;
}
