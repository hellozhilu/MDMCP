#include "defines.h"
#include "hungarian.h"
#include <algorithm>
using namespace std;

typedef struct ST_Pool{
	vector<PartitionType*> vec_ppt; //Each solution
	vector<int> vec_objval;			//Object value of each solution
	vector<int> vec_distance;		//Distances of each solution to the whole population
	int obj_max;					//The max objective value
	int obj_min; 					//The minimum objective value in the pool
	int dis_min;					//The minimum distance of each solution
	double obj_ave;					//The average object value
}PopulationType;


/*******************************Manage the population********************************/
void updatePopulation(PopulationType *pool){
	//Alias
	vector<PartitionType*> &pv = pool->vec_ppt;
	vector<int> &pobj = pool->vec_objval;
	vector<int> &pdis = pool->vec_distance;
	double ave_dis = 0.0;
	for (int i = 0; i < (int)pv.size(); i++){
		int dis_i2pool = MAX_VAL; //The minimum dis between two solution
		for (int j = 0; j < (int)pv.size(); j++){
			if (i == j)
				continue;
			int dis_i2j = calculateMaxMatch(pv[i]->pvertex, pv[i]->pbkt_size-1, pv[j]->pvertex, pv[j]->pbkt_size-1);
			if (dis_i2j < dis_i2pool)
				dis_i2pool = dis_i2j;
		}
		pdis[i] = dis_i2pool;
		ave_dis += dis_i2pool;
		if (pobj[i] > pool->obj_max)
			pool->obj_max = pobj[i];
		if (pobj[i] < pool->obj_min)
			pool->obj_min = pobj[i];
		if (pdis[i] < pool->dis_min)
			pool->dis_min = pdis[i];
	}
	ave_dis = ave_dis / param_pool_size;
}

int insertPopulationWhenFull2(PopulationType *pool, PartitionType *ppt, int objval, int gen){
	vector<PartitionType*> &pv = pool->vec_ppt;
	vector<int> &pobj = pool->vec_objval;
	vector<int> &pdis = pool->vec_distance;
	int replace_indv = -1;

	int minDst = nnode; //The minimum distance from current solution to the pool
	int idxMinDst = -1;

	int worst_obj = MAX_VAL;
	int worst_obj_indv = -1;
	for (int i = 0; i < (int)pv.size(); i++){
		int dist2p = calculateMaxMatch(ppt->pvertex, ppt->pbkt_size-1, pv[i]->pvertex, pv[i]->pbkt_size-1);
		if (dist2p < minDst){//minDst=is the minimum distance of the current solution to the pool
			minDst = dist2p;
			idxMinDst = i;
		}
		if (pobj[i] < worst_obj){
			worst_obj = pobj[i];
			worst_obj_indv = i;
		}
	}
	if (minDst > 0 && objval > worst_obj){
		replace_indv = worst_obj_indv;
		printf("replace %d because quality %d\n", replace_indv, worst_obj);
	}

	if (replace_indv != -1){
		copyPartition(pv[replace_indv], ppt, nnode);
		pobj[replace_indv] = objval;
		updatePopulation(pool);
	}
	return replace_indv;
}

int addPopulation(PopulationType *pool, PartitionType *ppt, int objval, int gen){
	vector<PartitionType*> &pv = pool->vec_ppt;
	vector<int> &pdis = pool->vec_distance;
	vector<int> &pobj = pool->vec_objval;
	int index = -1;
	if ((int)pv.size() < param_pool_size){/*pool initializing part*/
		int add = 1;
		if (pv.size() > 0){
			for (unsigned int i = 0; i < pv.size(); i++){
				int dis = calculateMaxMatch(pv[i]->pvertex, pv[i]->pbkt_size-1, ppt->pvertex, ppt->pbkt_size-1);
				if (dis == 0){
					printf("Solution %d is too close to the previous one %d, (dis %d)\n", objval, pobj[i], dis);
					add = 0;
					break;
				}
			}
		}
		if(add){
			index = pv.size();
			PartitionType *pptnew = allocatePartitionData(nnode);
			copyPartition(pptnew, ppt, nnode);
			pv.push_back(pptnew);
			pobj.push_back(objval);
			pdis.push_back(MAX_VAL);
			//If the pool is full, update all the values
			if ((int)pv.size() == param_pool_size)
				updatePopulation(pool);
		}
	}else { /*The pool is full, pool updating part*/
		index = insertPopulationWhenFull2(pool, ppt, objval, gen);
//#if (DEBUG_SWT)
		printf("Population: [");
		for (unsigned int i = 0; i < pv.size();i++){
			printf("%d-%d, ",pobj[i], pdis[i]);
		}
		printf("]\n");
//#endif
	}
	int obj_sum = 0;
	for (unsigned int i = 0; i < pv.size(); i++){
		obj_sum += pobj[i];
	}
	pool->obj_ave = obj_sum / (double)pv.size();
	return index;
}

PopulationType* allocatePopulation(){
	PopulationType *pool = new PopulationType;
	pool->obj_min = MAX_VAL;
	pool->obj_max = -MAX_VAL;
	pool->dis_min = MAX_VAL;
	pool->obj_ave = 0.0f;
	return pool;
}

/* release all the memories */
void disposePopulation(PopulationType *p){
	for (int i = 0; i < (int)p->vec_ppt.size();i++)
		disposePartition((PartitionType*)p->vec_ppt[i]);
	delete p;
}
/*******************************END of population manage*****************************/

/*******************************Merge-Divide Crossover*******************************/
void cliqueCover(int **we_mat, int **con_mat, int n, int *regrp, int *weight){
	set<int> candset;
	int *adjweight = new int[n];
	int grpid; //subset IDs start from id 1

	for (int i = 0; i < n; i++){//regrp[i] represents the current partially solution.
		regrp[i] = 0; //0 represents that current vertex is not assigned to a group yet.
	}
	grpid = 1;
	*weight = 0;
	while (1){
		//choose an unmarked vertex, assigned to new group
		int curvtx = -1;
		for (int i = 0; i < n; i++){
			if (!regrp[i]){
				curvtx = i;
				break;
			}
		}
		if (curvtx == -1)//all the vertices are labeled
			break;
		regrp[curvtx] = grpid;
		candset.clear();
		memset(adjweight, 0, sizeof(int) * n);
		for (int i = curvtx+1; i < n; i++){/*update current candset[]*/
			if (!regrp[i] && con_mat[curvtx][i] == -1){//vertex i is not assigned
				candset.insert(i);
				adjweight[i] = we_mat[curvtx][i];
			}
		}
		while (candset.size() > 0){
			int maxw = -MAX_VAL;
			/*Find the vertex with maximum edge weight to the current group
			 * Assign this vertex to the current group
			 * then update the adjacent edge weight*/
			int vmax = -1;
			for (set<int>::iterator itr = candset.begin(); itr != candset.end(); itr++){
				if (adjweight[(int)(*itr)] > maxw){
					maxw = adjweight[(int)(*itr)];
					vmax = (int)*itr;
				}
			}
			if (maxw < 0)//the expanding procedure for current subset aborts, then for a new subset
				break;
			candset.erase(vmax);
			regrp[vmax] = grpid;
			*weight += maxw;
//			printf("\nVtx %d grp %d, weight %d\n", vmax, grpid, *weight);
			for (int i = curvtx+1; i < n; i++){
				set<int>::iterator itr = candset.find(i);
				if (itr != candset.end()){
					if (con_mat[vmax][i] == 0){
						candset.erase(itr);
					}else{
						adjweight[i] += we_mat[vmax][i];
					}
				}
			}
		}
		grpid++;
	}
	for (int i = 0; i < n; i++){
		*weight += we_mat[i][i];
	}
	delete[] adjweight;
}

int findSet(int *s, int i){
	int r = s[i];
	while (s[r] != r){
		r = s[r];
	}
	return r;
}

int unionSet(int *s, int i, int j){
	int sameset = 0;
	int ri = findSet(s, i);
	int rj = findSet(s, j);
	if (ri < rj){
		s[rj] = ri;
	}else if (ri > rj){
		s[ri] = rj;
	}else{//ri=rj
		sameset = 1;
	}
	return sameset;
}

void MergeCrossover(PartitionType *p1, PartitionType *p2, int pval1, int pval2,
					PartitionType *child_part, int *child_val){
	int scale = nnode * param_shrink + rand() % 100; //the minimum vertices in the coarsen graph
	while (scale >= nnode)
		scale = nnode * param_shrink + rand() % 100;
	int **merge_sol = new int*[nnode];
	int *root = new int[nnode];	//the root of each vertex, the coarsen vertices are indicated in the same label (smallest)
	for (int i = 0; i < nnode; i++){
		merge_sol[i] = new int[nnode];
		memset(merge_sol[i], -1, sizeof(int) * nnode);
		merge_sol[i][i] = 1;
		root[i] = i;
	}
	int vrest = nnode;
	int fixcnt = 0; //The counter includes merge and divide process
	/*Merge and Divide edges*/
	while (vrest > scale){
		int i = rand() % nnode;
		int j = rand() % nnode;
		while (i == j)
			j = rand() % nnode;
		if (i > j)
			swap(i, j);
		if (merge_sol[i][j] == -1){//i and j is not merged now
			if (p1->pvertex[i] == p1->pvertex[j] && p2->pvertex[i] == p2->pvertex[j]){//i and j in the same set in both solution p1, p2
				merge_sol[i][j] = 1;
				int same = unionSet(root, i, j);
				if (!same)
					vrest--;
			}else if (p1->pvertex[i] != p1->pvertex[j] && p2->pvertex[i] != p2->pvertex[j]){//i and j in the different set
				merge_sol[i][j] = 0;
			}
			fixcnt++;
		}
	}
	printf("scale %d, fix rate %.4f \n", scale, (double)fixcnt * 2 / (nnode * (nnode - 1)));

	/*Reset the IDs of subsets and vertices*/
	int *flaggrp = new int[nnode];
	int *vtx_grp_id = new int[nnode]; //label which vertices are coarsen together
	memset(flaggrp, -1, sizeof(int) * (nnode));
	memset(vtx_grp_id, -1, sizeof(int) * (nnode));
	int gcnt = 0; //num of vertices of current coarsen graph
	for (int i = 0; i < nnode; i++){
		int ri = findSet(root, i);
		if (flaggrp[ri] == -1){
			flaggrp[ri] = gcnt;
			vtx_grp_id[i] = gcnt;
			gcnt++;
		}else{
			vtx_grp_id[i] = flaggrp[ri];
		}
	}

	/*Two martix to store the new coarsen graph
	 *They are symmetric vertices for the convenient of cplex_solve*/
	int **connectMtx = new int*[gcnt];
	int **weightMtx = new int*[gcnt];
	for (int i = 0; i < gcnt; i++){
		connectMtx[i] = new int[gcnt];
		weightMtx[i] = new int[gcnt];
		//or the edge are unfixed
		memset(connectMtx[i], -1, sizeof(int) * gcnt);
		memset(weightMtx[i], 0, sizeof(int) * gcnt);
		connectMtx[i][i] = 1;
	}
	for (int i = 0; i < nnode - 1; i++){
		for (int j = i+1; j < nnode; j++){//upper triangular matrix
			if (merge_sol[i][j] == 1){//in the same set
				assert(vtx_grp_id[i] == vtx_grp_id[j]);
			}else if (merge_sol[i][j] == 0){//in the different set
				assert(vtx_grp_id[i] != vtx_grp_id[j]);
				//fix the edge in coarse graph as conflict
				connectMtx[vtx_grp_id[i]][vtx_grp_id[j]] = 0;
				connectMtx[vtx_grp_id[j]][vtx_grp_id[i]] = 0;
			}
			weightMtx[vtx_grp_id[i]][vtx_grp_id[j]] += matrix[i][j];
			weightMtx[vtx_grp_id[j]][vtx_grp_id[i]] += matrix[i][j];
		}
	}
	//the weight of dialog vector is doubled
	for (int i = 0; i < gcnt; i++){
		weightMtx[i][i] = weightMtx[i][i]/ 2;
	}
	printf("coarse solution partition: %d\n", gcnt);

	/*greedy maximum edge-weight clique cover problem*/
	int *coarseSol = new int[gcnt];
	cliqueCover(weightMtx, connectMtx, gcnt, coarseSol, child_val);
 	int *pchild = new int[nnode];
	for (int i = 0; i < nnode; i++){
		pchild[i] = coarseSol[vtx_grp_id[i]];
	}

	buildPartition(child_part, pchild, nnode);

	//release memory
	delete[] pchild;
	delete[] coarseSol;
	for (int i = 0; i < gcnt; i++){
		delete[] connectMtx[i];
		delete[] weightMtx[i];
	}
	delete[] connectMtx;
	delete[] weightMtx;
	delete[] vtx_grp_id;
	delete[] flaggrp;
	delete[] root;
	for (int i = 0; i < nnode; i++){
		delete[] merge_sol[i];
	}
	delete[] merge_sol;
}
/****************************END of Merge-Divide Crossover***************************/

/**********************************MAIN memetic run**********************************/
void memeticRun(BestSolutionInfo *frt, int *total_gen, int max_gen, int max_seconds){
	double best_time = 0.0;
	clock_t start_time = clock();	//start timer
	int generation_cnt = 0;			//generaction counter

	/*Allocate memory*/
	PopulationType *population = allocatePopulation();		 //population pool
	PartitionType *child_ppt = allocatePartitionData(nnode); //(child_ppt, val) A trivial solution
	setLSEnvironment(matrix, nnode);						 //allocate the local search data

	/***************************************************/
	/* Calibrating the initial temperature */
	int child_val = generateInitSolution(matrix, nnode, child_ppt);
	setLSStart(child_ppt, child_val);
	double initTemp = calibrateTemp();

	/***************************************************/
	/* Building the initial pool */
	while((int)population->vec_ppt.size() < param_pool_size){
		int child_val = generateInitSolution(matrix, nnode, child_ppt);
		setLSStart(child_ppt, child_val);
		annealingSearch(initTemp, start_time, best_time);
		if (lsdata->fbest > frt->best_val){
			copyPartition(frt->best_partition, lsdata->ppt_best, nnode);
			frt->best_val = lsdata->fbest;
			frt->best_generation = generation_cnt;
			frt->best_foundtime = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			if (lsdata->fbest >= param_knownbest)
				goto end;
		}
		generation_cnt++;
		if ((double)(clock() - start_time)/ CLOCKS_PER_SEC >= max_seconds)
			goto end;
		addPopulation(population, lsdata->ppt_best, lsdata->fbest, generation_cnt);
	}
	printf("Best/Average value in the init population %d/%.3f\n", population->obj_max, population->obj_ave);

	/********************************************************/
	/* Memetic run: limit the up range of time or iteration */
	while (1){
		assert(param_pool_size > 0 && (int)population->vec_ppt.size() == param_pool_size);
		printf("\n------------The %dth generation-----------\n",generation_cnt);

		/****************************************************/
		/* Two parents selection and Merge-divide crossover */
		int idx1, idx2;
		idx1 = rand() % param_pool_size;
		idx2 = rand() % param_pool_size;
		while (idx1 == idx2) {
			idx2 = rand() % param_pool_size;
		}
		MergeCrossover(population->vec_ppt[idx1], population->vec_ppt[idx2],
					   population->vec_objval[idx1], population->vec_objval[idx2], child_ppt, &child_val);

		/****************************************************/
		/* Improved child solution by Simulated Annealing */
		setLSStart(child_ppt, child_val);
		annealingSearch(initTemp, start_time, best_time);
		printf("Child has been raised to by SA %d\n", lsdata->fbest);

		/****************************************************/
		/* Update the current best value after improvement */
		if (lsdata->fbest > frt->best_val){
			copyPartition(frt->best_partition, lsdata->ppt_best, nnode);
			frt->best_val = lsdata->fbest;
			frt->best_foundtime = (double)(clock() - start_time)/ CLOCKS_PER_SEC;
			frt->best_generation = generation_cnt;
			if (frt->best_val >= param_knownbest){
				goto end;
			}
		}

		/****************************************************/
		/* Pool updating */
		addPopulation(population, lsdata->ppt_best, lsdata->fbest, generation_cnt);
		printf("In %d generation, best val: %d, ave val: %.2f, min dis: %d\n",
			   generation_cnt, population->obj_max, population->obj_ave, population->dis_min);
		generation_cnt++;
		if ((double)(clock() - start_time)/ CLOCKS_PER_SEC >= max_seconds){
			goto end;
		}
	}
end:
	*total_gen = generation_cnt;
	disposeLSEnvironment();
	disposePartition(child_ppt);
	disposePopulation(population);
	return;
}
