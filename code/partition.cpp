#include "defines.h"
#include "hungarian.h"
const int MAX_CHAR = 10000;
const int MAX_VAL = 999999999;
const float CONST_E = 2.71828f;


void printPartition(PartitionType *ppt, int gnnode, FILE *fout){
	for (int i = 0; i < gnnode; i++){
		fprintf(fout, "%d ", ppt->pvertex[i]);
	}
	fprintf(fout, "\n");
}

/* Allocate space for partition solution */
PartitionType* allocatePartitionData(int gnnode){
	PartitionType *ppt = new PartitionType;
	ppt->pvertex = new int[gnnode];
	ppt->pcnt = new int[gnnode+1];
	ppt->ppos = new int[gnnode+1];
	ppt->pbkt = new int[gnnode+1];

	return ppt;
}

/* Build partition data structure from vpart list */
void buildPartition(PartitionType *ppt, int *vpart, int gnnode){
	for (int i = 0; i < gnnode+1; i++){
		ppt->ppos[i] = i;
		ppt->pbkt[i] = i;
	}
	ppt->pbkt_size = 1;
	memset(ppt->pcnt, 0, sizeof(int) * (gnnode+1));
	memset(ppt->pvertex, -1, sizeof(int) * gnnode);

	for (int i = 0; i < gnnode; i++){
		int pid = vpart[i];
		//ASSERTATION
		if (pid == EMPTY_IDX || pid > gnnode){
			printf("Invalid partition for node %d\n", pid);
			exit(0);
		}
		if (ppt->ppos[pid] >= ppt->pbkt_size){
			//ASSERTATION
			if (ppt->pbkt_size == gnnode+1){
				printf("The bucket is full, no new partition could be added\n");
				exit(0);
			}
			int end_pid = ppt->pbkt[ppt->pbkt_size];
			swapAry(ppt->pbkt, ppt->pbkt_size, ppt->ppos[pid]);
			swapAry(ppt->ppos, end_pid, pid);
			ppt->pbkt_size++;
			ppt->pcnt[pid]++;
		}else{
			ppt->pcnt[pid]++;
		}
	}
	memcpy(ppt->pvertex, vpart, sizeof(int)*gnnode);
}

void copyPartition(PartitionType *dest, PartitionType *source, int gnnode){
	dest->pbkt_size = source->pbkt_size;
	memcpy(dest->pbkt, source->pbkt, sizeof(int) * (gnnode+1));
	memcpy(dest->pcnt, source->pcnt, sizeof(int) * (gnnode+1));
	memcpy(dest->ppos, source->ppos, sizeof(int) * (gnnode+1));
	memcpy(dest->pvertex, source->pvertex, sizeof(int) * gnnode);
}

void disposePartition(PartitionType *ppt){
	delete[] ppt->ppos;
	delete[] ppt->pbkt;
	delete[] ppt->pcnt;
	delete[] ppt->pvertex;
	delete ppt;
}

void updatePartition(PartitionType *ppt, int v, int target){
	int oldpartition = ppt->pvertex[v];
	assert(target != EMPTY_IDX);
	if (oldpartition == target)
		printf("Meaningless move\n");
	ppt->pcnt[oldpartition]--;
	ppt->pcnt[target]++;
	//the old cluster of v is empty
	if (ppt->pcnt[oldpartition] == 0){
		int end_pid = ppt->pbkt[ppt->pbkt_size-1];
		swapAry(ppt->pbkt, ppt->pbkt_size-1, ppt->ppos[oldpartition]);
		swapAry(ppt->ppos, end_pid, oldpartition);
		ppt->pbkt_size--;
	}
	ppt->pvertex[v] = target;
}

int calculateDistance(int *p1, int *p2, int gnnode){
	int sum = 0;
	for (int i = 0; i < gnnode; i++){
		for (int j = i+1; j < gnnode; j++){
			if (p1[i] == p1[j]){
				if (p2[i] != p2[j])
					sum++;
			}else{
				if (p2[i] == p2[j])
					sum++;
			}
		}
	}
	return sum;
}

void generateRandList(int *randlist, int len){
	int idx = 0;
	assert(randlist != NULL);

	for (idx = 0; idx < len; idx++){
		randlist[idx] = idx;
	}
	for (idx = 0; idx < len; idx++){//swap
		int randid = rand() % len;
		int tmp = randlist[idx];
		randlist[idx] = randlist[randid];
		randlist[randid] = tmp;
	}
}

void swapAry(int *ary, int idx1, int idx2){
	int tmp = ary[idx1];
	ary[idx1] = ary[idx2];
	ary[idx2] = tmp;
}

int size_inter_section(vector<int> *v1, vector<int> *v2){
	sort(v1->begin(), v1->end());
	sort(v2->begin(), v2->end());

	/* Allocate a vector if the inter set need to return */
	vector<int> *rt = new vector<int>(v1->size() + v2->size());
	vector<int>::iterator it = set_intersection(v1->begin(), v1->end(), v2->begin(), v2->end(), rt->begin());
	rt->resize(it - rt->begin());
	int size = rt->size();
	delete rt;
	return size;
}

int calculateMaxMatch(int *p1, int n_part1, int *p2, int n_part2){
	vector<vector<int>* > group1(n_part1);
	vector<vector<int>* > group2(n_part2);
	int *par2idx = new int[nnode+1];
	int order_cnt = 0;
	for (int i = 0; i < nnode+1; i++)
		par2idx[i] = -1;
	for (int i = 0; i < nnode; i++){
		if(par2idx[p1[i]] == -1){
			par2idx[p1[i]] = order_cnt++;
			group1[par2idx[p1[i]]] =new vector<int>();
		}
		group1[par2idx[p1[i]]]->push_back(i);
	}
	assert(order_cnt == n_part1);
	order_cnt = 0;
	for (int i = 0; i < nnode+1; i++)
		par2idx[i] = -1;
	for (int i = 0; i < nnode; i++){
		if (par2idx[p2[i]] == -1){
			par2idx[p2[i]] = order_cnt++;
			group2[par2idx[p2[i]]] = new vector<int>();
		}
		group2[par2idx[p2[i]]]->push_back(i);
	}
	assert(order_cnt == n_part2);

	/* allocate a matrix */
	int row = max(n_part1, n_part2);
	int **cost_mat = new int*[row];
	for (int i = 0; i < row; i++){
		cost_mat[i] = new int[row];
	}

	for (int i = 0; i < row; i++){
		for (int j = 0; j < row; j++){
			if (i < n_part1 && j < n_part2){
				cost_mat[i][j] =  size_inter_section(group1[i], group2[j]);
			}else{
				cost_mat[i][j] = 0;
			}
		}
	}
	hungarian_problem_t p;
	int matrix_size = hungarian_init(&p, cost_mat, row, row, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;

	/* solve the assignment problem */
	hungarian_solve(&p);

	int sum = 0;
	for (int i = 0; i < matrix_size; i++){
		for(int j = 0; j < matrix_size; j++){
			sum += cost_mat[i][j] * p.assignment[i][j];
		}
	}

	/* free used memory */
	hungarian_free(&p);
	delete[] par2idx;
	for (int i = 0; i < row; i++){
		if (i < n_part1) delete group1[i];
		if (i < n_part2) delete group2[i];
		delete[] cost_mat[i];
	}
	return nnode-sum;
}

