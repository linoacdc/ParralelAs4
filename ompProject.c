#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/time.h>


//Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(int arr[], int arr2[], int l, int m, int r)
{
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	/* create temp arrays */
	//int* L2 = malloc(sizeof(int) * n1);
	//int* R2 = malloc(sizeof(int) * n2);
	//int* L = malloc(sizeof(int) * n1);
	//int* R = malloc(sizeof(int) * n2);
	int L[n1], R[n2], L2[n1], R2[n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++){
		L[i] = arr[l + i];
		L2[i] = arr2[l + i];
	}
	for (j = 0; j < n2; j++){
		R[j] = arr[m + 1 + j];
		R2[j] = arr2[m + 1 + j];
	}
	//for (i = 0; i < n1; i++)
	//	L2[i] = arr2[l + i];
	//for (j = 0; j < n2; j++)
	//	R2[j] = arr2[m + 1 + j];

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray 
	j = 0; // Initial index of second subarray 
	k = l; // Initial index of merged subarray 
	while (i < n1 && j < n2) {
		if (L[i] <= R[j]) {
			arr[k] = L[i];
			arr2[k]=L2[i];
			i++;
		}
		else {
			arr[k] = R[j];
			arr2[k]=R2[j];
			j++;
		}
		k++;
	}

	/* Copy the remaining elements of L[], if there
	   are any */
	while (i < n1) {
		arr[k] = L[i];
		arr2[k]=L2[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there
	   are any */
	while (j < n2) {
		arr[k] = R[j];
		arr2[k]=R2[j];
		j++;
		k++;
	}
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort_serial(int arr[],int arr2[], int l, int r)
{
	if (l < r) {
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		int m = l + (r - l) / 2;

		// Sort first and second halves 
		mergeSort_serial(arr,arr2, l, m);
		mergeSort_serial(arr,arr2, m + 1, r);

		merge(arr,arr2, l, m, r);
	}
}
void mergeSort_parallel_omp(int arr[],int arr2[], int l, int r, int threads){
	
	if(threads==1||r-l<threads){
		mergeSort_serial(arr,arr2,l,r);	
	}
	else if(threads>1&&l<r){
		int m = l + (r - l) / 2;
		#pragma omp parallel sections
		{
			#pragma omp section
			mergeSort_parallel_omp(arr,arr2,l,m,threads/2);
			#pragma omp section
			mergeSort_parallel_omp(arr,arr2,m+1,r,threads-threads/2);
		}
	merge(arr,arr2,l,m,r);	
	}
}



struct Queue {
	int front, rear, size;
	unsigned capacity;
	int* array;
};

void resetQueue(struct Queue* queue) {
	queue->front = queue->size = 0;
	queue->rear = queue->capacity - 1;
}

struct Queue* createQueue(unsigned capacity)
{
	struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
	queue->capacity = capacity;
	queue->front = queue->size = 0;

	// This is important, see the enqueue 
	queue->rear = capacity - 1;
	queue->array = (int*)malloc(queue->capacity * sizeof(int));
	return queue;
}

int isFull(struct Queue* queue)
{
	return (queue->size == queue->capacity);
}

int isEmpty(struct Queue* queue)
{
	return (queue->size == 0);
}

void enqueue(struct Queue* queue, int item)
{
	if (isFull(queue)) {
		printf("Queue full!\n");
		return;
	}
	queue->rear = (queue->rear + 1) % queue->capacity;
	queue->array[queue->rear] = item;
	queue->size = queue->size + 1;
}

int dequeue(struct Queue* queue)
{
	if (isEmpty(queue))
		return INT_MIN;
	int item = queue->array[queue->front];
	queue->front = (queue->front + 1) % queue->capacity;
	queue->size = queue->size - 1;
	return item;
}

//Returns pointer to the min
int calculateMin_parallel(int* array, int size,int num_threads) {
	printf("Threads=%d\n",num_threads);
	int chunk=size/num_threads;
	printf("Chunk: %d\n",chunk);
	int minVal = size+1;
	int min = -1;
	#pragma omp parallel
	{	

		int local_min=minVal;
		int local_index=min;
		#pragma omp for nowait schedule(dynamic,chunk)
		for (int i = 0; i < size; i++) {
			if (array[i] < local_min) {
				minVal = array[i];
				min = i;
		}
		#pragma omp critical
		{
			if(local_min<minVal){
				minVal=local_min;
				min=local_index;
		}
		}
		}
	}
	return min;
}

//Returns the index that has the minimum value
	int calculateMin(int* array, int size) {
		int minVal = INT_MAX;
                int min = 0;
                for (int i = 0; i < size; i++) {
                	if (array[i] < minVal) {
                        	minVal = array[i];
                        }
		}                                                                                                         return min;
   }


int* calculateDegrees(int* row_index,int n) {
	int* degrees = malloc(sizeof(int) * n);
	for (int i = 0; i < n; i++) {
		degrees[i] = row_index[i + 1] - row_index[i];
	}
	return degrees;
}

int* calculateDegrees_parallel_omp(int* row_index,int n,int num_threads) {
	int* degrees = malloc(sizeof(int) * n);
	int chunk=n/num_threads;
	if(chunk<=1)
		return calculateDegrees(row_index,n);
	else
		#pragma omp parallel for \
			shared(degrees,chunk,row_index) \
			schedule(static,chunk)
        		for (int i = 0; i < n; i++) {
				degrees[i] = row_index[i + 1] - row_index[i];
			}
        		return degrees;
}

void printArray(int A[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		printf("%d ", A[i]);
	printf("\n");
}





int* rCuthillMckee_parallel(int* row_index, int* col_index,int n,int num_threads) {

	struct Queue* queue = createQueue(n);
	int* permutationVector = malloc(sizeof(int) * n);
	int vSize = 0;
	int chunk=n/num_threads;
	int *visited=malloc(sizeof(int)*n);
	#pragma omp parallel for \
		shared(permutationVector,chunk) \
	        schedule(dynamic,chunk)
		for (int i = 0; i < n; i++) {
			permutationVector[i] = -1;
			visited[i]=0;
		}
	int *degrees = calculateDegrees_parallel_omp(row_index,n,num_threads);
	int *degOfNeighbours=malloc(sizeof(int)*n);
	permutationVector[0] = calculateMin_parallel(degrees,n,num_threads);
	vSize = 1;
	for (int i = 0; i < n-1; i++) {
		for (int j = row_index[permutationVector[i]]; j < row_index[permutationVector[i]+1]; j++) {
			//check which elements are already in the permutation vector 
			if(!visited[col_index[j]]){
				enqueue(queue, col_index[j]);
				visited[col_index[j]]=1;
				degOfNeighbours[j-row_index[permutationVector[i]]]=degrees[queue->array[j-row_index[permutationVector[i]]]];
			}
		}
		
		//for(int k = 0; k < queue->size; k++){
		//	degOfNeighbours[k]=degrees[queue->array[k]];
		//}
		

		mergeSort_parallel_omp(degOfNeighbours, queue->array, queue->front, queue->size-1,num_threads);
		

		while (!isEmpty(queue)) {
			permutationVector[vSize] = dequeue(queue);
			vSize++;
		}
		resetQueue(queue);

		
	}
	free(queue->array);
	free(visited);
	free(degrees);
	free(degOfNeighbours);

	//REVERSE
	int temp;

	for (int i = 0; i < n/2; i++) {
		temp=permutationVector[i];
		permutationVector[i]=permutationVector[n-i-1];
		permutationVector[n-i-1]=temp;
	}
	


	
	printf("vSize=%d\n",vSize);
	return permutationVector;

}







int main() {

	struct timeval start,end;
	int num_threads;
	double startTime, endTime;
	 #pragma omp parallel
	        {
			#pragma omp master
		        {
                        	num_threads = omp_get_num_threads();
	                }
	        }
	
	int n_r=1177;
	int n_c=18552;
	int n=1176;
	
	int *rows=malloc(sizeof(int)*n_r);
	int *cols=malloc(sizeof(int)*n_c);
	FILE *ptr;

	ptr = fopen("rowsEris.bin","rb");
	fread(rows,sizeof(int),n_r,ptr);
	fclose(ptr);
	//printArray(rows,n_r);
	
	ptr=fopen("colsEris.bin","rb");
	fread(cols,sizeof(int),n_c,ptr);
	fclose(ptr);
	//printArray(cols,n_c);
	
	
	int *vector2=malloc(sizeof(int)*n);
	gettimeofday(&start,NULL);
	vector2=rCuthillMckee_parallel(rows,cols,n,num_threads);	
	gettimeofday(&end,NULL);
	
	//printArray(vector2,n);
	printf("%ld\n",end.tv_usec-start.tv_usec+1000000*(end.tv_sec-start.tv_sec));

	return 0;

}

