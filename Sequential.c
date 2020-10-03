#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>

////////Classic merge sort changed so it can sort matrix arr2 based on matrix arr (sort indeces based on their degrees)

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(int arr[],int arr2[], int l, int m, int r)
{
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	/* create temp arrays */
	int L[n1], R[n2], L2[n1], R2[n2];
	
	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++)
		L[i] = arr[l + i];
	for (j = 0; j < n2; j++)
		R[j] = arr[m + 1 + j];

	for (i = 0; i < n1; i++)
		L2[i] = arr2[l + i];
	for (j = 0; j < n2; j++)
                R2[j] = arr2[m + 1 + j];

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
void mergeSort(int arr[],int arr2[], int l, int r)
{
	if (l < r) {
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		int m = l + (r - l) / 2;

		// Sort first and second halves 
		mergeSort(arr, arr2, l, m);
		mergeSort(arr, arr2, m + 1, r);

		merge(arr, arr2, l, m, r);
	}
}

//End of Merge Sort

////////////  Implementation of a classic queue
struct Queue {
	//Queue variables
	int front, rear, size;
	unsigned capacity;
	int* array;
};

//Function to reset queue pointers so in every iteration front,size=0
void resetQueue(struct Queue* queue) {
	queue->front = queue->size = 0;
	queue->rear = queue->capacity - 1;
}
//Queue constructor
struct Queue* createQueue(unsigned capacity)
{
	struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
	queue->capacity = capacity;
	queue->front = queue->size = 0;

	queue->rear = capacity - 1;
	queue->array = (int*)malloc(queue->capacity * sizeof(int));
	return queue;
}
//Functio to check if queue if full
int isFull(struct Queue* queue)
{
	return (queue->size == queue->capacity);
}
//Functio to check if queue is empty
int isEmpty(struct Queue* queue)
{
	return (queue->size == 0);
}

//Function to enqueue a new item
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


//Function to extract the first item of the queue
int dequeue(struct Queue* queue)
{
	if (isEmpty(queue))
		return INT_MIN;
	int item = queue->array[queue->front];
	queue->front = (queue->front + 1) % queue->capacity;
	queue->size = queue->size - 1;
	return item;
}

//Returns the index that has the minimum value
int calculateMin(int* array, int size) {
	int minVal = INT_MAX;
	int min = 0;
	for (int i = 0; i < size; i++) {
		if (array[i] < minVal) {
			minVal = array[i];
			min = i;
		}
	}
	return min;
}

//Function to calculate the degrees of the nodes, row_index is the row index matrix of the CSR sparse matrix format
int* calculateDegrees(int* row_index,int n) {
	int* degrees = malloc(sizeof(int) * n);
	for (int i = 0; i < n; i++) {
		degrees[i] = row_index[i + 1] - row_index[i];
	}
	return degrees;
}

//Plain old array printing for testing
void printArray(int A[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		printf("%d ", A[i]);
	printf("\n");
}

/*Main function where row_index is the row index matrix of the CSR sparse matrix format, col_index is the matrix that contains the column indeces of the non-zero elements of the sparse matrix in CSR format and n is the number of rows/columns of the matrix n*n*/
int* rCuthillMckee(int* row_index, int* col_index,int n) {
	//Initialization of the queue and the permutationVector(result)
	struct Queue* queue = createQueue(n);
	int* permutationVector = malloc(sizeof(int) * n);
	int vSize = 0;
	for (int i = 0; i < n; i++) {
		permutationVector[i] = n;
	}
	//Matrix N*1 that has 0 in the place i if the node i has not been visited yet, and 1 if it has
	int *visited=malloc(sizeof(int)*n);
	for(int i=0;i<n;i++){
		visited[i]=0;
	}
	//Initialize and find the degrees of all nodes of the matrix
	int *degrees = calculateDegrees(row_index,n);
	int *degOfNeighbours=malloc(sizeof(int)*n);
	//Initialize the permVector with the element with the minimum degree and put it in the visited matrix
	permutationVector[0] = calculateMin(degrees,n);
	visited[permutationVector[0]]=1;
	vSize = 1;
	
	//Iterate for every element of the permutation vector
	for (int i = 0; i < n; i++) {
		/*if(i==vSize && vSize<N){
			int *unvisited=malloc(sizeof(int)*(N-vSize));
			int *unvisitedDeg=malloc(sizeof(int)*(N-vSize));
			int size=0;
			for(int j=0;j<N;j++){
				if(!alreadyInVector(permutationVector,vSize,j)){
					unvisited[size]=j;
					unvisitedDeg[size]=degrees[j];	
				}
			}
			permutationVector[i]=unvisited[calculateMin(unvisitedDeg,N-vSize)];
			free(unvisited);
			free(unvisitedDeg);		
		}
		*/
		/*col[j] is the index of a neighbour node to the node we are examining, first we check if it has not yet been visited and if not, we enqueue it and mark it in the visited matrix. Also we put its degree on the array with the degrees of the neibhbours of the node under examination*/
		for (int j = row_index[permutationVector[i]]; j < row_index[permutationVector[i]+1]; j++) {
			if(!visited[col_index[j]]){
				enqueue(queue, col_index[j]);
				visited[col_index[j]]=1;
				degOfNeighbours[j-row_index[permutationVector[i]]]=degrees[queue->array[j-row_index[permutationVector[i]]]];
			}
			}
		//Now we sort both the degrees and our queue
		mergeSort(degOfNeighbours,queue->array, queue->front, queue->size-1);
		//And while the queue is not empty we dequeue into the permutation vector and increase its size
		while (!isEmpty(queue)) {
			permutationVector[vSize] = dequeue(queue);
			vSize++;
		}
		//We reset the queue to use it for the next iteration in order to avoid problems with indexing in mergesort
		resetQueue(queue);
	}
	//free allocated arrays
	free(queue->array);
	free(visited);
	free(degrees);
	free(degOfNeighbours);
	
	//Reverse the vector
	int temp;
	for (int i = 0; i < n/2;i++) {
		temp=permutationVector[i];
		permutationVector[i]=permutationVector[n-i-1];
		permutationVector[n-i-1]=temp;
	}
	printf("vSize= %d\n",vSize);
	return permutationVector;

}







int main() {
	
	struct timeval start,end;
	
	int n_r=1177;
	int n_c=18552;
	int n=1176;
 	
	int *rows=malloc(sizeof(int)*n_r);
	int *cols=malloc(sizeof(int)*n_c);
	FILE *ptr;
	
	ptr = fopen("rowsEris.bin","rb");
	fread(rows,sizeof(int),n_r,ptr);
	fclose(ptr);
	
	ptr=fopen("colsEris.bin","rb");
	fread(cols,sizeof(int),n_c,ptr);
	fclose(ptr);
	
	int *vector2=malloc(sizeof(int)*n);
	
	gettimeofday(&start,NULL);
	vector2=rCuthillMckee(rows,cols,n);	
	gettimeofday(&end,NULL);
	
	ptr=fopen("roadsPerm.bin","wb");
	fwrite(vector2,sizeof(int),n,ptr);
	fclose(ptr);
		
	printf("%ld\n",end.tv_usec-start.tv_usec+1000000*(end.tv_sec-start.tv_sec));
	//printArray(vector2,n);
	return 0;



}
