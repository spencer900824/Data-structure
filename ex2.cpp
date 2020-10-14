#include<iostream>
#include <cstdlib>
#include <cmath>
#include<ctime>
#define K_MIN 10
#define K_MAX 29

using namespace std;

int* arr = new int[1024*1024*1024];
// sortting alogorism ----------------------------


//http://alrightchiu.github.io/SecondRound/comparison-sort-insertion-sortcha-ru-pai-xu-fa.html
void insertionSort(int *arr, int size){

    for (int i = 1; i < size; i++) {
        int key = arr[i];
        int j = i - 1;
        while (key < arr[j] && j >= 0) {
            arr[j+1] = arr[j];
            j--;
        }
        arr[j+1] = key;
    }
}
//https://lakesd6531.pixnet.net/blog/post/344227318---------------------------------
void merge(int* arr, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
	int* L = new int[n1];
	int* R = new int[n2];
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
    delete [] L;
    delete [] R;
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(int* arr, int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}
//--------------------------------
int partition(int arr[], int low, int high)
{
    int pivot = arr[low];
    int i = low - 1, j = high + 1;

    while (true) {

        // Find leftmost element greater than
        // or equal to pivot
        do {
            i++;
        } while (arr[i] < pivot);

        // Find rightmost element smaller than
        // or equal to pivot
        do {
            j--;
        } while (arr[j] > pivot);

        // If two pointers met
        if (i >= j)
            return j;

        swap(arr[i], arr[j]);
    }
}

// Generates Random Pivot, swaps pivot with
// end element and calls the partition function
// In Hoare partition the low element is selected
// as first pivot
int partition_r(int arr[], int low, int high)
{
    // Generate a random number in between
    // low .. high
    srand(time(NULL));
    int random = low + rand() % (high - low);

    // Swap A[random] with A[high]
    swap(arr[random], arr[low]);

    return partition(arr, low, high);
}


/* The main function that implements
QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(int arr[], int low, int high)
{
    if (low < high) {

        /* pi is partitioning index,
        arr[p] is now
        at right place */
        int pi = partition_r(arr, low, high);

        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
//---------------------------------*/
void countingSort( int *arr , int size , int min , int max)
{
	int* tem = new int[1024*1024*1024];
	for( int i=0 ; i<size ; i++)
	{
		if( arr[i]>=min && arr[i]<=max)
			tem[ arr[i] ]++;
		else
		{
		cout << "ERROR at index " << i << ":" << arr[i]<<endl;
			return ;
		}
	}
	
	for( int t=min ,i=0 ; t<=max && i<size ; t++)
	{
		while( tem[t] > 0)
		{
			tem[t]--;
			arr[i++] = t;
		}
		
	}
	delete [] tem;
}


//https://www.geeksforgeeks.org/cpp-program-for-heap-sort/---------------------
// To heapify a subtree rooted with node i which is
// an index in arr[]. n is size of heap
void heapify(int arr[], int n, int i)
{
    int largest = i; // Initialize largest as root
    int l = 2 * i + 1; // left = 2*i + 1
    int r = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (l < n && arr[l] > arr[largest])
        largest = l;

    // If right child is larger than largest so far
    if (r < n && arr[r] > arr[largest])
        largest = r;

    // If largest is not root
    if (largest != i) {
        swap(arr[i], arr[largest]);

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// main function to do heap sort
void heapSort(int arr[], int n)
{
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        swap(arr[0], arr[i]);

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}


void almost_sorted_array( int* arr , int size)
{
	for( int i=0 ; i<size ; i++)
			{
				arr[i] = i;
			}
			int ran_cou=100;
			int rep=0;
			while( ran_cou>0 )
			{
				srand(ran_cou+rep);
				int i = rand()%size;
				if( arr[i] == i)
				{
					arr[i] = rand()%1000+1;
					ran_cou--;
				}
				else
					rep--;
			}
}

int main()
{
	clock_t start, end, total_time;

	
	for( int k = K_MIN ,size = (int)pow( 2.0,double(K_MIN) ) ; k <= K_MAX ; k++ , size*=2)// k = 10~KMAX
	{
		total_time = 0;
		for( int runtimes = 0 ; runtimes<10 ; runtimes++)// run ten times for each k
		{
			// set the elements of ALMOST SORTED array--------------------------
			almost_sorted_array(arr,size);
			
			/*set the array-----------------------
			srand(runtimes);
			for( int i=0 ; i<size ; i++)
			{
				arr[i] = rand()%1000+1;
			}
			*/
			start = clock();
			quickSort(arr ,0,size-1);// sort
			end = clock() - start;
			
			total_time+=end;
		}
		
		cout << "Sorting an almost array of 2^" << k <<" elements with quick sort cost "
			 << double(total_time) /( 10*CLOCKS_PER_SEC )<<" sec"<<endl;
	}

    delete [] arr;
    return 0;
}

