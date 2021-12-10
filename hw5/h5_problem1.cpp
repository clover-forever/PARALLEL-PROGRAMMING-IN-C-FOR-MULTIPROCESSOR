#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<climits>
#include <ctime>
using namespace std;

const int max = 99999; //max number in a sequence to sort

int thread_count;

void Count_sort(int a[], int n) { //count sort algorithm
	int count;
	int* temp = (int*)malloc(n*sizeof(int));

    #pragma omp parallel for private(count) num_threads(thread_count)
	for(int i = 0; i < n; i++) {
		count = 0;
		for(int j = 0; j < n; j++) {
			if (a[j] < a[i]) {
                count++; 
            }
            else if (a[j] == a[i] && j < i){
                count++;
            }
		}
		temp[count] = a[i];
	}

	memcpy(a, temp, n*sizeof(int));
	free(temp);
}

int main(int argc, char* argv[])
{
	int a[max], n; //a[max]:array, n:how many number in a input sequence

	printf("Please input how many number you want to sort: ");
	scanf("%d", &n);
	printf("Input the sequence you want to sort: ");
	for(int i = 0; i < n; i++){
        scanf("%d", a + i);
    }

	thread_count = atoi(argv[1]); //get number of threads

    double startwtime = time(NULL); //start count sort
	Count_sort(a, n);
    double endwtime = time(NULL); //end count sort

	for(int i = 0; i < n; i++){ //print out sorted array
        printf("%d ", a[i]);
    }
	printf("\n");
    printf("The execution time = %f\n", 1000*(endwtime - startwtime));
    //cout << "The execution time = "<< endwtime - startwtime <<endl ;
	return 0;
}