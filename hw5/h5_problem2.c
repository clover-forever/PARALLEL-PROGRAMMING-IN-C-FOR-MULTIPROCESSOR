#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<string.h>
#include<dirent.h>
#define MAX_FILES 1000 //max number of files
#define MAX_CHAR 1000 //max number of characters per line
#define MAX_KEY_COUNT 1000 //max number of keyword in key.txt

struct node{
	char* data;
	struct node* next;
};

int key_count = 0; //index of keyword in key.txt 
int num_of_key[MAX_KEY_COUNT] = {0}; //store our answer(count of keyword in key.txt)
char keyword[MAX_KEY_COUNT][MAX_CHAR]; //store keywords in key.txt

void insert_queue(char* line,struct node** queue_head,struct node** queue_tail){ //insert node to single shared queue
	struct node* tmp_node = NULL;
	tmp_node = malloc(sizeof(struct node));
	tmp_node->data = line;
	tmp_node->next = NULL;
	//insert a node to the tail of single shared queue
#pragma omp critical
	if(*queue_tail == NULL){ //empty queue
		*queue_head = tmp_node;
		*queue_tail = tmp_node;
	}
	else{ //not empty queue
		(*queue_tail)->next = tmp_node;
		*queue_tail = tmp_node;
	}
}

struct node* del_queue(struct node** queue_head,struct node** queue_tail,int my_rank){ //get node from single shared queue
	struct node* tmp_node = NULL;
	//get a node from the front of the single shared queue
	if(*queue_head == NULL){ //empty queue
		return NULL;
	}
#pragma omp critical
	{
		if(*queue_head == *queue_tail){ //single node queue
			*queue_tail = (*queue_tail)->next;
		}
		tmp_node = *queue_head;
		*queue_head = (*queue_head)->next;
	}
	return tmp_node;
}

void tokenize(char* line,int my_rank){ //tokenize from single shared queue and count keywords which we are interested in
	int i;
	char *word;
	word = strtok(line," ,.-\n\r");
	while(word != NULL){
		for(i = 0; i < key_count ; i++){
			if(strcasecmp(word,keyword[i])==0){ //upper and lower alphabet are the same in here
#pragma omp atomic
				num_of_key[i]++;
			}
		}
		word = strtok(NULL," ,.-\n\r");
	}
}

void readfile(FILE* file, struct node** queue_head, struct node** queue_tail, int my_rank){ //read file
	char* line = malloc(MAX_CHAR*sizeof(char));
	while(fgets(line, MAX_CHAR, file) != NULL){
		printf("Thread %d : %s", my_rank, line);
		insert_queue(line, queue_head, queue_tail);
		line = malloc(MAX_CHAR*sizeof(char));
	}
	fclose(file);
}

void process(int producer_count,int consumer_count,FILE* files[],int file_count){ //process producers and consumers
	int thread_count;
	thread_count = producer_count + consumer_count;
	struct node* queue_head = NULL;
	struct node* queue_tail = NULL;
	int producer_done = 0;
#pragma omp parallel num_threads(thread_count) default(none) shared(file_count,queue_head,queue_tail,files,producer_count,consumer_count,producer_done)
	{
		int my_rank = omp_get_thread_num(),i;
		if(my_rank < producer_count){ //producers
			for(i = my_rank; i < file_count; i = i + producer_count){
				readfile(files[i], &queue_head, &queue_tail, my_rank);
			}
#pragma omp critical
			producer_done++;
		}
		else{//consumers
			for(i = 0; i < 9999999; i++); //delay for insert_queue
			struct node* tmp_node;
			while(producer_done < producer_count)
			{
				tmp_node = del_queue(&queue_head,&queue_tail,my_rank);
				if(tmp_node != NULL){
#pragma omp critical
					tokenize(tmp_node->data, my_rank);
					free(tmp_node);
				}
			}
			while(queue_head != NULL){
				tmp_node = del_queue(&queue_head, &queue_tail, my_rank);
				if(tmp_node != NULL){
#pragma omp critical
					tokenize(tmp_node->data, my_rank);
					free(tmp_node);
				}
			}
		}
	}
}

int main(int argc,char* argv[]){
	int num_of_producer; //number of producer
	int num_of_consumer; //number of consumer
	if(atoi(argv[1]) == 1){
		num_of_producer = atoi(argv[1]); //1 producer
		num_of_consumer = atoi(argv[1]); //1 consumer
	}
	else{
		if(atoi(argv[1]) > 32){
			num_of_producer = 16;
			num_of_consumer = atoi(argv[1]) -16;
		}
		else{
			num_of_producer = atoi(argv[1]) / 2; //half of threads become producer
			num_of_consumer = atoi(argv[1]) - num_of_producer; //other threads become consumer
		}
	}

	int producer_count = num_of_producer; 
	int consumer_count = num_of_consumer;
	FILE* files[MAX_FILES];
	int file_count = 0;
	double startime = 0.0, endtime = 0.0;
	
	startime = omp_get_wtime();

	//read files in a directory
	DIR *dir;
	char dirname[100];
	sprintf(dirname,"/home/H54084078/article");
	dir = opendir(dirname);
	if(dir){
		struct dirent *entry;
		while((entry = readdir(dir)) != NULL){
			if(entry->d_name[0] == '.'){
				continue;
			}
			char arti[]="article/";
			strcat(arti, entry->d_name);
			files[file_count] = fopen(arti, "r");
			if(files[file_count++] == NULL){
				perror("Open file error");
				return(-1);
			}
		}
	}
	closedir(dir);

	printf("\nOur Keyword : ");
	FILE* fp;
	char line[100], *token;
	fp = fopen("key.txt","r");
	while(fgets(line, 999, fp) != NULL){
		token = strtok(line," ,.-\n\r");
		while(token != NULL){
			strcpy(keyword[key_count++], token);
			printf("%s ",keyword[key_count-1]);
			token = strtok(NULL," ,.-\n\r");
		}
	}
	printf("\n\n");
	fclose(fp);

	//producers read line from article and insert line as a node to single shared queue, consumers get node from queue and tokenize
	process(producer_count, consumer_count, files, file_count);
	
	endtime = omp_get_wtime();

#pragma omp barrier
	int i;
	printf("-----------Keyword Count Summary------------\n");
	for(i = 0; i < key_count; i++){
		printf("%s:%d\n", keyword[i], num_of_key[i]);
	}
	printf("-----------Keyword Count Summary------------\n");
	printf("Total time = %lf\n",endtime - startime);
	return 0;
}