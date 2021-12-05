#include <pthread.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "bmp.h"
#include <ctime>
#include <sstream>

using namespace std;

//定義平滑運算的次數
#define NSmooth 1000

/*********************************************************/
/*變數宣告：                                              */
/*  bmpHeader    ： BMP檔的標頭                           */
/*  bmpInfo      ： BMP檔的資訊                           */
/*  **BMPSaveData： 儲存要被寫入的像素資料                  */
/*  **BMPData    ： 暫時儲存要被寫入的像素資料               */
/*********************************************************/
BMPHEADER bmpHeader;                        
BMPINFO bmpInfo;
RGBTRIPLE **BMPSaveData = NULL;                                               
RGBTRIPLE **BMPData = NULL;                                                   

//thread variable
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

int thread_num; //number of thread
int local_height;       //local height for a thread to process
int counter1, counter2; //for busy waiting

/*********************************************************/
/*函數宣告：                                             */
/*  readBMP    ： 讀取圖檔，並把像素資料儲存在BMPSaveData*/
/*  saveBMP    ： 寫入圖檔，並把像素資料BMPSaveData寫入  */
/*  swap       ： 交換二個指標                           */
/*  **alloc_memory： 動態分配一個Y * X矩陣               */
/*********************************************************/
int readBMP( char *fileName);        //read file
int saveBMP( char *fileName);        //save file
void swap(RGBTRIPLE *a, RGBTRIPLE *b);
RGBTRIPLE **alloc_memory( int Y, int X );        //allocate memory
void *run(void *rank);

int main(int argc,char *argv[])
{
        /*********************************************************/
        /*變數宣告：                                             */
        /*  *infileName  ： 讀取檔名                             */
        /*  *outfileName ： 寫入檔名                             */
        /*  startwtime   ： 記錄開始時間                         */
        /*  endwtime     ： 記錄結束時間                         */
        /*********************************************************/
	char *infileName = "input.bmp";
        char *outfileName = "output_h4_problem1.bmp";
	double startwtime = 0.0, endwtime=0.0;

	//MPI_Init(&argc,&argv);
	
	//記錄開始時間
	startwtime = time(NULL);

	//讀取檔案
        if ( readBMP( infileName) )
                cout << "Read file successfully!!" << endl;
        else 
                cout << "Read file fails!!" << endl;

	//動態分配記憶體給暫存空間
        BMPData = alloc_memory( bmpInfo.biHeight, bmpInfo.biWidth);

        /*create thread*/

        stringstream ss;  
        ss << argv[1];  
        ss >> thread_num;  


        pthread_t thread_arr[thread_num];
        local_height = bmpInfo.biHeight / thread_num;
        counter1 = 0; counter2 = 0;

        pthread_mutex_init(&mutex, 0);
        //pthread_create
        for (int i = 0; i < thread_num; i++){
                pthread_create(&thread_arr[i], NULL, run, (void*)i);
        }

        //pthread_join
        for (int i = 0; i < thread_num; i++){
                pthread_join(thread_arr[i], NULL);
        }
 
 	//寫入檔案
        if ( saveBMP( outfileName ) )
                cout << "Save file successfully!!" << endl;
        else
                cout << "Save file fails!!" << endl;
	
	//得到結束時間，並印出執行時間
        endwtime = time(NULL);
    	cout << "The execution time = "<< endwtime - startwtime <<endl ;

	free(BMPData[0]);
 	free(BMPSaveData[0]);
 	free(BMPData);
 	free(BMPSaveData);
 	//MPI_Finalize();

    return 0;
}

/*********************************************************/
/* 讀取圖檔                                              */
/*********************************************************/
int readBMP(char *fileName)
{
	//建立輸入檔案物件
        ifstream bmpFile( fileName, ios::in | ios::binary );
 
        //檔案無法開啟
        if ( !bmpFile ){
                cout << "It can't open file!!" << endl;
                return 0;
        }
 
        //讀取BMP圖檔的標頭資料
    	bmpFile.read( ( char* ) &bmpHeader, sizeof( BMPHEADER ) );
 
        //判決是否為BMP圖檔
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }
 
       //讀取BMP的資訊
        bmpFile.read( ( char* ) &bmpInfo, sizeof( BMPINFO ) );
        
        //判斷位元深度是否為24 bits
        if ( bmpInfo.biBitCount != 24 ){
                cout << "The file is not 24 bits!!" << endl;
                return 0;
        }

        //修正圖片的寬度為4的倍數
        while( bmpInfo.biWidth % 4 != 0 )
        	bmpInfo.biWidth++;

        //動態分配記憶體
        BMPSaveData = alloc_memory( bmpInfo.biHeight, bmpInfo.biWidth);
        
        //讀取像素資料
    	//for(int i = 0; i < bmpInfo.biHeight; i++)
        //	bmpFile.read( (char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE));
	    bmpFile.read( (char* )BMPSaveData[0], bmpInfo.biWidth*sizeof(RGBTRIPLE)*bmpInfo.biHeight);
	
        //關閉檔案
        bmpFile.close();
 
        return 1;
 
}
/*********************************************************/
/* 儲存圖檔                                              */
/*********************************************************/
int saveBMP( char *fileName)
{
 	//判決是否為BMP圖檔
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }
        
 	//建立輸出檔案物件
        ofstream newFile( fileName,  ios:: out | ios::binary );
 
        //檔案無法建立
        if ( !newFile ){
                cout << "The File can't create!!" << endl;
                return 0;
        }
 	
        //寫入BMP圖檔的標頭資料
        newFile.write( ( char* )&bmpHeader, sizeof( BMPHEADER ) );

	//寫入BMP的資訊
        newFile.write( ( char* )&bmpInfo, sizeof( BMPINFO ) );

        //寫入像素資料
        //for( int i = 0; i < bmpInfo.biHeight; i++ )
        //        newFile.write( ( char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE) );
        newFile.write( ( char* )BMPSaveData[0], bmpInfo.biWidth*sizeof(RGBTRIPLE)*bmpInfo.biHeight );

        //寫入檔案
        newFile.close();
 
        return 1;
 
}


/*********************************************************/
/* 分配記憶體：回傳為Y*X的矩陣                           */
/*********************************************************/
RGBTRIPLE **alloc_memory(int Y, int X )
{        
	//建立長度為Y的指標陣列
        RGBTRIPLE **temp = new RGBTRIPLE *[ Y ];
	    RGBTRIPLE *temp2 = new RGBTRIPLE [ Y * X ];
        memset( temp, 0, sizeof( RGBTRIPLE ) * Y);
        memset( temp2, 0, sizeof( RGBTRIPLE ) * Y * X );

	//對每個指標陣列裡的指標宣告一個長度為X的陣列
        for( int i = 0; i < Y; i++){
                temp[ i ] = &temp2[i*X];
        }
 
        return temp;
 
}
/*********************************************************/
/* 交換二個指標                                          */
/*********************************************************/
void swap(RGBTRIPLE *a, RGBTRIPLE *b)
{
	RGBTRIPLE *temp;
	temp = a;
	a = b;
	b = temp;
}

//function for thread to process image smoothing
void *run(void *rank){
        long my_rank = (long)rank;
        for(int i = 0 ; i < NSmooth; i++){
                //busy waiting and a mutex
                pthread_mutex_lock(&mutex);
                counter1++;
                pthread_mutex_unlock(&mutex);
                while (counter1 < thread_num * (i + 1)) { 
                        ; 
                }
                
                //thread 0 doing swap
                //把像素資料與暫存指標做交換
                if (my_rank == 0){
                        swap(BMPSaveData, BMPData);
                }

                //other thread doing smoothing
                //busy waiting and a mutex
                //barrier
                pthread_mutex_lock(&mutex);
                counter2++;
                pthread_mutex_unlock(&mutex);
                while (counter2 < thread_num * (i + 1)) {
                        ;
                }

                //進行平滑運算

                for (int i = my_rank * local_height; i < (my_rank + 1) * local_height; i++){
                        for (int j = 0; j < bmpInfo.biWidth; j++)
                        {
                                /*********************************************************/
                                /*設定上下左右像素的位置                                 */
                                /*********************************************************/
                                int Top = i > 0 ? i - 1 : bmpInfo.biHeight - 1;
                                int Down = i < bmpInfo.biHeight - 1 ? i + 1 : 0;
                                int Left = j > 0 ? j - 1 : bmpInfo.biWidth - 1;
                                int Right = j < bmpInfo.biWidth - 1 ? j + 1 : 0;
                                /*********************************************************/
                                /*與上下左右像素做平均，並四捨五入                       */
                                /*********************************************************/
                                BMPSaveData[i][j].rgbBlue = (double)(BMPData[i][j].rgbBlue + BMPData[Top][j].rgbBlue + BMPData[Down][j].rgbBlue + BMPData[i][Left].rgbBlue + BMPData[i][Right].rgbBlue) / 5 + 0.5;
                                BMPSaveData[i][j].rgbGreen = (double)(BMPData[i][j].rgbGreen + BMPData[Top][j].rgbGreen + BMPData[Down][j].rgbGreen + BMPData[i][Left].rgbGreen + BMPData[i][Right].rgbGreen) / 5 + 0.5;
                                BMPSaveData[i][j].rgbRed = (double)(BMPData[i][j].rgbRed + BMPData[Top][j].rgbRed + BMPData[Down][j].rgbRed + BMPData[i][Left].rgbRed + BMPData[i][Right].rgbRed) / 5 + 0.5;
                        }
                }
        }
        return NULL;       
}
