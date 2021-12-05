#include "bmp.h"

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <string>


//定義平滑運算的次數
static const int NSmooth = 1000; // Number of smoothpProcess

static const int kProcessRoot = 0; // Process root rank

/*********************************************************/
/*變數宣告：                                             */
/*  bmpHeader    ： BMP檔的標頭                          */
/*  bmpInfo      ： BMP檔的資訊                          */
/*  **BMPSaveData： 儲存要被寫入的像素資料               */
/*  **BMPData    ： 暫時儲存要被寫入的像素資料           */
/*********************************************************/
BMPHEADER bmpHeader;         
BMPINFO bmpInfo;             
RGBTRIPLE **BMPSaveData = 0; 

/*********************************************************/
/*函數宣告：                                             */
/*  readBMP    ： 讀取圖檔，並把像素資料儲存在BMPSaveData*/
/*  saveBMP    ： 寫入圖檔，並把像素資料BMPSaveData寫入  */
/*  swap       ： 交換二個指標                           */
/*  **alloc_memory： 動態分配一個Y * X矩陣               */
/*********************************************************/
int ReadBMP(const char *fileName);       // read file
int SaveBMP(const char *fileName);       // save file
RGBTRIPLE **alloc_memory(int Y, int X);   // allocate 2D memory
void Delete2DMemory(RGBTRIPLE **memory); // delete 2D memory
void Smoothing(RGBTRIPLE &ans,
               RGBTRIPLE &a, RGBTRIPLE &b, RGBTRIPLE &c, RGBTRIPLE &d,
               RGBTRIPLE &e);            // smoothing image

int main(int argc, char *argv[])
{
  /*********************************************************/
/*變數宣告：                                             */
/*  *infileName  ： 讀取檔名                             */
/*  *outfileName ： 寫入檔名                             */
/*  startwtime   ： 記錄開始時間                         */
/*  endwtime     ： 記錄結束時間                         */
/*********************************************************/
  const char *infileName = "input.bmp";    // Input file
  const char *outfileName = "output.bmp"; // Output file

  // MPI Variable (basic)
  int local_rank, comm_size; // local_rank: my current rank,   comm_size: number of process

  int *global_bmp_size = 0, *global_bmp_displacements = 0;/*  *global_bmp_size :record each process's bmp size,  
                                                              *global_bmp_displacements:record each process's bmp displacements */

  // local bmp data
  int local_bmp_height;//local_bmp_height: Height of the local_bmp
  RGBTRIPLE **local_bmp_data = NULL;      // temporary store local bmp data
  RGBTRIPLE **local_bmp_save_data = NULL; // local bmp data

  // MPI Variable (communication temporary storage)
  RGBTRIPLE **local_bmp_upper_temp = NULL; // upper storage
  RGBTRIPLE **local_bmp_lower_temp = NULL; // lower storage

  // MPI Type
  MPI_Datatype type_mpi_rgbtriple; // MPI type: RGBTRIPLE

  double start_time, end_time; 

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

  // Create MPI Type
  MPI_Type_contiguous(3, MPI_UNSIGNED_CHAR, &type_mpi_rgbtriple);
  MPI_Type_commit(&type_mpi_rgbtriple);

  //MPI_Init(&argc,&argv);

	//記錄開始時間
	//startwtime = MPI_Wtime();

	//讀取檔案
  if (local_rank == kProcessRoot){
    if (ReadBMP(infileName)){
      std::cout << "Read file successfully!!" << std::endl;
      fflush(stdout);
    }
    else{
      std::cout << "Read file fails!!" << std::endl;
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (local_rank == kProcessRoot){
    start_time = MPI_Wtime();
  }

  // Broadcast   bmpInfo.biHeight   and   bmpInfo.biWidth   to each process
  MPI_Bcast(&bmpInfo.biHeight, 1, MPI_INT, kProcessRoot, MPI_COMM_WORLD);
  MPI_Bcast(&bmpInfo.biWidth, 1, MPI_INT, kProcessRoot, MPI_COMM_WORLD);

  //動態分配記憶體給暫存空間
  if (local_rank != kProcessRoot){
    BMPSaveData = alloc_memory(bmpInfo.biHeight, bmpInfo.biWidth);
  }

  // Calculate the 'local_bmp_height' and 'local_bmp_size 'of each process
  global_bmp_size = new int[comm_size];
  global_bmp_displacements = new int[comm_size];

  local_bmp_height = (bmpInfo.biHeight / comm_size) + ((bmpInfo.biHeight % comm_size >= local_rank + 1) ? 1 : 0);
  int local_bmp_size = bmpInfo.biWidth * local_bmp_height;

  // Gather each process
  MPI_Gather(&local_bmp_size, 1, MPI_INT, global_bmp_size, 1, MPI_INT, kProcessRoot, MPI_COMM_WORLD);

  // Calculate displacements
  if (local_rank == kProcessRoot){
    for (size_t i = 0; i < comm_size; ++i){
      global_bmp_displacements[i] = ((i != kProcessRoot) ? (global_bmp_displacements[i - 1] + global_bmp_size[i - 1]) : 0);
    }
  }

  // Allocate memory for the local bmp data and temporary storage
  local_bmp_data = alloc_memory((bmpInfo.biHeight / comm_size) + 1, bmpInfo.biWidth);
  local_bmp_save_data = alloc_memory((bmpInfo.biHeight / comm_size) + 1, bmpInfo.biWidth);

  local_bmp_upper_temp = alloc_memory(1, bmpInfo.biWidth);
  local_bmp_lower_temp = alloc_memory(1, bmpInfo.biWidth);

  // Scatterv to each process
  MPI_Scatterv(*BMPSaveData, global_bmp_size, global_bmp_displacements, type_mpi_rgbtriple, *local_bmp_save_data, local_bmp_size,
               type_mpi_rgbtriple, kProcessRoot, MPI_COMM_WORLD);

  //smooth
  for (int count = 0; count < NSmooth; ++count)
  {
    // evaluate partner
    int left_partner = ((local_rank != 0) ? (local_rank - 1) : (comm_size - 1));
    int right_partner = ((local_rank != comm_size - 1)  ? (local_rank + 1) : 0);

    // Send upper communication temporary storage
    MPI_Sendrecv(&local_bmp_save_data[local_bmp_height - 1][0], bmpInfo.biWidth, type_mpi_rgbtriple, right_partner,
                 1, &local_bmp_upper_temp[0][0], bmpInfo.biWidth, type_mpi_rgbtriple, left_partner,
                 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send lower communication temporary storage
    MPI_Sendrecv(&local_bmp_save_data[0][0], bmpInfo.biWidth, type_mpi_rgbtriple, left_partner, 1, &local_bmp_lower_temp[0][0],
                  bmpInfo.biWidth, type_mpi_rgbtriple, right_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //把像素資料與暫存指標做交換
    std::swap(local_bmp_save_data, local_bmp_data);

    //進行多次的平滑運算
    for (int i = 0; i < local_bmp_height; ++i){
      //進行平滑運算
      for (int j = 0; j < bmpInfo.biWidth; ++j){
        /*********************************************************/
				/*設定上下左右像素的位置                                 */
				/*********************************************************/
        int Top = ((i > 0) ? (i - 1) : (local_bmp_height - 1));
        int Down = ((i < local_bmp_height - 1) ? (i + 1) : 0);
        int Left = ((j > 0) ? (j - 1) : (bmpInfo.biWidth - 1));
        int Right = ((j < bmpInfo.biWidth - 1) ? (j + 1) : 0);

        /*********************************************************/
				/*與上下左右像素做平均，並四捨五入                       */
				/*********************************************************/
        if (i == 0){
          // Top pixel
          Smoothing(local_bmp_save_data[i][j], local_bmp_data[i][j], local_bmp_upper_temp[0][j],
                    local_bmp_data[Down][j], local_bmp_data[i][Left], local_bmp_data[i][Right]);
        }
        else if (i == local_bmp_height - 1){
          // Button pixel
          Smoothing(local_bmp_save_data[i][j], local_bmp_data[i][j], local_bmp_data[Top][j],
                    local_bmp_lower_temp[0][j], local_bmp_data[i][Left], local_bmp_data[i][Right]);
        }
        else{
          // Middle pixel
          Smoothing(local_bmp_save_data[i][j], local_bmp_data[i][j], local_bmp_data[Top][j],
                    local_bmp_data[Down][j], local_bmp_data[i][Left], local_bmp_data[i][Right]);
        }
      }
    }
  }

  // Gatherv to root process
  MPI_Gatherv(*local_bmp_save_data, local_bmp_size, type_mpi_rgbtriple, *BMPSaveData, global_bmp_size, global_bmp_displacements,
              type_mpi_rgbtriple, kProcessRoot, MPI_COMM_WORLD);

  // print out total time
  MPI_Barrier(MPI_COMM_WORLD);
  if (local_rank == kProcessRoot){
    end_time = MPI_Wtime();
    std::cout << "The execution time = " << (end_time - start_time) << std::endl;
  }

  //寫入檔案
  if (local_rank == kProcessRoot){
    if (SaveBMP(outfileName)){
      std::cout << "Save file successfully!!" << std::endl;
      fflush(stdout);
    }
    else{
      std::cout << "Save file fails!!" << std::endl;
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  //得到結束時間，並印出執行時間
        //endwtime = MPI_Wtime();
    	//cout << "The execution time = "<< endwtime-startwtime <<endl ;
  //free(BMPData[0]);
  //free(BMPSaveData[0]);
 	//free(BMPData);
 	//free(BMPSaveData);
  Delete2DMemory(BMPSaveData);
  Delete2DMemory(local_bmp_data);
  Delete2DMemory(local_bmp_save_data);
  Delete2DMemory(local_bmp_upper_temp);
  Delete2DMemory(local_bmp_lower_temp);
  delete[] global_bmp_size;
  delete[] global_bmp_displacements;

  MPI_Type_free(&type_mpi_rgbtriple);
  MPI_Finalize();

  return 0;
}

/*********************************************************/
/* 讀取圖檔                                              */
/*********************************************************/
int ReadBMP(const char *fileName)
{
  //建立輸入檔案物件
  std::ifstream bmpFile(fileName, std::ios::in | std::ios::binary);

  // If file cannot open
  if (!bmpFile)
  {
    std::cerr << "It can't open file!!" << std::endl;
    return 0;
  }

  //讀取BMP圖檔的標頭資料
  bmpFile.read(reinterpret_cast<char*>(&bmpHeader), sizeof(BMPHEADER));

  //判決是否為BMP圖檔
  if (bmpHeader.bfType != 0x4d42)
  {
    std::cerr << "This file is not .BMP!!" << std::endl;
    return 0;
  }

  //讀取BMP的資訊
  bmpFile.read(reinterpret_cast<char*>(&bmpInfo), sizeof(BMPINFO));

  //判斷位元深度是否為24 bits
  if (bmpInfo.biBitCount != 24)
  {
    std::cerr << "The file is not 24 bits!!" << std::endl;
    return 0;
  }

  //修正圖片的寬度為4的倍數
  while (bmpInfo.biWidth % 4 != 0)
  {
    bmpInfo.biWidth++;
  }

  //動態分配記憶體
  BMPSaveData = alloc_memory(bmpInfo.biHeight, bmpInfo.biWidth);

   //讀取像素資料
    	//for(int i = 0; i < bmpInfo.biHeight; i++)
        //	bmpFile.read( (char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE));
  bmpFile.read(reinterpret_cast<char*>(BMPSaveData[0]),
               bmpInfo.biWidth * sizeof(RGBTRIPLE) * bmpInfo.biHeight);

  //關閉檔案
  bmpFile.close();

  return 1;
}

/*********************************************************/
/* 儲存圖檔                                              */
/*********************************************************/
int SaveBMP(const char *fileName)
{
  //判決是否為BMP圖檔
  if (bmpHeader.bfType != 0x4d42)
  {
    std::cout << "This file is not .BMP!!" << std::endl;
    return 0;
  }

  //建立輸出檔案物件
  std::ofstream newFile(fileName, std::ios::out | std::ios::binary);

   //檔案無法建立
  if (!newFile)
  {
    std::cout << "The File can't create!!" << std::endl;
    return 0;
  }

  //寫入BMP圖檔的標頭資料
  newFile.write(reinterpret_cast<char*>(&bmpHeader), sizeof(BMPHEADER));

  //寫入BMP的資訊
  newFile.write(reinterpret_cast<char*>(&bmpInfo), sizeof(BMPINFO));

  //寫入像素資料
        //for( int i = 0; i < bmpInfo.biHeight; i++ )
        //        newFile.write( ( char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE) );
  newFile.write(reinterpret_cast<char*>(BMPSaveData[0]),
                bmpInfo.biWidth * sizeof(RGBTRIPLE) * bmpInfo.biHeight);

  //寫入檔案
  newFile.close();

  return 1;
}

/*********************************************************/
/* 分配記憶體：回傳為Y*X的矩陣                           */
/*********************************************************/
RGBTRIPLE **alloc_memory(int Y, int X)
{
  //建立長度為Y的指標陣列
  RGBTRIPLE **temp = new RGBTRIPLE *[Y];
  RGBTRIPLE *temp2 = new RGBTRIPLE[Y * X];
  memset(temp, 0, sizeof(RGBTRIPLE) * Y);
  memset(temp2, 0, sizeof(RGBTRIPLE) * Y * X);

  //對每個指標陣列裡的指標宣告一個長度為X的陣列
  for (int i = 0; i < Y; ++i)
  {
    temp[i] = &temp2[i * X];
  }

  return temp;
}

// delete 2D array memory.
void Delete2DMemory(RGBTRIPLE **memory)
{
  // Construct and initialize pointer
  delete[] memory[0];
  delete[] memory;
}

// Smooth
void Smoothing(RGBTRIPLE &ans, RGBTRIPLE &a, RGBTRIPLE &b,
               RGBTRIPLE &c, RGBTRIPLE &d, RGBTRIPLE &e){
  ans.rgbBlue  = static_cast<double>(a.rgbBlue  + b.rgbBlue  + c.rgbBlue +
                                     d.rgbBlue  + e.rgbBlue)  / 5 + 0.5f;
  ans.rgbGreen = static_cast<double>(a.rgbGreen + b.rgbGreen + c.rgbGreen +
                                     d.rgbGreen + e.rgbGreen) / 5 + 0.5f;
  ans.rgbRed   = static_cast<double>(a.rgbRed   + b.rgbRed   + c.rgbRed +
                                     d.rgbRed   + e.rgbRed)   / 5 + 0.5f;
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