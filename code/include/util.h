#pragma once
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "matrix.h"


//#define FILEPATH "/tmp/mmapped.bin"
//#define NUMINTS  (10)
//#define FILESIZE (NUMINTS * NUMINTS * sizeof(int))

// dump some data to memory to be used by the generated code
template <class T>
int dumpToMem(vector<T>& data, const char* filePath, int numOfElements) {
  int fd;
  int result;
  size_t size = numOfElements * sizeof(T);

  cout << "data size: " << data.size() << "\n";

  fd = open(filePath, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  if(fd == -1) {
    perror("Error opening file for writing");
    exit(EXIT_FAILURE);
  }

  result = lseek(fd, size-1, SEEK_SET);
  if(result == -1) {
    close(fd);
    perror("Error calling lseek() to 'stretch' the file");
    exit(EXIT_FAILURE);
  }

  result = write(fd, "", 1);
  if(result != 1) {
    close(fd);
    perror("Error writing last byte of the file");
    exit(EXIT_FAILURE);
  }

  T* map = static_cast<T*>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(map == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }

  // TODO: vectorize
  for(int i = 0; i < numOfElements; i++) 
    map[i] = data[i];

  fd = open(filePath, O_RDWR);
  if(fd == -1) {
    perror("Error opening file for reading");
    exit(EXIT_FAILURE);
  }

  T* map_r = static_cast<T*>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(map_r == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }

  if(munmap(map, size) == -1 || munmap(map_r, size) == -1)
    perror("Error un-mmapping the file");

  close(fd);
    return 0;
}

void flattenDumpToMem(Part* LCSR, const char* filePath, const char* filePathRows) {
  int fd;
  int result;

  vector<int>& rowPtrLCSR = LCSR->getRowPtr();
  vector<int>& colIdxLCSR = LCSR->getColIdx();
  int numVals = LCSR->getNNZs();
  int numRows = LCSR->getRows();

  cout << "numRows: " << numRows << " numVals: " << numVals << "\n";
  cout << "rowPtr size: " << rowPtrLCSR.size() << " colIdx size: " << colIdxLCSR.size() << "\n";

  size_t size = (numVals) * sizeof(int);

  // colIdxLCSR
  fd = open(filePath, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  if(fd == -1) {
    perror("Error opening file for writing");
    exit(EXIT_FAILURE);
  }

  result = lseek(fd, size-1, SEEK_SET);
  if(result == -1) {
    close(fd);
    perror("Error calling lseek() to 'stretch' the file");
    exit(EXIT_FAILURE);
  }

  result = write(fd, "", 1);
  if(result != 1) {
    close(fd);
    perror("Error writing last byte of the file");
    exit(EXIT_FAILURE);
  }

  int* flatData = static_cast<int (*)>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(flatData == MAP_FAILED) {
    close(fd);
    perror("Error mflatDataping the file");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < numVals ; i++)
    flatData[i] = colIdxLCSR[i]; 

  // read back
  fd = open(filePath, O_RDWR);
  if (fd == -1) {
    perror("Error opening file for reading");
    exit(EXIT_FAILURE);
  }

  int* flatData_r = static_cast<int (*)>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(flatData_r == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }

  // rowPtrLCSR
  size = (numRows) * sizeof(int);

  fd = open(filePathRows, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  if(fd == -1) {
    perror("Error opening file for writing");
    exit(EXIT_FAILURE);
  }

  result = lseek(fd, size-1, SEEK_SET);
  if(result == -1) {
    close(fd);
    perror("Error calling lseek() to 'stretch' the file");
    exit(EXIT_FAILURE);
  }

  result = write(fd, "", 1);
  if(result != 1) {
    close(fd);
    perror("Error writing last byte of the file");
    exit(EXIT_FAILURE);
  }

  int* flatDataRows = static_cast<int (*)>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(flatDataRows == MAP_FAILED) {
    close(fd);
    perror("Error mflatDataping the file");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < numRows ; i++)
    flatDataRows[i] = rowPtrLCSR[i]; 

  // read back
  fd = open(filePath, O_RDWR);
  if (fd == -1) {
    perror("Error opening file for reading");
    exit(EXIT_FAILURE);
  }

  int* flatData_rr = static_cast<int (*)>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(flatData_rr == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }


  if(munmap(flatData, size) == -1 || munmap(flatData_r, size) == -1 || munmap(flatDataRows, size) == -1 || munmap(flatData_rr, size) == -1)
    perror("Error un-mmapping the file");

  close(fd);
} 

// This is a temporary function to read the output of calculateWorkload.py
void readWorkloadPartition(vector<vector<int>>& partitionList) {
  string fName("partition_932.csv");
  ifstream file(fName, ifstream::in);
  string line;

  // Iterate through each line and split the content using delimeter
  while(getline(file, line)){
    vector<int> vec;
    stringstream values(line);
    string val;
    while(getline(values,val,',')) {
      size_t found=val.find("\r");
      if (found!=string::npos)
        val.erase(found,1);

      if(val.empty())
        break;

      int value = 0;
      stringstream ss(val);
      ss >> value;
      vec.push_back(value);
    }
    partitionList.push_back(vec);
  }

  file.close();
}


