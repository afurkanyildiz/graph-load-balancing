#include <iostream>
#include <libufget.h>
#include <fstream>
#include <memory>
#include "matrix.h"
#include "analyzer.h"
//#include "util.h"

using namespace std;

int main(int argc, char *argv[]) {
  uf_collection_t *collection;
  int matID;

  collection = uf_collection_init();

  if (argc != 2) {
    cout << "Usage: programName matID\n";
    return 0;
  }

  matID = std::stoi(std::string(argv[1]));
  printf("UF Sparse Matrix Collection Database contains %d matrices.\n", uf_collection_num_matrices(collection));

  Matrix* matrixCSR = new Matrix(collection, matID);
  int successCSR = matrixCSR->convertToCSR();
  printf("success (0) for COO --> CSR:%d\n", successCSR);

  Matrix* matrixCSC = new Matrix(collection, matID);
  int successCSC = matrixCSC->convertToCSC();
  printf("success (0) for COO --> CSC:%d\n", successCSC);

  matrixCSR->extractLCSR();
//  cout << "################################################################ Printing after extractLCSR\n";
//  matrixCSR->print();

//  matrixCSC->print();
  matrixCSC->extractLCSC();
// cout << " ############################################################### Printing after extractLCSC\n";
//  matrixCSC->print();

  Analyzer* analyzer = new Analyzer(matrixCSR, matrixCSC);
   analyzer->buildLevels();
   analyzer->printLevelTable();
//   analyzer->printDAG();
//   analyzer->printValues();
//   analyzer->printLevels();
//   analyzer->printLevelTable();
//   analyzer->calculateFLOPS();
//   analyzer->reportBefore();
//   analyzer->printDependencies();

  return 0;
}
