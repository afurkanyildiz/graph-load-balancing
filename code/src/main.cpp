#include <iostream>
#include <libufget.h>
#include <fstream>
#include <memory>
#include "matrix.h"
#include "rewrite.h"
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

  matrixCSC->extractLCSC();

  Analyzer* analyzer = new Analyzer(matrixCSR, matrixCSC);
 // analyzer->buildRowHist();
  analyzer->buildLevels();
//  analyzer->printLevelTable();
//  analyzer->printDAG();
//  analyzer->printValues();
//  analyzer->printLevels();
  analyzer->calculateFLOPS();
  analyzer->report(string ("BEFORE"));

  #ifdef REWRITE_ENABLED
    //analyzer->printDependencies();

    RewriteByCostMap* rewritingStrategy = new RewriteByCostMap(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, analyzer);

    //RewriteByThreeCriteria* rewritingStrategy = new RewriteByThreeCriteria(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, analyzer);

    rewritingStrategy->shiftRowsUp();
    vector<int> emptyLevels;
    rewritingStrategy->findEmptyLevels(emptyLevels);
    analyzer->correctAfterRewritingStrategy(emptyLevels);

    analyzer->report(string ("AFTER"));
    rewritingStrategy->printRowsToBeRewritten();

  #endif

  #ifdef REWRITE_ENABLED
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), rewritingStrategy, analyzer);
  #else
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), analyzer);
  #endif

//  rewriter.rewrite();
   
  #ifdef REWRITE_ENABLED
//    rewritingStrategy->printRowsToBeRewritten();
//    rewritingStrategy->printRewritingMap();
//    analyzer->report(string ("AFTER"));
  #endif

  return 0;
}
