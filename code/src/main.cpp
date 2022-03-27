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

//  matrixCSC->print();
  matrixCSC->extractLCSC();
  // cout << " ############################################################### Printing after extractLCSC\n";
//  matrixCSC->print();

  Analyzer* analyzer = new Analyzer(matrixCSR, matrixCSC);
 // analyzer->buildRowHist();
  analyzer->buildLevels();
//  analyzer->printLevelTable();
//  analyzer->printDAG();
//  analyzer->printValues();
  //analyzer->printRowHist();
//  analyzer->printLevels();
//  cout << "reporting before\n";
//  analyzer->printLevelTable();
  analyzer->calculateFLOPS();
//  analyzer->reportBefore();

  #ifdef REWRITE_ENABLED
    //analyzer->printDependencies();

  cout << "rewrite enabled\n";
    // MUST be called before the rewriting strategy
    analyzer->calculateLevelsToBeRewritten();
    RewriteByLevel* rewritingStrategy = new RewriteByLevel(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, 
                                        analyzer->getLevels(), analyzer->getLevelTable(), analyzer->getDAG(), analyzer->getFlopsBelowAvg());
    //RewriteByCostMap* rewritingStrategy = new RewriteByCostMap(matrixCSR->getL()->getRows() - 1, StartPoint::BottomUp, 
    //                                    analyzer->getLevels(), analyzer->getLevelTable(), analyzer->getDAG(), analyzer->getFlopsPerLevel(),
    //                                    analyzer->getAvgFLOPSPerLevel(), analyzer->getFlopsBelowAvg());

    rewritingStrategy->shiftRowsUp();
//    analyzer->printLevelTable();
    vector<int> emptyLevels;
    rewritingStrategy->findEmptyLevels(emptyLevels);
    analyzer->correctAfterRewritingStrategy(emptyLevels);

    // TODO: we dont need to calculate anymore
//    analyzer->calculateFLOPS();
    analyzer->reportBefore();
//    analyzer->printLevelTable();
    // no need for RewriteByLevel    
  #endif

  #ifdef REWRITE_ENABLED
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), rewritingStrategy, analyzer);
  #else
    Rewrite rewriter(matrixCSR, matrixCSC, std::string(matrixCSR->getUF_matrix()->name), analyzer);
  #endif

  rewriter.rewrite();
   
  #ifdef REWRITE_ENABLED
//    rewritingStrategy->printRowsToBeRewritten();
//    rewritingStrategy->printRewritingMap();
    analyzer->reportAfter();
  #endif

  return 0;
}
