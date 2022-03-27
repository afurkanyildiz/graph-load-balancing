#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include "matrix.h"

// TODO: make user defined
  // should be tuned according to DOP or some other metric(s)
  // REWRITE_UP >= REWRITE_DEPTH
//  #define REWRITE_UP      3       // upper bound to start rewriting
  #define REWRITE_DEPTH   0
  #define REWRITE_PERCENT 5


using namespace std;

typedef pair< vector<int>,vector<int> > Connectivity;
typedef vector<Connectivity> DAG;

// TODO: alignment
// TODO: delete matrices CSR & CSC
class Analyzer {
  private:
    int* levels;
    vector< vector<int> > levelTable;
    DAG dag;
    vector< vector<double>> values;
    int numOfLevels;
    bool singleLoopRows;  // if at least 2 rows exist with nnzs > 10 in a level, loops will be merged, CSR loop will be used
    // 2*nnzs-rows FLOPS in total
    double avgFLOPSPerLevel;
    vector<int> flopsPerLevel;
    map<int,double> flopsBelowAvg;
    map<int,double> flopsAboveAvg;
      // levels to avoid altering levelTable by rewriting while traversing it.
      vector<int> flopsPerLevelRewrite;

    Matrix* matrixCSR;
    Matrix* matrixCSC;


  public:
   Analyzer(Matrix* matrixCSR, Matrix* matrixCSC)
    : matrixCSR(matrixCSR), matrixCSC(matrixCSC) {}

   // TODO: any clearing to data structures needed?
    ~Analyzer() {
      matrixCSR = nullptr;
      matrixCSC = nullptr;

      if(levels != nullptr) {
        free(levels);
        levels = nullptr;
      }
    }

    int* getLevels();
    DAG& getDAG();
    int getNumOfLevels();
    bool getSingleLoopRows();
    vector< vector<int>>& getLevelTable();
    vector< vector<double>>& getValues();
    vector<int>& getFlopsPerLevel();
    map<int,double>& getFlopsBelowAvg();
    map<int,double>& getFlopsAboveAvg();
    double getAvgFLOPSPerLevel();
    void separateRows(int levelNum, int rowStartIndex, int rowEndIndex, vector<int>& loopedRows, vector<int>& unrolledRows);

      typedef map<int,pair<int,int>> ToBeRewritten;

      vector<int>& getFlopsPerLevelRewrite();
      int findMaxLevelOfPredecessors(int row);
      void correctAfterRewritingStrategy(vector<int>& emptyLevels);

    void buildLevels();
    void calculateFLOPS();
      void calculateLevelsToBeRewritten();

    void printLevels();
    void printLevelTable();
    void printDAG();
    void printValues();
    void printFLOPSPerLevel();
    void reportBefore();
      void printFLOPSPerLevelRewrite();
      void reportAfter();
      void printDependencies();
};
