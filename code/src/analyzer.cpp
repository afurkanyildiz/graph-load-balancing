#include <limits.h>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <typeinfo>
#include "analyzer.h"

int* Analyzer::getLevels() {
  return levels;
}

int Analyzer::getNumOfLevels() {
  return numOfLevels;
}

vector<vector<int>>& Analyzer::getLevelTable() {
  vector<vector<int>>& ref = levelTable;

  return ref;
}

DAG& Analyzer::getDAG() {
  DAG& ref = dag;

  return ref;
}

vector<vector<double>>& Analyzer::getValues() {
  vector<vector<double>>& ref = values;

  return ref;
}

vector<int>& Analyzer::getFlopsPerLevel() {
  vector<int>& ref = flopsPerLevel;

  return ref;
}

map<int,double>& Analyzer::getFlopsBelowAvg() {
  map<int,double>& ref = flopsBelowAvg;

  return ref;
}

map<int,double>& Analyzer::getFlopsAboveAvg() {
  map<int,double>& ref = flopsAboveAvg;

  return ref;
}

double Analyzer::getAvgFLOPSPerLevel() {
  return avgFLOPSPerLevel;
}

bool Analyzer::getSingleLoopRows() {
  return singleLoopRows;
}

#ifdef REWRITE_ENABLED
  vector<int>& Analyzer::getFlopsPerLevelRewrite() {
    vector<int>& ref = flopsPerLevelRewrite;

    return ref;
  }

  int Analyzer::findMaxLevelOfPredecessors(int row) {
    int maxLevel = 0;
    vector<int>& parents = dag[row].first;

    #pragma omp parallel for reduction(max:maxLevel) 
    for(auto& parent : parents) {
      if(maxLevel < levels[parent]) 
        maxLevel = levels[parent];
    }
  
  /*  if(maxLevel == 0) {
      cout << "maxLevel is 0. Algorithm failed\n";
      cout << row << " level: " << levels[row] << "\n";
      cout << rowPtrL[row] <<  " to " << rowPtrL[row + 1] - 1 << "\n";
      for(int j = rowPtrL[row]; j < rowPtrL[row + 1] - 1; j++) 
        cout << row <<  " " << colIdxL[j] << " " << levels[colIdxL[j]] << "\n";  
    }*/
  
    return maxLevel;
  }

  void Analyzer::correctAfterRewritingStrategy(vector<int>& emptyLevels) {
    numOfLevels -= emptyLevels.size();

    for(auto& level : emptyLevels) {
      levelTable.erase(levelTable.begin()+level);
      flopsPerLevel.erase(flopsPerLevel.begin()+level);
  //    flopsPerLevelRewrite.erase(flopsPerLevelRewrite.begin()+level);
    }

 //   cout << "levelTable size: " << levelTable.size() << "\n",
    cout << "after numOfLevels: " << numOfLevels << "\n";
/*    cout << "emptyLevels:\n";
    for(auto& emptyLevel : emptyLevels)
      cout << emptyLevel << "\n";*/

    flopsPerLevelRewrite.reserve(numOfLevels);
    //fill(flopsPerLevelRewrite.begin(), flopsPerLevelRewrite.begin() + numOfLevels, 0);
    for(int i = 0 ; i < numOfLevels ; i++)
      flopsPerLevelRewrite.push_back(0);

    flopsPerLevelRewrite.shrink_to_fit();
  //  avgFLOPSPerLevel = (double)accumulate(flopsPerLevelRewrite.begin(), flopsPerLevelRewrite.end(), 0)/numOfLevels;
  //  cout << "avg. FLOPS per level after rewrite : " << avgFLOPSPerLevel << "\n";
  }
#endif

void Analyzer::buildLevels() {
  if(matrixCSC != nullptr && matrixCSR != nullptr) {
    Part* L = matrixCSC->getL();
    vector<int>& rowPtrL = L->getRowPtr();
    vector<int>& colIdxL = L->getColIdx();
    vector<double>& valsL = L->getVals();
    int cols = L->getCols();
    int rows = L->getRows();

    vector<int>& rowPtrLCSR = matrixCSR->getL()->getRowPtr();
//    vector<double>& valsLCSR = matrixCSR->getL()->getVals();

    int succLevels = posix_memalign((void **)&levels, 32, (cols-1)*sizeof(int));
    CHK_MEM_ERR(succLevels,"levels");
    dag.reserve(cols-1);
    values.reserve(cols-1);

    #pragma vector aligned
    std::fill(levels, levels + cols-1, 0);

    cout << "rows: " << rows << " cols: " << cols << "\n";
    for(int i = 0; i < cols-1; i++) {
      pair< vector<int>,vector<int> > connectivity = make_pair(vector<int>(),vector<int>());
      dag.push_back(connectivity);
      values.push_back(vector<double>());
    }

    // create levels array:
    // index: row, value: level that the row is in
    int maxLevel = levels[0];
    //for(int i = 0; i < colsLCSC-1; i++) {
    for(int i = 0; i < cols-1; i++) {
      int currLevel = levels[rowPtrL[colIdxL[i]]];

 //     if(currLevel == 0) {
 //       cout << "filling values for " << rowPtrL[colIdxL[i]] << " with " << valsL[colIdxL[i]] << "\n";
 //       vector<double>& rowValues = values[rowPtrL[colIdxL[i]]];
 //       rowValues.push_back(valsL[colIdxL[i]]);
 //     }

      Connectivity& connectivity = dag[i];
      vector<int>& children = connectivity.second;

      for(int j = colIdxL[i]+1 ; j < colIdxL[i+1]; j++) {
        if(levels[rowPtrL[j]] <= currLevel) {
          levels[rowPtrL[j]] = currLevel+1;
          if(levels[rowPtrL[j]] > maxLevel)
              maxLevel = levels[rowPtrL[j]];
        }
     
        children.push_back(rowPtrL[j]);

        Connectivity& connectivityChild = dag[rowPtrL[j]];
        vector<int>& parents = connectivityChild.first;
        parents.push_back(i);

//        cout << "next " << rowPtrL[j] << " with " << valsL[j] << "\n";
        vector<double>& rowValues = values[rowPtrL[j]];
        rowValues.push_back(valsL[j]);
      }
    }
       
    // fill in the value for the diagonal
    for(int i = 0; i < cols-1; i++) {
      vector<double>& rowValues = values[rowPtrL[colIdxL[i]]];
      rowValues.push_back(valsL[colIdxL[i]]);
    }

    // fill in the parents of the last row in DAG
    Part* LCSR = matrixCSR->getL();
    //vector<int>& rowPtrLCSR = LCSR->getRowPtr();
    vector<int>& colIdxLCSR = LCSR->getColIdx();

    Connectivity& connectivity = dag[cols-1];
    vector<int>& parents = connectivity.first;
    for(int j = rowPtrLCSR[cols-1]; j < rowPtrLCSR[cols-1 + 1] - 1; j++)
      parents.push_back(colIdxLCSR[j]);


    // create levelTable
    levelTable.reserve(maxLevel+1);
    for(int i = 0 ; i < maxLevel+1; i++)
      levelTable.push_back(vector<int>());

    for(int k = 0 ; k < cols-1 ; k++) {
      vector<int>& level = levelTable[levels[k]];
      level.push_back(k);
    }

    numOfLevels = maxLevel+1;

    levelTable.shrink_to_fit();
  } else {
    printf("Wrong format for level building. It has to be CSC\n");
  }
}

/*void Analyzer::buildEfficientCSR() {
  Part* L = matrixCSR->getL();
  vector<int>& rowPtrL = L->getRowPtr();
  vector<int>& colIdxL = L->getColIdx();
  int rows = L->getRows();
  cout << "rows: " << rows << "\n";

  efficientCSR.reserve(rows);
  for(int i = 0 ; i < rows; i++)
    efficientCSR.push_back(vector<int>());

  for(int i = 0; i < rows ; i++) {
    vector<int>& row = efficientCSR[i];
    for(int j = rowPtrL[i]; j < rowPtrL[i + 1] - 1; j++)
      row.push_back(j);
  }
}*/

void Analyzer::separateRows(int levelNum, int rowStartIndex, int rowEndIndex, vector<int>& loopedRows, vector<int>& unrolledRows) {
  // if row length (a.k.a number of parents < 10) put them in unrolledRows
  // rest is in loopedRows, will be executed in 1 loop as level-set
  vector<int>::iterator it = loopedRows.begin();
  while(it != loopedRows.end()) {
    if(dag[*it].first.size() < 5) {
      unrolledRows.push_back(*it);
      loopedRows.erase(it);
    } else it++;
  }
}

void Analyzer::calculateFLOPS() {
  /*Part* L = matrixCSR->getL();
  vector<int>& rowPtrL = L->getRowPtr();
  vector<int>& colIdxL = L->getColIdx();*/

  if(!flopsPerLevel.empty())
    flopsPerLevel.clear();

  flopsPerLevel.reserve(numOfLevels);

//  cout << "FLOPS per row:\n";
  singleLoopRows = false;
  for(int i = 0; i < numOfLevels; i++) {
    int indegreeSum = 0;
    int lengthOverThreshold = 0; // threshold = 5
    vector<int>& level = levelTable[i];
//    cout << "level, " << i << "\n";
    for(auto& row : level) {
      if(dag[row].first.size() > 4)
        lengthOverThreshold++;
      //cout << (dag[row].first.size() << 1) + 1 << "," ;
      indegreeSum += dag[row].first.size()+1;
      //indegreeSum += (rowPtrL[row + 1] - rowPtrL[row]);
    }
//    cout << "\n";
    if(lengthOverThreshold > 0)
      singleLoopRows = true;

    flopsPerLevel.push_back((indegreeSum << 1) - level.size());

    indegreeSum = 0;
 //   cout << "singleLoopRows: " << singleLoopRows << "\n";
  }

  flopsPerLevel.shrink_to_fit();

  double totalFLOPSPerLevel = accumulate(flopsPerLevel.begin(), flopsPerLevel.end(), 0.0);
  avgFLOPSPerLevel = totalFLOPSPerLevel/numOfLevels;
}

#ifdef REWRITE_ENABLED
  void Analyzer::calculateLevelsToBeRewritten() {
    if(flopsBelowAvg.empty() && flopsAboveAvg.empty()) {
      for(int i = 0; i < numOfLevels; i++) {
        if(flopsPerLevel[i] < avgFLOPSPerLevel)
           flopsBelowAvg[i] = flopsPerLevel[i];
        else 
           flopsAboveAvg[i] = flopsPerLevel[i];
      }
    } else
      cout << "already calculated\n";

    cout << "avg. FLOPS per Level: " << avgFLOPSPerLevel << "\n"; 
    cout << "flopsBelowAvg:\n";
    for(auto& level : flopsBelowAvg)
      cout << level.first << " : " << level.second << "\n";
  
    cout << "flopsAboveAvg:\n";
    for(auto& level : flopsAboveAvg)
      cout << level.first << " : " << level.second << "\n";
  }
#endif

/*void Analyzer::printRowHist() {
  int rows = matrixCSR->getNumOfRows();
  cout << "row histogram:\n";
  for(int i = 0 ; i < rows ; i++)
    cout << i << " : " << rowHist[i] << "\n";
  cout << "\n";
}*/

void Analyzer::printLevels() {
  if(matrixCSC != nullptr) {
    int rows = matrixCSC->getNumOfRows();

    cout << "num of rows: " << rows << "\n";
    for(int i = 0 ; i < rows ; i++)
      cout << i << ": " << levels[i] << "\n";
    cout << "\n";
  } else {
    printf("Wrong format while printing levels. It has to be CSC\n");
  }
}

void Analyzer::printLevelTable() {
  if(matrixCSC != nullptr) {
    for(auto& level : levelTable) {
      for(auto& row : level)
        cout << row << ", ";
      cout << "\n";
    }
  } else {
    printf("Wrong format while printing level table.\n");
  }
}

void Analyzer::printLevelSizes(){
  cout << "num levels: " << levelTable.size() << "\n";
  cout << "new level sizes:\n";
  if(matrixCSC != nullptr) {
    for(auto& level : levelTable)
      cout << level.size() << "\n";
  } else {
    printf("Wrong format while printing level table.\n");
  }
}

void Analyzer::printDAG() {
  cout << "DAG:\n";
  for(int i = 0 ; i < dag.size(); i++) {
    cout << "row " << i << ":\n";
    Connectivity& connectivity = dag[i];
    vector<int>& parents = connectivity.first;
    vector<int>& children = connectivity.second;

    cout << "parents:\n";
    for(auto& parent : parents)
      cout << parent << ", ";
    cout << "\n\n";

    cout << "children:\n";
    for(auto& child : children)
      cout << child << ", ";
    cout << "\n\n";
  }
}

void Analyzer::printValues() {
  cout << "Values:\n";
  for(int i = 0 ; i < values.size(); i++) {
    cout << "row " << i << ":\n";
    for(auto& value : values[i]) {
      cout << value << ", ";
    cout << "\n";
    }

    cout << "\n\n";
  }
}

void Analyzer::printFLOPSPerLevel() {
  cout << "FLOPS per Level:\n";
  for(auto& level : flopsPerLevel)
    cout << level << "\n";
  cout << "\n";
}

void Analyzer::reportBefore() {
  double totalFLOPSPerLevel = accumulate(flopsPerLevel.begin(), flopsPerLevel.end(), 0.0);
  avgFLOPSPerLevel = totalFLOPSPerLevel/numOfLevels;
  cout << "total FLOPS: " << totalFLOPSPerLevel << "\n";
  cout << "avg. FLOPS per level: " << avgFLOPSPerLevel << "\n";
  printFLOPSPerLevel();
}


#ifdef REWRITE_ENABLED
  void Analyzer::printFLOPSPerLevelRewrite() {
    cout << "FLOPS per Level:\n";
    for(auto& level : flopsPerLevelRewrite)
      cout << level << "\n";
    cout << "\n";
  }

  void Analyzer::reportAfter() {
    double totalFLOPSPerLevel = accumulate(flopsPerLevelRewrite.begin(), flopsPerLevelRewrite.end(), 0.0);
    avgFLOPSPerLevel = totalFLOPSPerLevel/numOfLevels;
    cout << "total FLOPS: " << totalFLOPSPerLevel << "\n";
    cout << "avg. FLOPS per level: " << avgFLOPSPerLevel << "\n";
    printFLOPSPerLevelRewrite();
  }

  // helper function to validate if rewriting is done correctly
  void Analyzer::printDependencies() {
    Part* L = matrixCSR->getL();
  /*  vector<int>& rowPtrL = L->getRowPtr();
    vector<int>& colIdxL = L->getColIdx();*/
    int rows = L->getRows();
    cout << "rows: " << rows << "\n";
  
    for(int i = 0; i < rows-1 ; i++) {
      vector<int>& parents = dag[i].first;
      cout << "row " << i << " has " << parents.size() << " dependencies:\n";
      for(auto& parent : parents) {
        cout << parent << " ,";
      cout << "\n";
      }
    }
  }
#endif

