#pragma once
#include <algorithm>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <iostream>
#include <omp.h>
#include <shared_mutex>
#include <thread>
#include "analyzer.h"

#define REWRITE_UP    3
#define REWRITE_LOW   5
#define REWRITE_RATIO 0.05
//#define UPDATE_GRAPES  

// add random
enum class StartPoint {
  TopDown, BottomUp   };

enum class Scope   {
  Depth, Selective };

using namespace std;

typedef map<int,pair<int,int>> ToBeRewritten;
typedef map<int,set<int>> RewritingMap;
typedef map<int,vector<int>> InitialParents;

// indegree & outdegree vectors of a row: dependencies and children
typedef pair< vector<int>,vector<int> > Connectivity;
typedef vector<Connectivity> DAG;


// TODO: alignment
class RewritingStrategy {
  protected:
    // <row,<oldLevel,newLevel>>
    ToBeRewritten  toBeRewritten;
    RewritingMap rewritingMap;
    InitialParents initialParents;
    int numOfRewrittenRows;
    int* levels;
    vector< vector<int> >& levelTable;
    DAG& dag;

    int rows;
    StartPoint startPoint;
    Scope scope;


  public:
    RewritingStrategy(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) :
                      levels(levels), levelTable(levelTable), dag(dag) {
      this->rows = rows;
      this->startPoint = startPoint;
      this->scope = Scope::Selective;
    }

    virtual void policy() = 0;

    ToBeRewritten& getToBeRewritten() {
      ToBeRewritten& ref = toBeRewritten;

      return ref;
    }

    RewritingMap& getRewritingMap() {
      RewritingMap& ref = rewritingMap;

      return ref;
    }

    void printRowsToBeRewritten() {
      cout << "Rows to be rewritten:\n";
      for(auto& row : toBeRewritten)
        cout << "row : " << row.first << " at level " << row.second.first << " to " << row.second.second << "\n";
      cout << "\n";
    }

    void printRewritingMap() {
      cout << "Rewriting map for rows to be rewritten\n";
      for(auto& row : rewritingMap) {
        cout << "row: " << row.first << "\n";
        for(auto& item : row.second)
          cout << item << ",";
        cout << "\n";
      }

      cout << "\n\n";
    }

    void printInitialParents() {
      cout << "Initial parents:\n";
      for(auto& row : initialParents) {
        cout << "row : " << row.first << "\n";
        for(auto& parent : row.second)
          cout << parent << ",";
        cout << "\n";
      }

      cout << "\n\n";
    }

    bool isRewritten(int row) {
      return toBeRewritten.find(row) != toBeRewritten.end();
    }

    bool isTopDown() {
      return (startPoint == StartPoint::TopDown);
    }

    bool isScopeSelective() {
      return (scope == Scope::Selective);
    }

    set<int>& getRewritingMapOf(int row) {
      set<int>& ref = rewritingMap[row];

      return ref;
    }

    // duplicate function with Analyzer class
    int findMaxLevelOfPredecessors(int row) {
      int maxLevel = 0;
      vector<int>& parents = dag[row].first;

      #pragma omp parallel for reduction(max:maxLevel)
       for(auto& parent : parents) {
         if(maxLevel < levels[parent])
           maxLevel = levels[parent];
       }

       return maxLevel;
    }

    void shiftRowsUp() {
      if(startPoint == StartPoint::TopDown) {
        for(auto& row : toBeRewritten) {
       //   cout << "shifting up " << row.first << "\n";
          initialParents[row.first] = vector<int>(dag[row.first].first);
          createRewritingMapFor(row.first, toBeRewritten[row.first].second);
          updateLevelOf(row.first, toBeRewritten[row.first].second);
          #ifdef UPDATE_GRAPES
            updateGrapeRows(row.first);
          #endif
        }
      } else {
     //   cout << "bottom up\n";
        for(auto it = toBeRewritten.rbegin() ; it != toBeRewritten.rend() ; ++it) {
     //     cout << "shifting up " << it->first << "\n";
          initialParents[it->first] = vector<int>(dag[it->first].first);
          createRewritingMapFor(it->first, toBeRewritten[it->first].second);
          updateLevelOf(it->first, toBeRewritten[it->first].second);
          #ifdef UPDATE_GRAPES
            updateGrapeRows(it->first);
          #endif
        }
      }
    }

    void updateLevelOf(int row, int newLevel) {
      int oldLevel = levels[row];
      levels[row] = newLevel;
  
      vector<int>& oLevel = levelTable[oldLevel];
      vector<int>& nLevel = levelTable[newLevel];
  
      vector<int>::iterator it = lower_bound(oLevel.begin(), oLevel.end(), row);
      // TODO: to optimize erase (due to shift) swap with the last element and erase it if the order doesnt matter
      if(it != oLevel.end())
        oLevel.erase(it);
      else
        cout << "couldnt erase row " << row << " at level " << oldLevel << "\n";
  
      it = upper_bound(nLevel.begin(), nLevel.end(), row);
      if(it == nLevel.end())
        nLevel.push_back(row);
      else
        nLevel.insert(it,row);
    }

    virtual void updateGrapeRows(int row) {
  //    cout << "update grape rows for " << row << "\n";
      vector<int>& children = dag[row].second;
     
      vector<int> toBeErased;
      for(int i = 0 ; i < children.size(); i++) {
        int childLevel = levels[children[i]];

//        cout << "child in question: " << children[i] << "\n";
        // child is already rewritten
        if(childLevel < levels[row]+1) {
          // parent (i.e. row) is not a parent of this child anymore (due to rewriting) remove it from children list
          if(find(dag[children[i]].first.begin(), dag[children[i]].first.end(), row) == dag[children[i]].first.end()) {
  //          cout << children[i] << " will be erased\n";
            toBeErased.push_back(i);
          }
        } else {
          int maxLevel = findMaxLevelOfPredecessors(children[i]);
   //       cout << "child : " << children[i] << " @ level: " << childLevel << " maxLevel of predecessors: " << maxLevel << "\n";
    
          // this is a grape to be shifted up
          int targetLevel = maxLevel+1;
          if(targetLevel < childLevel) {
            updateLevelOf(children[i], targetLevel);
            updateGrapeRows(children[i]);
          }
        }
      }
    }

    void createRewritingMapFor(int row, int targetLevel) {
     /*  cout << "createRewritingMapFor:\n";
       cout << "row: " << row << " targetLevel: " << targetLevel << "\n";*/
   
       // TODO: dag[row].first will have duplicates how about rewritten as the new dag[row].first
       vector<int> predsWithMaxLevel;
       int maxLevel = findPredecessorsWithMaxLevel(dag[row].first, predsWithMaxLevel);
       rewritingMap[row] = set<int>();
   
   //    cout << "maxLevel: " << maxLevel << " targetLevel: " << targetLevel << "\n";

       if(maxLevel == targetLevel)  { // rewriting to the previous level
         for(auto& pred : predsWithMaxLevel) {
     //      cout << "pred: " << pred << "\n";
           expandPredsWith(pred, dag[row].first, row);
         }
       }

       while(maxLevel > targetLevel) {
   //      cout << "maxLevel: " << maxLevel << " targetLevel: " << targetLevel << "\n";
         for(auto& pred : predsWithMaxLevel) {
  //         cout << "pred: " << pred << "\n";
           expandPredsWith(pred, dag[row].first, row);
         }
   
         rewritingMap[row].insert(predsWithMaxLevel.begin(), predsWithMaxLevel.end());
         predsWithMaxLevel.clear();
         maxLevel = findPredecessorsWithMaxLevel(dag[row].first, predsWithMaxLevel);
      }
         
      rewritingMap[row].insert(predsWithMaxLevel.begin(), predsWithMaxLevel.end());
    }

    void expandPredsWith(int row, vector<int>& preds, int child) {
      auto it = lower_bound(preds.begin(),preds.end(),row);
      // Replace row with first pred. of it instead of erasing and causing a shift.
      // Then push back the rest since order doesnt matter.
      if(it != preds.end()) {
        vector<int>& parents = dag[row].first;

       /* cout << "parents of " << row << ":\n";
        for(auto& parent: parents)
          cout << parent << ", ";
        cout << "\n";

        cout << "inserting " << parents[0] << " to slot " << it-preds.begin() << " in place of " << preds[it-preds.begin()] << "\n"; */
        if(!parents.empty()) { 
          //cout << "parents empty for row " << row << "\n";
          preds[it-preds.begin()] = parents[0];
          preds.insert(preds.end(),parents.begin()+1,parents.end());
        }

        sort(preds.begin(), preds.end());
        auto it = unique(preds.begin(), preds.end());
        preds.resize(std::distance(preds.begin(),it));

        vector<int>& children = dag[row].second;
  /*      cout << "searching for " << child << " in children of " << row << " to remove\n";
        for(auto& child: children)
          cout << child << ", ";
        cout << "\n";  */
        it = find(children.begin(), children.end(), child);
        if(it != children.end())
        children.erase(it);

      /*  cout << "new preds:\n";
        for(auto& newPred: preds)
          cout << newPred << ", ";
        cout << "\n";*/
      }
      else
        cout << "row " << row << " doesnt exist in predecessors\n";
    }

    int findPredecessorsWithMaxLevel(vector<int>& preds, vector<int>& predsWithMaxLevel) {
      int maxLevel = 0;
      #pragma omp parallel for reduction(max:maxLevel)
      for(int j = 0; j < preds.size(); j++) {
        if(maxLevel < levels[preds[j]])
          maxLevel = levels[preds[j]];
      }
  
      for(auto& pred : preds)
        if(levels[pred] == maxLevel)
          predsWithMaxLevel.push_back(pred);
  
      return maxLevel;
    }


    // after rewriting strategy is applied, some levels will be emptied.
    // reflect the changes to related data structures (e.g. levelTable)
    void findEmptyLevels(vector<int>& emptyLevels) {
      int numOfLevels = levelTable.size();
      cout << "before numOfLevels: " << numOfLevels << "\n";
      for(vector<vector<int>>::reverse_iterator rit = levelTable.rbegin() ; rit != levelTable.rend() ; ++rit) {
        if(rit->empty()) {
       //   cout << "empty level: " << numOfLevels-1-(rit-levelTable.rbegin()) << "\n";
          emptyLevels.push_back(numOfLevels-1-(rit-levelTable.rbegin()));
        }
      }
    }

    vector<int>& getInitialParentsOf(int row) {
      vector<int>& ref = initialParents[row];

      return ref;
    }
};

class RewriteLongerThan : public RewritingStrategy {
  public:
    RewriteLongerThan(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) : 
                      RewritingStrategy(rows,startPoint,levels,levelTable,dag) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy() {
      if(startPoint == StartPoint::TopDown)
        for(int k = 0 ; k < rows ; k++) {
          if(dag[k].first.size() >= REWRITE_UP &&
              levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }
      else {
        for(int k = rows-1 ; k >= 0 ; k--) {
          if(dag[k].first.size() >= REWRITE_UP &&
             levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }
      }

      // we intend to rewrite rows * REWRITE_RATIO rows. The matrix might not have
      // that many rows with the criteria (>= REWRITE_UP) we set.
      numOfRewrittenRows = toBeRewritten.size();
      cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
    }
};

class RewriteShorterThan : public RewritingStrategy {
  public:
    RewriteShorterThan(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) : 
                      RewritingStrategy(rows,startPoint,levels,levelTable,dag) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy() {
      if(startPoint == StartPoint::TopDown)
        for(int k = 0 ; k < rows ; k++) {
          if(dag[k].first.size() <= REWRITE_UP &&
             levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }
      else
        for(int k = rows-1 ; k >= 0 ; k--) {
          if(dag[k].first.size() <= REWRITE_UP &&
             levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }

      numOfRewrittenRows = toBeRewritten.size();
      cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
    }
};

class RewriteInBetween : public RewritingStrategy {
  public:
    RewriteInBetween(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) : 
                      RewritingStrategy(rows,startPoint,levels,levelTable,dag) {
      numOfRewrittenRows = rows * REWRITE_RATIO;
      cout << "intended num. of rewritten rows: " << numOfRewrittenRows << "\n";

      policy();
    }

    void policy() {
      if(startPoint == StartPoint::TopDown)
        for(int k = 0 ; k < rows ; k++) {
          if(dag[k].first.size() >= REWRITE_UP && dag[k].first.size() <= REWRITE_LOW &&
             levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }
      else {
        for(int k = rows-1 ; k >= 0 ; k--) {
          if(dag[k].first.size() >= REWRITE_UP && dag[k].first.size() <= REWRITE_LOW &&
             levels[k] > 1 && numOfRewrittenRows-- > 0)
             toBeRewritten[k] = make_pair(levels[k],-1);
        }
      }

      numOfRewrittenRows = toBeRewritten.size();
      cout << "actual num. of rewritten rows: " << toBeRewritten.size() << "\n";
    }
};

class RewriteByCostMap : public RewritingStrategy {
  protected:
    // cost map: <level, <row, rowCost>>
    // will only keep the costs of the newly inserted rows
    map< int,map<int,int> > costMap;
    vector<int> levelCost;

    double avgCostPerLevel;
    map<int,double> flopsBelowAvg;


  public:
    RewriteByCostMap(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag, vector<int>& flopsPerLevel, double avgFLOPSPerLevel, map<int,double> flopsBelowAvg) :
                      RewritingStrategy(rows, startPoint, levels, levelTable, dag) {
       // copy the contents so that we can work on it
       levelCost = flopsPerLevel;

       this->avgCostPerLevel = avgFLOPSPerLevel;
       this->flopsBelowAvg = flopsBelowAvg;

       cout << "before policy\n";
       policy();

       flopsPerLevel = levelCost;

      /* for(int i = 0; i < levelTable.size(); i++) {
         vector<int>& level = levelTable[i];
         cout << "level " << i << ":\n";
         for(auto& row : level)
           cout << row << ", ";

         if(costMap.find(i) != costMap.end()) {
           map<int, int>& rows = costMap[i];
           for(auto& row : rows)
             cout << row.first << ", ";
           cout << "\n";
         }
       }*/

       cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
       printRowsToBeRewritten();
    }
    
    void policy() {
      // we can drag a level's cost up to the cost of level with the maximum cost in the original graph 
      //int maxCostPerLevel = *max_element(begin(levelCost),end(levelCost));

      // choose a level to start the process
//      auto rit_next = flopsBelowAvg.rbegin(); rit_next++;
      for(auto rit = flopsBelowAvg.rbegin(); rit != flopsBelowAvg.rend(); rit++) {
    //    createCostMapMin(rit->first, rit_next->first);
    //    rit_next++;
        createCostMapMax(rit->first, flopsBelowAvg.begin()->first, avgCostPerLevel);
      }

      prune(avgCostPerLevel);
      //pruneReverse(avgCostPerLevel);


 /*      cout << "levelCost:\n";
       for(auto& flop : levelCost)
         cout << flop<< " ";
       cout << "\n";*/

       // printing just the rows at each level in the costMap
     /*  for(auto it = costMap.begin(); it !=  costMap.end(); ++it) {
         cout << "level " << it->first << ":\n";
         map<int, int> rows = it->second;
         for(auto& row : rows)
           cout << row.first << ", ";
         cout << "\n";
       }*/
    }

    // 2 ways:
    // move a row up the furthest it can go   : createCostMapMax
    // move rows up level by level            : createCostMapMin
    void createCostMapMax(int level, int levelEnd, int maxCostPerLevel) {
      vector<int>& levelRows = levelTable[level];

  //    cout << "level: " << level << "\n";
      for(auto& row : levelRows) {
         vector<int> predsWithMaxLevel;
    
 //        cout << "row: " << row << "\n";
         vector<int> parents = dag[row].first;
         int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

         while(maxLevel >= levelEnd) {
           for(auto& pred : predsWithMaxLevel)
             expandPredsWith(pred, parents, row);

           int numOfPreds = parents.size();
           int costRow = numOfPreds <= 4 ? (numOfPreds << 1) + 1 : (numOfPreds << 1);

           if(flopsBelowAvg.find(maxLevel) != flopsBelowAvg.end()) {
             levelCost[maxLevel] += costRow;

             map<int,int>& newLevel = costMap[maxLevel];
             newLevel[row] = costRow;

   //          cout << "adding to level " << maxLevel << "\n";
           }

           predsWithMaxLevel.clear();
           int oldLevel = maxLevel;
           maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
           if(oldLevel == maxLevel) break;
         }
      }
    }

    // version of expandPredsWith that does not alter the dag
    void expandPredsWith(int row, vector<int>& preds, int child) {
      auto it = lower_bound(preds.begin(),preds.end(),row);
      // Replace row with first pred. of it instead of erasing and causing a shift.
      // Then push back the rest since order doesnt matter.
      if(it != preds.end()) {
        vector<int>& parents = dag[row].first;

        if(!parents.empty()) { 
          preds[it-preds.begin()] = parents[0];
          preds.insert(preds.end(),parents.begin()+1,parents.end());
        }

        sort(preds.begin(), preds.end());
        auto it = unique(preds.begin(), preds.end());
        preds.resize(std::distance(preds.begin(),it));
      }
      else
        cout << "row " << row << " doesnt exist in predecessors\n";
    }

    void prune(int avgCostPerLevel) {
      for(auto it = flopsBelowAvg.begin(); it != flopsBelowAvg.end(); it++) {
        // TODO: add index distance as a stopping criteria
        map<int,int>& targetLevel = costMap[it->first];

        int currCost = it->second;
        map<int,int>::iterator itr = targetLevel.end();
        for(auto rewrittenRow = targetLevel.begin(); rewrittenRow != targetLevel.end(); rewrittenRow++) {
          if((currCost + rewrittenRow->second) < avgCostPerLevel) {
            currCost += rewrittenRow->second;
            toBeRewritten[rewrittenRow->first] = make_pair(levels[rewrittenRow->first],it->first);
          } else {
            itr = rewrittenRow;
            break;
          }
        }

        // now remove the remaining from the current level
        if(itr != targetLevel.end()) {
          for( auto it_target = itr; it_target != targetLevel.end(); it_target++)
            levelCost[it->first] -= it_target->second;

          targetLevel.erase(itr, targetLevel.end());
        }

        // now remove these rows from the following levels
        for(auto iter = next(it); iter != flopsBelowAvg.end() ; iter++) {
          map<int,int>& currLevel = costMap[iter->first];
          for(auto iter_targetLevel = targetLevel.begin();  iter_targetLevel != targetLevel.end(); iter_targetLevel++) {
            auto it_erase = currLevel.find(iter_targetLevel->first);
            if(it_erase != currLevel.end()) {
              currLevel.erase(it_erase);
              levelCost[iter->first] -= it_erase->second;
            }
          }
        } // for
      }
    }

    void pruneReverse(int avgCostPerLevel) {
      for(auto it = flopsBelowAvg.begin(); it != flopsBelowAvg.end(); it++) {
        // TODO: add index distance as a stopping criteria
        map<int,int>& targetLevel = costMap[it->first];

/*        cout << "level: " << it->first << "\n";
        cout << "costMap:\n";
        for(auto& row : targetLevel)
          cout << row.first << ", ";
        cout << "\n";*/

        int currCost = it->second;
        map<int,int>::reverse_iterator itr = targetLevel.rend();
        for(auto rewrittenRow = targetLevel.rbegin(); rewrittenRow != targetLevel.rend(); rewrittenRow++) {
          if((currCost + rewrittenRow->second) < avgCostPerLevel) {
 //           cout << "row " << rewrittenRow->first << " stays\n";
            currCost += rewrittenRow->second;
            toBeRewritten[rewrittenRow->first] = make_pair(levels[rewrittenRow->first],it->first);
          } else {
            itr = rewrittenRow;
            break;
          }
        }

        if(itr != targetLevel.rend())
          cout << "stopped at row " << itr->first << "\n";

        // now remove the remaining from the current level
        auto pivot = targetLevel.begin(); 
        if(itr != targetLevel.rend()) {
          for( auto it_target = targetLevel.rbegin(); it_target != itr; it_target++)
            levelCost[it->first] -= it_target->second;

          int dist = distance(begin(targetLevel), itr.base())-1;
          advance(pivot, dist);

          targetLevel.erase(targetLevel.begin(), pivot);
        }

  /*      cout << "after erasing costMap:\n";
        for(auto& row : targetLevel)
          cout << row.first << ", ";
        cout << "\n";*/

        // now remove these rows from the following levels
        for(auto iter = next(it); iter != flopsBelowAvg.end() ; iter++) {
          map<int,int>& currLevel = costMap[iter->first];
          for(auto iter_targetLevel = pivot;  iter_targetLevel != targetLevel.end(); iter_targetLevel++) {
            auto it_erase = currLevel.find(iter_targetLevel->first);
            if(it_erase != currLevel.end()) {
              currLevel.erase(it_erase);
              levelCost[iter->first] -= it_erase->second;
            }
          }
        } // for
      }
    }
    void createCostMapMin(int level, int lastLevelSeen) {
      // do for original rows of that level + rows in the costMap
      vector<int>& levelRows = levelTable[level];
      map<int,int>& costMapRows = costMap[level];

      for(auto& row : levelRows) {
         vector<int> predsWithMaxLevel;
     
         vector<int> parents = dag[row].first;
         int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

         while(maxLevel >= lastLevelSeen) {
           for(auto& pred : predsWithMaxLevel)
             expandPredsWith(pred, parents, row);
  
           int numOfPreds = parents.size();
           int costRow = numOfPreds <= 4 ? (numOfPreds << 1) + 1 : (numOfPreds << 1);
  
           if(flopsBelowAvg.find(maxLevel) != flopsBelowAvg.end() &&
              maxLevel != level) {
             levelCost[maxLevel] += costRow;
 
             map<int,int>& newLevel = costMap[maxLevel];
             newLevel[row] = costRow;
           }
         
           predsWithMaxLevel.clear();
           maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
           if(lastLevelSeen == 0) break;
         }
      }

      for(auto& costMapRow : costMapRows) {
        int row = costMapRow.first;
         vector<int> predsWithMaxLevel;
     
         vector<int> parents = dag[row].first;
         int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

         while(maxLevel >= lastLevelSeen) {
           for(auto& pred : predsWithMaxLevel)
             expandPredsWith(pred, parents, row);
  
           int numOfPreds = parents.size();
           int costRow = numOfPreds <= 4 ? (numOfPreds << 1) + 1 : (numOfPreds << 1);
  
           if(flopsBelowAvg.find(maxLevel) != flopsBelowAvg.end() &&
              maxLevel != level) {
             levelCost[maxLevel] += costRow;
  
             map<int,int>& newLevel = costMap[maxLevel];
             newLevel[row] = costRow;
           }

           predsWithMaxLevel.clear();
           maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
           if(lastLevelSeen == 0) break;
         }
      }
    }
};

// methods using RewriteByLevel
// must traverse rowsToBeRewritten in REVERSE order
// hence RewriteByLevel is BottomUp by nature
class RewriteByLevel : public RewritingStrategy {
  protected:
    map<int,vector<int>> levelLengthHistogram;
    vector<int> levelsToBeRewritten;

  public:
    RewriteByLevel(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag, map<int,double>& flopsBelowAvg) : 
                      RewritingStrategy(rows, startPoint, levels, levelTable, dag) {
      prepare(flopsBelowAvg);    
//      policy();

//        lung2_7levels_to8th();
        lung2_9levels_to10th();

        cout << "levelsToBeRewritten:\n";
        for(auto& level : levelsToBeRewritten)
          cout << level << ", ";
        cout << "\n";

        cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
        printRowsToBeRewritten();
    }

    
    // currently we're merging all levels with the min. # of rows for simplicity
    void policy() override {
      auto it = levelLengthHistogram.begin();
      for(auto& level : it->second)
        levelsToBeRewritten.push_back(level);

      // if there's only 1 level with the min. # of rows, check for levels with the second min. # of rows
      if(levelsToBeRewritten.size() == 1) {
        it++;
        for(auto& level : it->second)
          levelsToBeRewritten.push_back(level);
      }

      sort(levelsToBeRewritten.begin(),levelsToBeRewritten.end());

/*      cout << "levelsToBeRewritten:\n";
      for(auto& level : levelsToBeRewritten)
        cout << level << ", ";
      cout << "\n";*/

      // auto itTargetLevel = levelsToBeRewritten.rbegin();
      auto itTargetLevel = levelsToBeRewritten.begin();
      // skip the first level since it's our targetLevel
      auto itrback = levelsToBeRewritten.rend();
      itrback--;
      for(auto it = levelsToBeRewritten.rbegin() ; it != itrback; it++) {
        // together with commented line before this for loop where itTargetLevel is set to levelsToBeRewritten.rbegin()
        // this if condition makes each level to be shifted to the previous level in levelsToBeRewritten
        /*if(++itTargetLevel == levelsToBeRewritten.rend())
          break; */

        for(auto& row : levelTable[*it])
          toBeRewritten[row] = make_pair(levels[row],*itTargetLevel); 
      }

      // TODO: just to debug. remove this!
/*      toBeRewritten.erase(19);
      toBeRewritten.erase(28);
      toBeRewritten[29].second = 11;  // targetLevel is a level to be rewritten*/
  /*    toBeRewritten.erase(1);
//      toBeRewritten.erase(28);
      toBeRewritten[28].second = 9;  // targetLevel is a level to be rewritten
      toBeRewritten[19].second = 8;  // targetLevel is a level to be rewritten
      toBeRewritten.erase(29);  */

        cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
        printRowsToBeRewritten();
    }

    void prepare(map<int,double>& flopsBelowAvg) {
      // create a histogram of num. of levels with each num. of rows in flopsBelowAvg
      // merge all levels with the minimum num. of rows within themselves 
      for(auto it = flopsBelowAvg.begin() ; it != flopsBelowAvg.end(); ++it) {
        vector<int>& entry = levelLengthHistogram[levelTable[it->first].size()]; // key to histogram: level's size
        entry.push_back(it->first);
      }

    /*  cout << "levelLengthHistogram:\n";
      for(auto& levelSize : levelLengthHistogram) {
        vector<int>& levels = levelSize.second;
        cout << "size: " << levelSize.first << "\n";
        for(auto& level : levels)
          cout << level << ", ";
        cout << "\n";
      }*/
    }

    void lung2_9levels_to10th() {
      // lung2 - 9 levels written to 10th
      levelsToBeRewritten.push_back(1);
      levelsToBeRewritten.push_back(2);
      int targetRow = 12;
      for(int i = 3; i < 452; i++) {
        if(i != targetRow)
          levelsToBeRewritten.push_back(i);
        else
          targetRow+=10;
      }
     
      levelsToBeRewritten.push_back(453);
      levelsToBeRewritten.push_back(454);
      levelsToBeRewritten.push_back(455);
      levelsToBeRewritten.push_back(477);
      levelsToBeRewritten.push_back(478);

      for(auto& row : levelTable[478])
         toBeRewritten[row] = make_pair(478,476);
      for(auto& row : levelTable[477])
         toBeRewritten[row] = make_pair(477,476);
      for(auto& row : levelTable[455])
         toBeRewritten[row] = make_pair(455,452);
      for(auto& row : levelTable[454])
         toBeRewritten[row] = make_pair(454,452);
      for(auto& row : levelTable[453])
         toBeRewritten[row] = make_pair(453,452);

      targetRow = 442;
      for(int i = 451; i > 11; i--) {
        if(targetRow != i)
          for(auto& row : levelTable[i])
            toBeRewritten[row] = make_pair(i,targetRow);
        else
          targetRow-=10;
      }

      for(int i = 11; i > 0; i--)
        for(auto& row : levelTable[i])
           toBeRewritten[row] = make_pair(i,0);
    }

    void lung2_7levels_to8th() {
      // lung2 - 7 levels written to 8th
      levelsToBeRewritten.push_back(1);
      levelsToBeRewritten.push_back(2);
      levelsToBeRewritten.push_back(3);
      int targetRow = 10;
      for(int i = 5; i < 452; i++) {
        if(i != targetRow)
          levelsToBeRewritten.push_back(i);
        else
          targetRow+=8;
      }

      levelsToBeRewritten.push_back(453);
      levelsToBeRewritten.push_back(454);
      levelsToBeRewritten.push_back(455);
      levelsToBeRewritten.push_back(478);

      for(auto& row : levelTable[478])
         toBeRewritten[row] = make_pair(478,477);
      for(auto& row : levelTable[455])
         toBeRewritten[row] = make_pair(455,452);
      for(auto& row : levelTable[454])
         toBeRewritten[row] = make_pair(454,452);
      for(auto& row : levelTable[453])
         toBeRewritten[row] = make_pair(453,452);

      targetRow = 444;
      for(int i = 451; i > 3; i--) {
        if(targetRow != i)
          for(auto& row : levelTable[i])
            toBeRewritten[row] = make_pair(i,targetRow);
        else
          targetRow-=8;
      }
          
      for(auto& row : levelTable[3])
         toBeRewritten[row] = make_pair(3,0);
      for(auto& row : levelTable[2])
         toBeRewritten[row] = make_pair(2,0);
      for(auto& row : levelTable[1])
         toBeRewritten[row] = make_pair(1,0);
    }

    void updateGrapeRows(int row) override {
 //     cout << "update grape rows for " << row << "\n";
      vector<int>& children = dag[row].second;
     
      vector<int> toBeErased;
      for(int i = 0 ; i < children.size(); i++) {
        int childLevel = levels[children[i]];

  //      cout << "child in question: " << children[i] << "\n";
        // child is already rewritten
        if(childLevel < levels[row]+1) {
          // parent (i.e. row) is not a parent of this child anymore (due to rewriting) remove it from children list
          if(find(dag[children[i]].first.begin(), dag[children[i]].first.end(), row) == dag[children[i]].first.end()) {
 //           cout << children[i] << " will be erased\n";
            toBeErased.push_back(i);
          }
        } else {
          int maxLevel = findMaxLevelOfPredecessors(children[i]);
  //        cout << "child : " << children[i] << " @ level: " << childLevel << " maxLevel of predecessors: " << maxLevel << "\n";
    
          // this is a grape to be shifted up
          int targetLevel = maxLevel+1;
          if(targetLevel < childLevel) {
  //          cout << "targetLevel: " << targetLevel << " < chilLevel " << childLevel << "\n";
            bool update = true;
//            vector<int>& levelsToBeRewritten = levelLengthHistogram.begin()->second;
            bool isALevelToBeRewritten = (find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), targetLevel) != levelsToBeRewritten.end());
            if(isALevelToBeRewritten) {
  //            cout << children[i] << " will be in a level to be rewritten: " << targetLevel << "\n";
                  
              while(find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), targetLevel) != levelsToBeRewritten.end() &&
                    targetLevel > childLevel)
                targetLevel++;
                  
    //          cout << " targetLevel adjusted: " << targetLevel <<  " childLevel: " << childLevel << "\n";

              // already scheduled for rewriting, update
              if(toBeRewritten.find(children[i]) != toBeRewritten.end()) {
                // cout << children[i] << " is already scheduled to be rewritten\n";
                // if the new level is a level to be rewritten, we cannot shift the row there
                // instead just shift it to the nearest level that wont be rewritten
                if((targetLevel == childLevel) || (targetLevel == toBeRewritten[children[i]].second)) {
                  toBeRewritten.erase(children[i]);    
                  update = false;
                } else if(targetLevel < toBeRewritten[children[i]].second) {
                  // cout << "targetLevel : " << targetLevel << " < " << toBeRewritten[children[i]].second << " (scheduled level\n)";  
                  toBeRewritten.erase(children[i]);    

                } else if(targetLevel > toBeRewritten[children[i]].second) {
                  toBeRewritten[children[i]].first = targetLevel;
                  // cout << "new level assignment:\n";
                  // cout << toBeRewritten[children[i]].first << " : " << toBeRewritten[children[i]].second << "\n";
                }
              } else {
                if(targetLevel == childLevel)
                  update = false;
              }
            } else {
              if(toBeRewritten.find(children[i]) != toBeRewritten.end())
                if(targetLevel <= toBeRewritten[children[i]].second)
                  toBeRewritten.erase(children[i]);    
            }

            if(update) {
              updateLevelOf(children[i], targetLevel);
              updateGrapeRows(children[i]);
            }
          }
        }
      }

      for(int i = 0 ; i < toBeErased.size() ; i++)
        toBeErased.erase(toBeErased.begin() + i);
    }

    /*bool isALevelToBeRewritten(int level) {
      return (find(levelsToBeRewritten.begin(), levelsToBeRewritten.end(), level) != levelsToBeRewritten.end());
    }

    vector<int>& getLevelsToBeRewritten() {
      vector<int>& ref = levelsToBeRewritten;

      return ref;
    }*/
};

// rewrite rows on critical path (CP)
// TODO: implement
class RewriteOnCP : public RewritingStrategy {
  protected:
    vector< vector<int> > paths;
    shared_mutex mutex;

  public:
    RewriteOnCP(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) : 
                      RewritingStrategy(rows, startPoint, levels, levelTable, dag) {
      policy();
    }

    void policy() {
      int pathWithMaxLength = 0;
      int pathLength = 0;
      for(int i = 0 ; i < paths.size(); i++) {
        vector<int>& path = paths[i];
        if(pathLength < path.size()) {
          pathLength = path.size();
          pathWithMaxLength = i;
        }
      }
       
      vector<int>& criticalPath = paths[pathWithMaxLength];
      numOfRewrittenRows = criticalPath.size();
      for(auto& row : criticalPath)
        toBeRewritten[row] = make_pair(levels[row],-1); 
    }

    void traversePath(vector<int>& path, int row) {
      int parentMaxLevel = row;
      int maxLevel = 0;
      if(levels[row] > 0)
        for(auto& parent : dag[row].first) {
          if(levels[parent] > maxLevel) {
            maxLevel = levels[parent];
            parentMaxLevel = parent; 
          }
        }

      path.push_back(parentMaxLevel);

      if(levels[parentMaxLevel] > 0)
        traversePath(path, parentMaxLevel);
    }

    void buildCriticalPath(vector<int>& lastLevelRows) {
      #pragma omp parallel for schedule(static)
      for(int i = 0 ; i < lastLevelRows.size() ; i++) {
        int row = lastLevelRows[i];
        vector<int> path;
        path.push_back(row);
        traversePath(path, row);

        if(path.back() != 0)
          cout << "path doesnt end in 0\n";

        mutex.lock();
        paths.push_back(path);
        mutex.unlock();
      }
    }
};

// rewrite rows with the smallest # of in-degree at selectec levels
// TODO: implement
class RewriteLowestIndegree : public RewritingStrategy {
  public:
    RewriteLowestIndegree(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag) : 
                      RewritingStrategy(rows,startPoint,levels,levelTable,dag) {
      policy();
    }

    void policy() {
    }
};


class RewriteByThreeCriteria : public RewritingStrategy {
  protected:
    vector<int> levelCost;
    map<int,int> levelSizeBelowAvg;

    double avgCostPerLevel;
    double avgNumRowsPerLevel;
    double avgIndegrePerRow;
    
    map<int,double> flopsBelowAvg;
    vector<int> numOfLevels;

    public:
      RewriteByThreeCriteria(int rows, StartPoint startPoint, int* levels, vector< vector<int> >& levelTable, DAG& dag, vector<int>& flopsPerLevel, double avgFLOPSPerLevel, map<int,double> flopsBelowAvg) :
                      RewritingStrategy(rows, startPoint, levels, levelTable, dag) {
          levelCost = flopsPerLevel;

          // get # of rows for levels with flopsBelowAvg
          for(auto& level: flopsBelowAvg)
            levelSizeBelowAvg[level.first] = levelTable[level.first].size();

          this->avgCostPerLevel = avgFLOPSPerLevel;
          this->flopsBelowAvg = flopsBelowAvg;

          cout << "Before Policy\n";
          policy();
          
          cout << "num. of rewritten rows: " << toBeRewritten.size() << "\n";
          printRowsToBeRewritten();

      }

      void policy() {
        AvgRowsPerLevel();
        AvgInPerRow();

        cout << "old level sizes:\n";
        for(auto& level : levelSizeBelowAvg)
          cout << level.first << ", " << level.second << "\n";

        Rewrite(avgNumRowsPerLevel,avgCostPerLevel,avgIndegrePerRow);
      
        cout << "new level sizes:\n";
        for(auto& level : levelSizeBelowAvg)
          cout << level.first << ", " << level.second << "\n";
      }

    int AvgInPerRow(){
      double sum=0;
      for(int i = 0; i < dag.size(); i++){
        Connectivity& connectivity = dag[i];
        vector<int>& parents = connectivity.first;
        sum += parents.size();
      }

      //out<<dag.size()<<"\n";
      avgIndegrePerRow = ceil(sum / dag.size());        
      cout<<"Total number of indegree:"<< sum <<"\n";
      cout<<"Average Indegree Per row:" << avgIndegrePerRow <<"\n";
      
      return avgIndegrePerRow;
    } 

    int AvgRowsPerLevel() { 
      double rowSum = 0;
      cout<<"Number of levels below average cost:"<< flopsBelowAvg.size()<< "\n";
      // for(auto& level : flopsBelowAvg){
      //   //if(level.second <= avgCostPerLevel)
      //   numOfLevels.push_back(level.first);
      //  }
       //cout<< "Number of levels below average cost:"<<numOfLevels.size() <<"\n";
       //cout<< levelTable[2].size();

      for(auto it = flopsBelowAvg.begin(); it !=flopsBelowAvg.end(); it++){
          //cout<< levelTable[it->first].size()<<"\n";
          rowSum += levelTable[it->first].size();
          //cout << "RowSum:" << rowSum<<"\n";      
      }  

      cout << "Total number of rows:" << rowSum << "\n";
      avgNumRowsPerLevel = ceil(rowSum / flopsBelowAvg.size());
      cout << "Average number of rows per level at levels lower than Average Cost: "<< avgNumRowsPerLevel <<"\n";

      return avgNumRowsPerLevel;
    } 

    void Rewrite(double avgNumRowsPerLevel, double avgCostPerLevel, double avgIndegrePerLevel) {
      auto it = flopsBelowAvg.begin();
      int targetLevel = it->first;
      it++;
      for(; it != flopsBelowAvg.end(); it++){
        vector<int>& levelRows = levelTable[it->first];
        cout << "level: " << it->first << " level size: " << levelRows.size() << "\n";
        cout << "target level: " << targetLevel << " with size: " << levelTable[targetLevel].size() << "\n";

        for(auto& row : levelRows){
          vector<int> predsWithMaxLevel;
          vector<int> parents = dag[row].first;

          int newLevel;
          int maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);

          cout << "\nrow: " << row << "\n";
          cout << "max level: " << maxLevel << "\n";
          cout<< "indegree : "<< parents.size() <<"\n";

          while(maxLevel >= targetLevel) {
           for(auto& pred : predsWithMaxLevel)
             expandPredsWith(pred, parents, row);

           predsWithMaxLevel.clear();
           int oldLevel = maxLevel;

           // level with size > ARL cannot be a target level
           // skip levels with cost > avgCostPerLevel (i.e. not in flopsBelowAvg)
           cout << "original maxLevel size: " << levelTable[maxLevel].size() << "\n";
           int maxLevelSize = levelTable[maxLevel].size();
           if(flopsBelowAvg.find(maxLevel) != flopsBelowAvg.end())
             maxLevelSize = levelSizeBelowAvg[maxLevel];
          
           while(maxLevelSize >= avgNumRowsPerLevel ||
   //              levelTable[maxLevel].size() == 0 ||
                 flopsBelowAvg.find(maxLevel) == flopsBelowAvg.end()) {
   //          cout << "max level: " << maxLevel << " with size: " << levelTable[maxLevel].size() << " with cost: "<< flopsBelowAvg[maxLevel] << "\n";
             maxLevel = findPredecessorsWithMaxLevel(parents, predsWithMaxLevel);
             for(auto& pred : predsWithMaxLevel)
               expandPredsWith(pred, parents, row);
             predsWithMaxLevel.clear();

             maxLevelSize = levelTable[maxLevel].size();
             if(flopsBelowAvg.find(maxLevel) != flopsBelowAvg.end())
               maxLevelSize = levelSizeBelowAvg[maxLevel];
           }

           cout << "levelcost maxlevel: " << levelCost[maxLevel] << " avgCostPerLevel: " << avgCostPerLevel << " indegrees: " << parents.size() << " AIR: " << avgIndegrePerLevel << "\n";
           // we found a candidate level
           if(levelCost[maxLevel] < avgCostPerLevel && 
              parents.size() <= avgIndegrePerLevel) {
             cout << "new maxLevel: " << maxLevel << "\n";
                // TODO: what if parents.size() > avgIndegrePerLevel
                // maybe use distance between instances
             int numOfPreds = parents.size();
             int costRow = numOfPreds <= 4 ? (numOfPreds << 1) + 1 : (numOfPreds << 1);
             cout << "indegree: " << parents.size() << " costRow: " << costRow << " oldLevel: " << oldLevel << "\n\n";
             levelCost[maxLevel] += costRow;
             levelSizeBelowAvg[maxLevel]++;
             levelSizeBelowAvg[it->first]--;
             if(levelSizeBelowAvg[it->first] == 0) {
               it++;
               targetLevel = it->first;
             }
             toBeRewritten[row] = make_pair(levels[row],maxLevel);
             break;
           }

           if(oldLevel == maxLevel) break;
         } // while
        } // for
      }
    }

    void expandPredsWith(int row, vector<int>& preds, int child) {
      auto it = lower_bound(preds.begin(),preds.end(),row);
      // Replace row with first pred. of it instead of erasing and causing a shift.
      // Then push back the rest since order doesnt matter.
      if(it != preds.end()) {
        vector<int>& parents = dag[row].first;

        if(!parents.empty()) { 
          preds[it-preds.begin()] = parents[0];
          preds.insert(preds.end(),parents.begin()+1,parents.end());
        }

        sort(preds.begin(), preds.end());
        auto it = unique(preds.begin(), preds.end());
        preds.resize(std::distance(preds.begin(),it));
      }
      else
        cout << "row " << row << " doesnt exist in predecessors\n";
    }
};
