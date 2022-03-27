#include <ios>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "rewrite.h"
#include "util.h"

using namespace std;

void Rewrite::writeFunc(std::ofstream &stream, vector< vector<int> >& tracker, int size, int part, int maxNumOfThreads, int headerCounter) {
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";

  vector<int>& partStarts   = tracker[0];
  vector<int>& threadCounts = tracker[1];

/*  cout << "size: " << size << "\n";
  cout << "tracker table:\n";
  for(int i = 0; i < size ; i++)
    cout << i << ": " << partStarts[i] << ", " << threadCounts[i] << "\n";
  cout << "\n"; */

/*  for(int i = partStarts[0]; i < partStarts[size-1]+threadCounts[size-1]; i++)
    stream << "#include \"calculate" << to_string(i) << ".h\"\n";
  stream << "\n";*/

  stream << "#include \"calculators" + to_string(headerCounter) + ".h\"\n\n";

  int sign;
  if((part+1) * TABLE_SIZE > signatureLevel.size())
    sign = accumulate(signatureLevel.begin() + part * TABLE_SIZE, signatureLevel.end(), 0);
  else
    sign = accumulate(signatureLevel.begin() + part * TABLE_SIZE, signatureLevel.begin() + ((part+1) * TABLE_SIZE), 0);
  if(sign)
    stream << "void run" << part << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";
  else
    stream << "void run" << part << "(double* x) {\n";

  int counter = partStarts[0];
  if (maxNumOfThreads > 1)
    stream << "\n#pragma omp parallel num_threads(" << maxNumOfThreads << ")\n{\n";

  int levelCounter = part * TABLE_SIZE;
  int threadCounter = 0;
//  cout << "size: " << size << "\n";
  for(int j = 0 ; j < size ; j++) {
//    vector<int>& line = tracker[j];
  
    stream << "\n";
    if(threadCounts[j] > 1) {
      stream << "  #pragma omp sections // " << j << ", " << to_string(threadCounts[j]) << "\n"
             << "  {\n";

 //     cout << "counter: " << counter << " threadCounts[j]: " << threadCounts[j] << " signature size: " << signature[levelCounter].size() << "\n";
      for(int i = counter ; i < counter+threadCounts[j] ; i++) {
        stream << "    #pragma omp section\n";
  //      cout << "levelCounter: " << levelCounter << " threadCounter: " << threadCounter << "\n";
        if(signature[levelCounter][threadCounter] == 0)
          stream << "    { calculate" << i << "(x); }\n";
        else
          stream << "    { calculate" << i << "(x, b, parents, values, rowPtr, rowIndices); }\n";

        threadCounter++;
      }
      stream << "  }\n";
    } else {
       if(signature[levelCounter][0] == 0) {
         stream << "  #pragma omp single\n"
                << "  calculate" << counter << "(x); \n";
       } else {
         stream << "  #pragma omp single\n"
                << "  calculate" << counter << "(x, b, parents, values, rowPtr, rowIndices);\n";
       }

    }

    counter += threadCounts[j];
    levelCounter++;
    threadCounter = 0;
  }

  if (maxNumOfThreads > 1)
    stream << "  }\n";

  stream << "}\n";
    
  stream.close();
}

void Rewrite::allocateMemory(std::ostream &stream, int numOfParts) {
  Part* L = matrixCSR->getL();
  int rows = L->getRows();
  int vals = L->getNNZs();

  stream << "#include <stdlib.h>\n"
         << "#include <omp.h>\n"
         << "#include <math.h>\n"
         << "#include <stdio.h>\n"
         << "#include \"util.h\"\n\n";

  stream << "#define FILEPATH \"/tmp/" + fileName + ".bin\"\n"
         << "#define X_SIZE " +  to_string(rows-1) + "\n"
         << "#define FILESIZE (X_SIZE * sizeof(double))\n";

 if(analyzer->getSingleLoopRows()) {
   stream << "#define FILEPATH_B \"/tmp/" + fileName + "_b.bin\"\n"
          << "#define FILEPATH_PARENTS \"/tmp/" + fileName + "_parents.bin\"\n"
          << "#define FILEPATH_ROWPTR \"/tmp/" + fileName + "_rowPtr.bin\"\n"
          << "#define FILEPATH_ROWINDICES \"/tmp/" + fileName + "_rowIndices.bin\"\n"
          << "#define FILEPATH_VALUES \"/tmp/" + fileName + "_vals.bin\"\n"
          << "#define FILESIZE_P (PARENTS_SIZE * sizeof(int))\n"
          << "#define FILESIZE_ROWPTR (X_SIZE * sizeof(int))\n"
          << "#define FILESIZE_ROWINDICES (" + to_string(rowIndices.size()) + " * sizeof(int))\n"
          << "#define FILESIZE_VAL (PARENTS_SIZE * sizeof(double))\n"
          << "#define PARENTS_SIZE " +  to_string(vals) + "\n\n"; 
 }

 stream << "#include <unistd.h>\n"
        << "#include <time.h>\n"
        << "\n"
        << "struct timespec t1, t2;\n"
        << "\n"
        << "#define TIME(...) fprintf(stdout, __VA_ARGS__)\n"
        << "#define START() clock_gettime(CLOCK_MONOTONIC, &t1);\n"
        << "#define STOP(event) \\\n"
        << " do {  \\\n"
        << "     double timeSpent = 0.0; \\\n"
        << "     struct timespec temp; \\\n"
        << "     clock_gettime(CLOCK_MONOTONIC, &t2); \\\n"
        << "     if (((&t2)->tv_nsec-(&t1)->tv_nsec)<0) { \\\n"
        << "         temp.tv_sec = (&t2)->tv_sec-(&t1)->tv_sec-1; \\\n"
        << "         temp.tv_nsec = 1000000000+(&t2)->tv_nsec-(&t1)->tv_nsec; \\\n"
        << "     } else { \\\n"
        << "         temp.tv_sec = (&t2)->tv_sec-(&t1)->tv_sec; \\\n"
        << "         temp.tv_nsec = (&t2)->tv_nsec-(&t1)->tv_nsec; \\\n"
        << "     } \\\n"
        << "     \\\n"
        << "     timeSpent = (temp.tv_nsec) / 1e6; \\\n"
        << "     timeSpent = timeSpent + temp.tv_sec*1000; \\\n"
        << "     TIME(\"%s: %4.2f ms\\n\",event,timeSpent); \\\n"
        << "    } while (0)\n\n";


  for(int i = 0; i < numOfParts; i++) {
    int sign;
    if((i+1) * TABLE_SIZE > signatureLevel.size())
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.end(), 0);
    else
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.begin() + ((i+1) * TABLE_SIZE), 0);

    if(sign == 0)
      stream << "void  run" << i << "(double* x);\n"; 
    else
      stream << "void run" << i << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";
  }

  stream << "\n";

  stream << "int main() {\n"
         << "  int fd;\n"                 // holds x (x_r & x_w)
         << "  int result;\n\n"
         << "  double* x;\n"
         << "  fd = open(FILEPATH, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);\n"
         << "  check_err_open_lseek(fd);\n"
         << "  result = lseek(fd, FILESIZE-1, SEEK_SET);\n"
         << "  check_err_open_lseek(result);\n"
         << "  result = write(fd, \"\", 1);\n"
         << "  check_err_write(result, fd);\n\n"
         << "  x = mmap(0, FILESIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);\n"
         << "  check_fail_mmap(x, fd);\n\n";

 if(analyzer->getSingleLoopRows()) {
   stream << "  int fd_b;\n"             // holds b vector
          << "  int fd_parents;\n"       // holds parents vector
          << "  int fd_vals;\n"          // holds values vector
          << "  int fd_rowPtr;\n"        // holds rowPtr vector
          << "  int fd_rowIndices;\n\n"  // holds rowIndices vector
          << "  double* b;\n"
          << "  int* parents;\n"
          << "  int* rowPtr;\n"
          << "  int* rowIndices;\n"
          << "  double* values;\n\n"
          << "  fd_b = open(FILEPATH_B, O_RDWR);\n"
          << "  fd_parents = open(FILEPATH_PARENTS, O_RDWR);\n"
          << "  fd_rowPtr = open(FILEPATH_ROWPTR, O_RDWR);\n"
          << "  fd_rowIndices = open(FILEPATH_ROWINDICES, O_RDWR);\n"
          << "  fd_vals = open(FILEPATH_VALUES, O_RDWR);\n\n"
          << "  check_err_open_lseek(fd_b);\n"
          << "  check_err_open_lseek(fd_parents);\n"
          << "  check_err_open_lseek(fd_rowPtr);\n"
          << "  check_err_open_lseek(fd_rowIndices);\n"
          << "  check_err_open_lseek(fd_vals);\n\n"
          << "  parents = mmap(0, FILESIZE_P, PROT_READ, MAP_SHARED, fd_parents, 0);\n"
          << "  values = mmap(0, FILESIZE_VAL, PROT_READ, MAP_SHARED, fd_vals, 0);\n"
          << "  rowPtr = mmap(0, FILESIZE_ROWPTR, PROT_READ, MAP_SHARED, fd_rowPtr, 0);\n"
          << "  rowIndices = mmap(0, FILESIZE_ROWINDICES, PROT_READ, MAP_SHARED, fd_rowIndices, 0);\n"
          << "  b = mmap(0, FILESIZE, PROT_READ, MAP_SHARED, fd_b, 0);\n\n"
          << "  check_fail_mmap(b, fd_b);\n"
          << "  check_fail_mmap((double*)rowPtr, fd_rowPtr);\n"
          << "  check_fail_mmap((double*)rowIndices, fd_rowIndices);\n"
          << "  check_fail_mmap((double*)parents, fd);\n"
          << "  check_fail_mmap(values, fd);\n\n";
 }

  stream << "START();\n";   
  for(int i = 0; i < numOfParts; i++) {
    int sign;
    if((i+1) * TABLE_SIZE > signatureLevel.size())
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.end(), 0);
    else
      sign = accumulate(signatureLevel.begin() + i * TABLE_SIZE, signatureLevel.begin() + ((i+1) * TABLE_SIZE), 0);

    if(sign == 0)
      stream << "  run" << i << "(x);\n"; 
    else
      stream << "  run" << i << "(x, b, parents, values, rowPtr, rowIndices);\n"; 
  }

  stream << "STOP(\"execution\");\n\n";   

  stream << "  int errCnt = 0;\n"
         << "  for (int i = 0; i < " << (L->getRows()-1) << " ; i++) {\n"
         << "    if(((fabs(1.0000000000 - x[i])/1.0000000000) >= 1e-2) && fabs(x[i]) >= 0) {\n"
         << "      errCnt++;\n"
         << "      printf(\"x[%d]: %.5f\\n\",i,x[i]);\n"
         << "    }\n"
         << "  }\n\n";

  stream << "  if(!errCnt)\n"
         << "    printf(\"Chainbreaker passed!\\n\");\n"
         << "  else\n"
         << "    printf(\"Chainbreaker failed! errCnt:%d\\n\",errCnt);\n\n";

  if(analyzer->getSingleLoopRows())
    stream << "  if(munmap(b, FILESIZE) == -1 || \n \
                    munmap(parents, FILESIZE_P) == -1 || munmap(values, FILESIZE_VAL) == -1 ||  \n \
                    munmap(rowPtr, FILESIZE) == -1 || munmap(rowIndices, FILESIZE) == -1)\n"
           << "    perror(\"Error un-mmapping the file\");\n\n"
           << "  close(fd_b);\n"
           << "  close(fd_parents);\n"
           << "  close(fd_rowPtr);\n"
           << "  close(fd_rowIndices);\n"
           << "  close(fd_vals);\n";

  stream << "  if(munmap(x, FILESIZE) == -1)\n"
         << "    perror(\"Error un-mmapping the file\");\n\n"
         << "  close(fd);\n"
         << "  return 0;\n"
         << "}";
}

void Rewrite::writeUtil() {
  std::ofstream stream(fileName + "/util.h");
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";

  stream << "#include <sys/types.h>\n"
         << "#include <sys/stat.h>\n"
         << "#include <unistd.h>\n"
         << "#include <fcntl.h>\n"
         << "#include <sys/mman.h>\n\n"
         << "void check_err_open_lseek(int result) {\n"
         << "  if(result == -1) {\n"
         << "    perror(\"Error opening file/lseek for reading\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n"
         << "void check_fail_mmap(double* map, int fd) {\n"
         << "  if(map == MAP_FAILED) {\n"
         << "    close(fd);\n"
         << "    perror(\"Error mmapping the file\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n"
         << "void check_err_write(int result, int fd) {\n"
         << "  if(result != 1) {\n"
         << "    close(fd);\n"
         << "    perror(\"Error writing last byte of the file\");\n"
         << "    exit(EXIT_FAILURE);\n"
         << "  }\n"
         << "}\n\n";

  stream.close();
}

#ifdef REWRITE_ENABLED
int Rewrite::rewriteRow(int row, vector<double> &b, string& currRow, int rewriteDepth) {
  int flops = 0;

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  if(parents.empty()) {
    if(rewriteDepth == REWRITE_DEPTH)
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow.append("b[" + to_string(row) + "] / " + to_string(rowValues.back()));
      currRow.append(to_string(b[row]) + " / " + to_string(rowValues.back()));
      flops++;
    }
  } else {
    currRow += "((" + to_string(b[row]) + " - (";
    flops++;
    for(int j = 0 ; j < parents.size() ; j++) {
      flops++;
      if(rewriteDepth == REWRITE_DEPTH) {
        currRow += to_string(rowValues[j]) + "*x[" + to_string(parents[j]) + "]";
      } else {
        currRow += to_string(rowValues[j]) + "*(";
        flops += rewriteRow(parents[j], b, currRow, rewriteDepth+1);
        currRow += string(")");
      }

      if(j != parents.size() - 1) {
        currRow += string(" + ");
        flops++;
      }
    }

    currRow.append(string(")) / ") + to_string(rowValues.back()) + ")");
    flops++;
  }

  //cout << "row: " << row << " added " << flops << "\n";
  return flops;
}

int Rewrite::dumpString(int row, stringstream& currRow, set<int>& rewritten, map<int,double>& multipliers) {
  int flops = (multipliers.size() << 1);

//  cout << "dumpString multipliers[-1]: " << multipliers[-1] << "\nmultipliers size: " << multipliers.size() << "\n";
  currRow << multipliers[-1];
  multipliers.erase(-1);

  for(auto it = multipliers.begin() ; it != multipliers.end(); it++)
    if(rewritten.find(it->first) == rewritten.end())
      currRow << " - " << it->second <<  " * x[" << it->first <<  "]";

  currRow << ";\n";
//  cout << "generated string: " << currRow << "\n";

  return flops;
}

// doesnt need to return flops
// put -1 as key for constant to multipliers
// map<int,double> multipliers;  // the constant will have the key -1
// rewritingMultiplicant is 1 (this is for rewritten 
void Rewrite::rewriteRow2(int row, vector<double> &b, set<int>& rewritten, map<int,double>& multipliers, double rewritingMultiplicant, char sign) {

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  //cout.fixed;
  //cout.precision(10);
//  cout << "row: " << row << " rewritingMultiplicant: " << rewritingMultiplicant <<  " sign: " << (int)sign << "\n";
  if(!parents.empty()) {
    for(int j = 0 ; j < parents.size() ; j++) {
      auto it = rewritten.find(parents[j]);
      if(it != rewritten.end()) {
//        cout << fixed << setprecision(10) << "rewritten row: " << parents[j] << " rowValues[" << j << "]: " << rowValues[j] << " rowValues.back(): " << rowValues.back() << "\n";
        rewriteRow2(parents[j], b, rewritten, multipliers, rewritingMultiplicant * (rowValues[j] / rowValues.back()), !sign);
      } else {
 //       cout << fixed << setprecision(10) << " sign: " << (int)sign << " rowValues[" << j  << "]: " << rowValues[j] << " rowValues.back(): " << rowValues.back() <<  " for x" << parents[j] << ": " << rewritingMultiplicant * (rowValues[j] /rowValues.back()) << "\n"; 
//        if(multipliers.find(parents[j]) != multipliers.end())
          if((int)sign > 0) {
 //          cout << " adding to multipliers[" << parents[j] << "]: " << multipliers[parents[j]] << "\n";
           multipliers[parents[j]] +=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
          } else {
 //          cout << " removing from multipliers[" << parents[j] << "]: " << multipliers[parents[j]] << "\n";
           multipliers[parents[j]] -=  rewritingMultiplicant * (rowValues[j] /rowValues.back());
          }
  //      else
    //      multipliers[parents[j]] = rewritingMultiplicant * (rowValues[j] /rowValues.back());
      }
    }
  
//    cout << fixed << setprecision(10) << " sign: " << (int)sign << " b[" <<  row << "] (" << b[row] << ") / " << rowValues.back()  << "  rowValues.back() to the constant slot\n\n";
 //   cout << "row: " << row << " rewritingMultiplicant: " << rewritingMultiplicant << "\n";
    if((int)sign > 0)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  } else { // parents empty means I rewritten to level 0.
//    cout << "row: " << row << " rewritingMultiplicant: " << rewritingMultiplicant << " sign: " << (int)sign << "\n";
//    cout << fixed << setprecision(10) << "setting b[" << row << "]: " << b[row] << " / rowValues.back(): " << rowValues.back() << " for x" << row << "\n\n";
    if((int)sign > 0)
      multipliers[-1] += rewritingMultiplicant * (b[row] / rowValues.back());
    else
      multipliers[-1] -= rewritingMultiplicant * (b[row] / rowValues.back());
  }

//  cout << fixed << setprecision(10) << "multipliers[-1]: " << multipliers[-1] << "\n";
}

int Rewrite::rewriteRow(int row, vector<double> &b, string& currRow, set<int>& rewritten) {
  int flops = 0;

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<double>& rowValues = values[row];

  vector<int>& parents = rewritingStrategy->isRewritten(row) ? rewritingStrategy->getInitialParentsOf(row) : dag[row].first;

  if(parents.empty()) {
    auto it = rewritten.find(row);
    if(it == rewritten.end())
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow.append("b[" + to_string(row) + "] / " + to_string(rowValues.back()));
      //currRow.append(to_string(b[row]) + " / " + to_string(rowValues.back()));
      currRow.append(to_string(b[row] / rowValues.back()));
      flops++;
    }
  } else {
    auto it = rewritten.find(row);
    if(it == rewritten.end())
      currRow.append("x[" + to_string(row) + "]");
    else {
      //currRow += "((b[" + to_string(row) + "] - (";
      currRow += "((" + to_string(b[row]) + " - (";
      flops++;
      for(int j = 0 ; j < parents.size() ; j++) {
        flops++;
        auto it = rewritten.find(parents[j]);
        if(it != rewritten.end()) {
            currRow += to_string(rowValues[j]) + "*(";
     //       cout << "calling for parent " << parents[j] << "\n";
            flops += rewriteRow(parents[j], b, currRow, rewritten);
            currRow += string(")");
        } else
            currRow += to_string(rowValues[j]) + "*x[" + to_string(parents[j]) + "]";
  
        if(j != parents.size() - 1) {
          currRow += string(" + ");
          flops++;
        }
      }
  
      currRow.append(string(")) / ") + to_string(rowValues.back()) + ")");
      flops++;
    }
  }

  return flops;
}

// rewrite the whole level
int Rewrite::rewriteLevel(int level, vector<double> &b, string& currRow, std::ostream &stream) {
  int flops = 0;
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();
  vector<vector<int>>& levelTable = analyzer->getLevelTable();

  for(auto& row : levelTable[level]) {
    set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
    //int row = level[i];
    vector<int>& parents = dag[row].first;
    vector<double>& rowValues = values[row];

    if(parents.empty()) {
      stream << "  x[" << row << "] = " << b[row] << " / " << rowValues.back() << ";\n";
      flops++;
      continue;
    }

    ToBeRewritten& toBeRewritten = rewritingStrategy->getToBeRewritten();
    int* levels = analyzer->getLevels();

    stream << "  x[" << row << "] = (" << b[row] << " - (";
    flops++;

    string currRow;
    vector<int>& initialParents = rewritingStrategy->getInitialParentsOf(row);
    for(int j = 0 ; j < initialParents.size() ; j++) {
      currRow += to_string(rowValues[j]) + "*";
      flops++;

      if(rewritingStrategy->isScopeSelective()) {
        set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
        flops += rewriteRow(initialParents[j], b, currRow, rewritten);
      } else
        flops += rewriteRow(initialParents[j], b, currRow, 0);

      if(j != initialParents.size() - 1) {
        currRow += string(" + ");
        flops++;
      }
    }

    stream << currRow;
    stream << ")) / " << rowValues.back() << ";\n";
    flops++;
   // flops += rewriteRow(row, b, currRow, targetLevel);
  }

  return flops;
}
#endif

/*int Rewrite::writeLevel(vector<int>& level) {
  vector<int>& rowPtrL =  matrixCSR->getL()->getRowPtr();
  // TODO: these will be string to stream except for level.size()
  #pragma omp parallel for schedule(static)
  for(int i = 0 ; i < level.size() ; i++) {
      // we'll denote mmaped rowValues and flatData as
      //                     values    and parents
      stream << "  double x" << row << " = 0;\n"
             << "  for(int j = " << rowPtrL[row] << "; j < " << (rowPtrL[row + 1] - 1) << " ; j++)\n"
             << "    x" << row << " += values[j] * x_r[parents[j]];\n\n"; 
      stream << "  x_w[" << row << "] = (b[" << row << "]-x" << row << ")/values[" << (rowPtrL[row + 1]-1) << "];\n\n";
  
      flops += ((rowPtrL[row + 1] - rowPtrL[row] - 1) << 1) + 1;
  }
}*/

int Rewrite::writePart(int levelNum, int rowStartIndex, int rowEndIndex, int levelPart, vector<double> &b) {
  int flops = 0;

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  vector<int>& level = levelTable[levelNum];

  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();

  vector<int>& rowPtrL =  matrixCSR->getL()->getRowPtr();

  std::ofstream stream(fileName + "/calculate" + to_string(levelPart) + ".c");
  //stream.precision(10);
  stream.precision(5);
  //stream << scientific;

  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";

  vector<int> unrolledRows;
  vector<int> loopedRows(level.begin() + rowStartIndex, level.begin() + rowEndIndex);
  #ifdef REWRITE_ENABLED
    // TODO: this is not efficient, construct rewritten rows per level from the beginning
    vector<int> rewrittenRows;
  #endif

  analyzer->separateRows(levelNum, rowStartIndex, rowEndIndex, loopedRows, unrolledRows);
 
  if(!loopedRows.empty()) {
  //if(loopedRows.size() > 1) {
    /*cout << "level: " << levelNum << " startIndex: " << rowStartIndex << " endIndex: " << rowEndIndex << "\n";
    for(auto& row : loopedRows)
      cout << row << ", ";
    cout << "\n";*/

    if(rowStartIndex == 0) {
      signatureLevel.push_back(1);
      signature.push_back(vector<int>());
    } else
      signatureLevel[levelNum] = 1;

    signature[levelNum].push_back(1);  // if we OR signature[levelNum] == signatureLevel[levelNum]

    // remove rewritten rows from loopedRows, fill in rewrittenRows
    vector<int>::iterator it = loopedRows.begin();
    while(it != loopedRows.end()) {
      #ifdef REWRITE_ENABLED
        if(rewritingStrategy->isRewritten(*it)) {
          rewrittenRows.push_back(*it);
          it = loopedRows.erase(it);
        } else {
      #endif
        flops += ((dag[*it].first.size()) << 1) + 1;
        it++;
      #ifdef REWRITE_ENABLED
       }
      #endif
    }
  } else {
    if(rowStartIndex == 0) {
      signatureLevel.push_back(0);
      signature.push_back(vector<int>());
    }
      
    signature[levelNum].push_back(0);
  }


  #ifdef REWRITE_ENABLED
    if(!unrolledRows.empty()) {
      // remove rewritten rows from unrolledRows, fill in rewrittenRows
      vector<int>::iterator it = unrolledRows.begin();
      while(it != unrolledRows.end()) {
        if(rewritingStrategy->isRewritten(*it)) {
          rewrittenRows.push_back(*it);
          it = unrolledRows.erase(it);
        } else
          it++;
      }
    }
  #endif

  rowIndices.insert(rowIndices.end(), loopedRows.begin(), loopedRows.end());
  vector<int>& currLevelStartIndex = startIndex[levelNum];
  if(currLevelStartIndex.empty()) // this will work even if loopedRows is empty since start & end will be the same for the loop
    currLevelStartIndex.push_back(startIndex[levelNum-1].back() + loopedRows.size());
  else
    currLevelStartIndex.push_back(currLevelStartIndex.back() + loopedRows.size());

  /*cout << "rowIndices:\n";
  for(auto& index : rowIndices)
    cout << index << ", ";
  cout << "\n";*/


  if(!loopedRows.empty()) {
    // TODO: we actually can remove rowIndices from loopedRows with size == 1. That'll complicate things though.
    stream << "void calculate" << levelPart << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices) {\n";

    if(loopedRows.size() > 1) {
      if(currLevelStartIndex.size() == 1)
        stream << "  for(int i = " << startIndex[levelNum-1].back() << " ; i < " << currLevelStartIndex.back() << " ; i++) {\n";
      else
        stream << "  for(int i = " << currLevelStartIndex[currLevelStartIndex.size()-2] << " ; i < " << currLevelStartIndex.back() << " ; i++) {\n";

      stream << "    int row = rowIndices[i];\n";
    } else
        stream << "    int row = " << loopedRows[0] << ";\n";
          
    stream << "    double xi = 0;\n"
           << "    for (int j = rowPtr[row]; j < rowPtr[row+1]-1; j++)\n"
           << "      xi += values[j] * x[parents[j]];\n\n"
           << "    x[row] = (b[row]-xi)/values[rowPtr[row+1]-1];\n";

    if(loopedRows.size() > 1)
      stream << "  }\n\n";
  }

  if(!unrolledRows.empty()) {
    if(loopedRows.empty())
      stream << "void calculate" << levelPart << "(double* x) {\n";

    for(auto& row : unrolledRows) {

      vector<int>& parents = dag[row].first;
      vector<double>& rowValues = values[row];

      if(parents.empty()) {
  //      stream << setprecision(10) << "  x[" << row << "] = " << b[row] / rowValues.back() << ";\n";
        stream << std::scientific << "  x[" << row << "] = " << b[row] / rowValues.back() << ";\n";
        flops++;
      } else if(parents.size() == 1) {
  //      stream << setprecision(10) << "  x[" << row << "] = (" << b[row] << "-(" << rowValues[0] << ") * x[" << parents[0] << "])/" << rowValues.back() << ";\n";
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-(" << rowValues[0] << ") * x[" << parents[0] << "])/" << rowValues.back() << ";\n";
        flops+= 3;
      } else if(parents.size() == 2) {
  /*      stream << setprecision(10) << "  x[" << row << "] = (" << b[row] << "-((" \  */
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]))/" << rowValues.back() << ";\n";
        flops+= 5;
      } else if(parents.size() == 3) {
   /*     stream << setprecision(10) << "  x[" << row << "] = (" << b[row] << "-((" \  */
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]+(" \
               << rowValues[2] << ") * x[" << parents[2] << "]))/" << rowValues.back() << ";\n";
        flops+= 7;
      } else if(parents.size() == 4) {
     /*   stream << setprecision(10) << "  x[" << row << "] = (" << b[row] << "-((" \  */
        stream << std::scientific << "  x[" << row << "] = (" << b[row] << "-((" \
               << rowValues[0] << ") * x[" << parents[0] << "]+(" \
               << rowValues[1] << ") * x[" << parents[1] << "]+(" \
               << rowValues[2] << ") * x[" << parents[2] << "]+(" \
               << rowValues[3] << ") * x[" << parents[3] << "]))/" << rowValues.back() << ";\n";
        flops+= 9;
      }
    }
  }

  #ifdef REWRITE_ENABLED
  if(!rewrittenRows.empty()) {
    if(loopedRows.empty() && unrolledRows.empty())
      stream << "void calculate" << levelPart << "(double* x) {\n";

    ToBeRewritten& toBeRewritten = rewritingStrategy->getToBeRewritten();
    int* levels = analyzer->getLevels();

    for(auto& row : rewrittenRows) {
      if(!rewritingStrategy->isRewritten(row))
        cout << "row: " << row << " doesnt appear in rewrittenStrategy's list\n"; 


      vector<int>& parents = dag[row].first;
      vector<double>& rowValues = values[row];

      // ORIGINAL IMPLEMENTATION
      /*  stream << "  x[" << row << "] = (" << b[row] << " - (";
      flops++;
      */

      stream << "  x[" << row << "] = ";

      //string currRow;
      vector<int>& initialParents = rewritingStrategy->getInitialParentsOf(row);
      set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
      map<int,double> multipliers;
      for(auto& row : rewritten)
        multipliers[row] = 0;
      multipliers[-1] = 0;
      double rewritingMultiplicant = 1;

      //for(int j = 0 ; j < initialParents.size() ; j++) {
        // ORIGINAL IMPLEMENTATION
        //currRow += to_string(rowValues[j]) + "*";

        //rewriteRow2(initialParents[j], b, rewritten, multipliers, rewritingMultiplicant, 1);
        rewriteRow2(row, b, rewritten, multipliers, rewritingMultiplicant, 1);
     // }
      
      stringstream currRow;
      flops += dumpString(row, currRow, rewritten, multipliers);
      stream << currRow.str();
     // stream << currRow;


     // ORIGINAL IMPLEMENTATION
     /* string currRow;
      vector<int>& initialParents = rewritingStrategy->getInitialParentsOf(row);
      for(int j = 0 ; j < initialParents.size() ; j++) {
        currRow += to_string(rowValues[j]) + "*";
        flops++;

        if(rewritingStrategy->isScopeSelective()) {
          set<int>& rewritten = rewritingStrategy->getRewritingMapOf(row);
          flops += rewriteRow(initialParents[j], b, currRow, rewritten);
        } else
          flops += rewriteRow(initialParents[j], b, currRow, 0);

        if(j != initialParents.size() - 1) {
          currRow += string(" + ");
          flops++;
        }
      }

      stream << currRow;
      stream << ")) / " << rowValues.back() << ";\n";
      flops++;
      */
    }
  }
  #endif

  stream << "}\n\n";
  stream.close();

  return flops;
}

int Rewrite::writeMakefile(int execBlockCnt) {
  std::ofstream stream(fileName + "/Makefile");
  if(!stream.is_open()) {
    std::cout << "Cannot open output file!\n";
    return -1;
  }

  stream << "CXX      := clang\n"
  //<< "CXXFLAGS := -O3 -ffast-math -funroll-loops -march=native -I/kuacc/apps/llvm-omp/include\n"
  << "CXXFLAGS := -O3 -flto=thin -ffast-math -march=native \n"
  << "LDFLAGS  := -fopenmp\n"
  << "SRC      := $(wildcard calculate*.c)\n"
  << "OBJECTS  := $(SRC:%.c=%.o)\n"
  << "SRC_RUN      := $(wildcard run*.c)\n"
  << "OBJECTS_RUN  := $(SRC_RUN:%.c=%.o)\n\n"
  << "%.o: %.c\n"
  << "\t$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ -c $<\n\n"
  << "all: $(OBJECTS) $(OBJECTS_RUN)\n"
  << "\t$(CXX) $(CXXFLAGS) -c main.c \n"
  << "\t$(CXX) $(CXXFLAGS) $(OBJECTS) $(OBJECTS_RUN) main.o $(LDFLAGS) -o chainbreaker\n\n";

  stream << "clean:\n"
         << "\t$(RM) $(OBJECTS) $(OBJECTS_RUN) main.o chainbreaker";

  stream << "\n\n";

  stream.close();
  return 0;
}

int Rewrite::writeHeader(int headerCounter, int partCounterStart, int partCounterEnd) {
  std::ofstream stream(fileName + "/calculators" + to_string(headerCounter) + ".h");
  if(!stream.is_open()) {
    std::cout << "Cannot open output file!\n";
    return -1;
  }

  int startingLevel = headerCounter * TABLE_SIZE;
  int counter = 0;
  stream << "#pragma once\n";
  for(int i = partCounterStart; i < partCounterEnd; i++) {
    if(signature[startingLevel][counter] == 0)
      stream << "void calculate" << to_string(i) << "(double* x);\n";
    else
      stream << "void calculate" << to_string(i) << "(double* x, double* b, int* parents, double* values, int* rowPtr, int* rowIndices);\n";

    counter++;
    if(counter == signature[startingLevel].size()) {
      startingLevel++;
      counter = 0;
    }
  }

  stream.close();
  return 0;
}

int Rewrite::balanceLevel(int toBeBalanced, int& workloadPerThread) {
    int numThreads = NUM_THREADS;
    workloadPerThread = toBeBalanced+NUM_THREADS-1;
    workloadPerThread /= NUM_THREADS;

 //     cout << "toBeBalanced size: " << toBeBalanced << "\n";
    if(workloadPerThread < LOWER_BOUND) {
      numThreads = toBeBalanced+LOWER_BOUND-1;
      numThreads /= LOWER_BOUND;
      workloadPerThread = toBeBalanced+numThreads-1;
      workloadPerThread /= numThreads;
 //       cout << "updated numThreads: " << numThreads << "\n";
 //       cout << "updated workloadPerThread: " << workloadPerThread << "\n";
    } else
      numThreads = (const int)ceil((double)(toBeBalanced/workloadPerThread));

    return numThreads;
}

#ifdef BALANCE_FLOPS
int Rewrite::balanceRows(vector<int>& level) {
//    int flopsPerThread = (const int)ceil((double)flopsPerLevel[j]/NUM_THREADS);
  
  // rowDistBelowAvg and rowdistAboveAvg keep row distribution to threads. rowDistBelowAvg
  // cooporate with flopSum which keeps track of sum of flops in rows assigned to a thread.
  // [GOAL]: distribute rows to threads where sum of flops assigned to each thread
  // is equally likely. A row can be assigned to 1 thread at most, multiple rows can
  // be assigned to a thread. We couldn't use a map instead of these 3 vectors since
  // sum of flops assigned to a thread does not need to satisfy uniqueness propoerty.
  // As needed new slots can be opened until we hit NUM_THREADS. Upon hitting, if there's
  // still a row unassigned, the one sitting at the end of rowDistBelowAvg is used : this
  // could be improved but overhead will become too costly if we keep rowDistBelowAvg and flopSum
  // sorted and try selecting the best candidate row to add on top.
  //

  // TODO: optimize for levels with all rows with same length: use BALANCE_ROWS for these
  // TODO: refactor this section

  cout << "\nnum of rows in level " << j << " : " << level.size() << "\n";
  cout << "flopsPerLevel: " << flopsPerLevel[j] << "\n";

  int avgFlopsPerLevel = flopsPerLevel[j]/level.size();
  cout << "avgFlopsPerLevel: " << avgFlopsPerLevel << "\n";

  vector<vector<int>> rowDistBelowAvg;
  vector<vector<int>> rowDistAboveAvg;
  vector<int> flopSum(NUM_THREADS,0);

  DAG& dag = analyzer->getDAG();

  for(int k = 0 ; k < level.size(); k++) {
    //int rowFlops = ((rowPtrL[level[k] + 1] - rowPtrL[level[k]]) << 1) - 1;
    vector<int>& parents = dag[level[k]].first;
    int rowFlops = ((parents.size()+1) << 1) - 1;
    cout << "row: " << level[k] << " with flops: " << rowFlops << "\n";
    if((rowDistAboveAvg.size() + rowDistBelowAvg.size()) > NUM_THREADS) {
      /*if(rowDistBelowAvg.empty()) {
        if(rowFlops <= avgFlopsPerLevel) {
        }
      } else*/
      // find the min
      int minIndex = 0;
      for(int m = 0 ; m < flopSum.size(); m++)
        if(flopSum[m] < flopSum[minIndex])
          minIndex = m;

      if(rowFlops <= avgFlopsPerLevel) {
        vector<int>& toBeMerged = rowDistBelowAvg[minIndex];
         toBeMerged.push_back(level[k]);
        if(rowFlops + flopSum[minIndex] > avgFlopsPerLevel) {
          cout << "rowDistBelowAvg size: " << rowDistBelowAvg.size() << "\n";
          rowDistAboveAvg.push_back(toBeMerged);
          if(rowDistBelowAvg.size() == 1) {
            rowDistBelowAvg.erase(rowDistBelowAvg.end() - 1);
            flopSum.erase(flopSum.end() - 1);
          }  else {
            rowDistBelowAvg.erase(rowDistBelowAvg.begin() + minIndex);
            flopSum.erase(flopSum.begin() + minIndex);
          }
        }
      } else {
        // find the second min
        int sndMinIndex = 0;
        for(int m = 0 ; m < flopSum.size(); m++)
          if(flopSum[m] < flopSum[sndMinIndex] && m != minIndex)
            sndMinIndex = m;

        vector<int>& toBeMerged = rowDistBelowAvg[minIndex];
        vector<int>& toBeMerged2 = rowDistBelowAvg[sndMinIndex];
        toBeMerged.insert(toBeMerged.end(), toBeMerged2.begin(), toBeMerged2.end());

        if(flopSum[minIndex] + flopSum[sndMinIndex] < avgFlopsPerLevel)
          flopSum[minIndex] += flopSum[sndMinIndex];
        else {   
          rowDistAboveAvg.push_back(toBeMerged);

          flopSum.erase(flopSum.begin() + minIndex);
          rowDistBelowAvg.erase(rowDistBelowAvg.begin() + minIndex);
          if(minIndex < sndMinIndex)
            sndMinIndex--;
        }
          
        flopSum.erase(flopSum.begin() + sndMinIndex);
        rowDistBelowAvg.erase(rowDistBelowAvg.begin() + sndMinIndex);
      }
    } else {
      if(rowFlops <= avgFlopsPerLevel) {
        if(rowDistBelowAvg.empty()) {
          rowDistBelowAvg.push_back(vector<int>(1,level[k])); 
          flopSum[0] = rowFlops;
        } else {
          int i = 0;
          while((i < rowDistBelowAvg.size()) && (flopSum[i] + rowFlops > avgFlopsPerLevel))
            i++;
       
          // we couldn't find a suitable candidate. Use the last slot in rowDistBelowAvg
          // and move it to rowDistAboveAvg. not the best approach but keeps from increasing
          // the overhead
          if(i != rowDistBelowAvg.size()) { 
            flopSum[i] += rowFlops;
            vector<int>& newSlot = rowDistBelowAvg[i];
            newSlot.push_back(level[k]);
          } else {
            if(rowDistBelowAvg.empty())
              rowDistAboveAvg.push_back(vector<int>()); 
            else
              rowDistAboveAvg.push_back(rowDistBelowAvg.back()); 

            vector<int>& newSlot = rowDistAboveAvg.back();
            newSlot.push_back(level[k]);
     //       cout << "rowDistAboveAvg has " << newSlot[0] << "\n";

            if(!rowDistBelowAvg.empty()) {
              rowDistBelowAvg.erase(rowDistBelowAvg.end()-1);
              flopSum.erase(flopSum.end()-1);
            }
          }
        }
      } else {
         rowDistAboveAvg.push_back(vector<int>()); 
         vector<int>& newSlot = rowDistAboveAvg.back();
         newSlot.push_back(level[k]);
      }
    }
  }

  cout << "numThreads after BALANCE_FLOPS: " << rowDistAboveAvg.size() << ", " << rowDistBelowAvg.size() << "\n";

  for(int m = 0 ; m < rowDistBelowAvg.size() ; m++) {
    cout << m << ":\n";
    for(auto& slot : rowDistBelowAvg[m])
      cout << slot << ", ";
    cout << "\n";
  }

  for(int m = 0 ; m < rowDistAboveAvg.size() ; m++) {
    cout << rowDistBelowAvg.size() + m << ":\n";
    for(auto& slot : rowDistAboveAvg[m])
      cout << slot << ", ";
    cout << "\n";
  }

  flopSum.shrink_to_fit();
  numThreads = rowDistAboveAvg.size() + rowDistBelowAvg.size();
  vector<vector<int>> rowDist(rowDistBelowAvg);
  rowDist.insert(rowDist.end(), rowDistAboveAvg.begin(),rowDistAboveAvg.end());
}
#endif

#ifdef BALANCE_ROWS_CV
// calculate coeff. of variation(CV) = stdev/mean
void Rewrite::collectFLOPS(vector<int>& level, vector<int>& flopsLevel) {
  DAG& dag = analyzer->getDAG();

 for(auto& row : level)
  flopsLevel.push_back(((dag[row].first.size()+1) << 1) - 1);

 /*for(auto& row : level)
   cout << row << ", ";
 cout << "\n";*/
}

double Rewrite::calculateCV(int numOfSum, int start, vector<int>& level) {
  /*for(auto& i : level)
    cout << i << ", ";
  cout << "\n";*/

  double mean = accumulate(level.begin(),level.end(),0.0)/(double)numOfSum;

  double sum = 0;
  for(int i = start; i < numOfSum; i++ )
    sum += (level[i] - mean) * (level[i] - mean);

  sum /= (numOfSum-1);

 double CV = sqrt(sum)/mean;
// cout << "start: " << start << " numOfSum: " << numOfSum << "\n";
// cout << "CV: " << CV << "\n";

 return CV;
}

void Rewrite::advanceTillBinLimit(vector<int>& level, int start, int& end, int binLimit) {
  int sum = accumulate(level.begin()+start, level.begin()+end, 0);

  while(end < level.size() && sum < binLimit) {
    sum += level[end];
    end = end+1;
  }
}

void Rewrite::binLastElement(vector<int>& level, int& start, int& end, vector<int>& rowDist, vector<int>& rowDistSum) {
  // if we reached the end consider whether to add the last element to the bin or not
  double CV_start_end_1 = 0;
  int currSum = 0;

  if(start == level.size()-1) {
    rowDist.push_back(1);
    rowDistSum.push_back(level.back());
    start++;

   /* cout << "rowDist:\n";
    for(auto& i : rowDist)
      cout << i << ", ";
    cout << "\n";

    cout << "rowDistSum:\n";
    for(auto& i : rowDistSum)
      cout << i << ", ";
    cout << "\n";*/

    return;
  }

  if(end == start + 1)
    currSum = accumulate(level.begin()+start, level.begin()+end, 0);
  else
    currSum = accumulate(level.begin()+start, level.begin()+(end-1), 0);

//  cout << "checking binLastElement:\n";
//  cout << "start: " << start << " end: " << end << " currSum: " << currSum << "\n";
  vector<int> currRowDistSum(rowDistSum);
  currRowDistSum.push_back(currSum);

/*  cout << "currRowDistSum:\n";
  for(auto& i : currRowDistSum)
    cout << i << ", ";
  cout << "\n";

  cout << "end: " << end << " level.size(): " << level.size() << "\n";*/

 // if(end == level.size())
 //   currRowDistSum.push_back(level[end-1]);

  if(currRowDistSum.size() > 1)
    CV_start_end_1 = calculateCV(currRowDistSum.size(), 0, currRowDistSum);

//  cout << "CV_start_end_1: " << CV_start_end_1 << "\n";
 
  if(end == start + 1)
    end++;
    
  currSum = accumulate(level.begin()+start, level.begin()+end, 0);

// currSum += *(level.end());
// cout << "currSum: " << currSum << " *level.end(): " << *(level.end()) << "\n";
  
  vector<int> currRowDistSumEnd(rowDistSum);
  currRowDistSumEnd.push_back(currSum);

  double CV_start_end = 0;
  if(currRowDistSum.size() > 1)
    CV_start_end = calculateCV(currRowDistSumEnd.size(), 0, currRowDistSumEnd);

/*  cout << "currSum: " << currSum << "\n";  

  cout << "currRowDistSumEnd:\n";
  for(auto& i : currRowDistSumEnd)
    cout << i << ", ";
  cout << "\n";

  cout << "CV_start_end: " << CV_start_end << "\n";*/
  
  if(CV_start_end_1 < CV_start_end) {
    rowDist.push_back(end-1-start);    
    rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+(end-1), 0));    
    start = end-1;
  } else {
    rowDist.push_back(end-start);    
    rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+end, 0));    
    start = end;
    end = end+1;
  }

/*  cout << "rowDist:\n";
  for(auto& i : rowDist)
    cout << i << ", ";
  cout << "\n";

  cout << "rowDistSum:\n";
  for(auto& i : rowDistSum)
    cout << i << ", ";
  cout << "\n";*/
}

double Rewrite::balanceRowsCV(vector<int>& level, vector<int>& rowDist, vector<int>& rowDistSum, double totalCV) {
//  double totalCV = calculateCV(level.size(), 0, level);
  //int binLimit = (*max_element(begin(level), end(level))) << 1;
  int binLimit = (*max_element(begin(level), end(level)));
  int start = 0; int end = 1; 

/*  cout << "balanceRowsCV: " << totalCV << ", ";
  cout << "CV before: " << totalCV << "\n";
  cout << "binLimit: " << binLimit << "\n";*/

  if(totalCV > 0.05) {
    // construct the first bin. This is not ideal but we need an initial bin and dont know the rest of the values yet
    // hence we fill the first bin no matter what
    advanceTillBinLimit(level, start, end, binLimit);
//    cout << "new end: " << end << "\n";

    if(end == start + 1){
      rowDist.push_back(1);    
      rowDistSum.push_back(level[start]);
      start = end;
      end++;
    } else {
      rowDist.push_back(end-1-start);    
      rowDistSum.push_back(accumulate(level.begin()+start, level.begin()+(end-1), 0));    
      start = end-1;
    }

/*    cout << "new start: " << start << " new end: " << end << "\n";

    cout << "rowDist:\n";
    for(auto& i : rowDist)
      cout << i << ", ";
    cout << "\n";

    cout << "rowDistSum:\n";
    for(auto& i : rowDistSum)
      cout << i << ", ";
    cout << "\n";*/

    //while(end <= level.size()) {
    while(start < level.size()) {
      advanceTillBinLimit(level, start, end, binLimit);
//      cout << "new end: " << end << "\n";
      binLastElement(level, start, end, rowDist, rowDistSum);
    }
  } else {
    // if variation is already < 0.05 do not bother running the algorithm
    rowDist.push_back(level.size());    
    rowDistSum.push_back(accumulate(level.begin(), level.end(), 0));    
  }

  double totalCVAfter = totalCV;
  if(rowDistSum.size() > 1)
    totalCVAfter = calculateCV(rowDistSum.size(), 0, rowDistSum);

  // if we increased the CV after running the algorithm, revert back
  // this can happen due to non-ideal condition for constructing the
  // first bin
  if((totalCVAfter > totalCV) && (level.size() <= NUM_THREADS)) {
    cout << totalCV << "\n";
    return totalCV;
  }

  cout << totalCVAfter << "\n";
  return totalCVAfter;
}

void Rewrite::reduceNumOfRows(vector<int>& level, vector<int>& rowDistReduced, vector<int>& rowDistSumReduced) {
  int partitionSize = (int)(ceil((double)level.size() / NUM_THREADS));
  int i = 0;

/*  cout << "reduceNumOfRows: levels \n";
  for(auto& i : level)
    cout << i << ", ";
  cout << "\n";

  cout << "partitionSize: " << partitionSize << "\n";*/
  while((i + partitionSize) < level.size()) {
    int sum = accumulate(level.begin()+i, level.begin()+i+partitionSize, 0);
    
    rowDistSumReduced.push_back(sum);
    rowDistReduced.push_back(partitionSize);
    i += partitionSize;
  }

  if(level.size()-i >= 1) {
    int sum = accumulate(level.begin()+i, level.end(), 0);
    rowDistSumReduced.push_back(sum);
    rowDistReduced.push_back(level.size()-i);
  }

/*  cout << "rowDistReduced:\n";
  for(auto& i : rowDistReduced)
    cout << i << ", ";
  cout << "\n";

  cout << "rowDistSumReduced:\n";
  for(auto& i : rowDistSumReduced)
    cout << i << ", ";
  cout << "\n";*/
}
#endif

int Rewrite::rewriteExecutor(vector<double> &b, vector<double> &x) {
  Part *L = matrixCSR->getL();
  writeUtil();

  vector<vector<int>>& levelTable = analyzer->getLevelTable();
  auto& flopsPerLevel = analyzer->getFlopsPerLevel();

  for(int i = 0 ; i < flopsPerLevel.size(); i++)
    startIndex.push_back(vector<int>());
  startIndex[0].push_back(0);

  #ifdef REWRITE_ENABLED
    auto& flopsPerLevelRewrite = analyzer->getFlopsPerLevelRewrite();
  #endif

  int rowStartIndex = 0;
  int partCounter = 0;
  int headerCounter = 0;
  int partCounterStart = 0;

  int maxNumOfThreads = 1;

  // Tracker keeps start of each part & num of threads
  // part: workload of a thread
  vector< vector<int> > tracker;
  for(int i = 0 ; i < 2; i++)
    tracker.push_back(vector<int>(TABLE_SIZE,0));

//  #pragma omp parallel for private(rowStartIndex) schedule(static)

  // assign only 1 thread to level 0. Thread creating overhead is too much for such
  // easy & fast computation
  int flops = 0;
 // cout << "level: " << 0 << "\n";
    #ifdef BALANCE_ROWS
      int flopsBefore = flops;
      flops += writePart(0, 0, levelTable[0].size(), partCounter++, b);
      cout << flops-flopsBefore << ",";
    #elif defined BALANCE_ROWS_CV
      vector<vector<int>> partitionList;
//      readWorkloadPartition(partitionList);  
//      cout << "Reading workload partition\n";
      flops += writePart(0, 0, levelTable[0].size(), partCounter++, b);
 /*   cout << "Reading workload partition\n";
      for(auto& level : partitionList) {
        cout << "count: " << level.size() << "\n";
        for(auto& workload : level)
          cout << workload << ", ";
      cout << "\n";
      }
      cout << "\n";  */
    #else
      //vector<int>& rows = rowDist[i];
      //flops += writePart(rows, 0, rows.size(), partCounter++, b);
    #endif

  vector<int>& partStarts =   tracker[0];
  vector<int>& threadCounts = tracker[1];
  partStarts[0] = partCounter-1;
  threadCounts[0] = 1;

  int size = (0 == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE; 

  if((0 == (TABLE_SIZE-1)) || (0 == flopsPerLevel.size()-1)) {
 //   cout << "j: " << 0 << " part: " << 0 << "\n"; 
    std::ofstream stream(fileName + "/run0.c");
    writeHeader(0, 0, 1);
    writeFunc(stream, tracker, size, 0, 1, headerCounter);
  }

//  cout << "CV before & after\n";
  for(int j = 1 ; j < flopsPerLevel.size(); j++) {
    int numThreads = NUM_THREADS;
    int workloadPerThread;
    vector<int>& level = levelTable[j];
    vector<int> flopsLevel;
    int flops = 0;
    #ifdef BALANCE_ROWS
      numThreads = balanceLevel(level.size(), workloadPerThread);
    #elif defined BALANCE_ROWS_CV
      collectFLOPS(level, flopsLevel);
      vector<int> flopsLevelDist;
      //cout << "level: " << j << " num of rows: " << flopsLevel.size() << "\n";
      //cout << "num of rows: " << flopsLevel.size() << "\n";
      //cout << j << ",";
      if(flopsLevel.size() > NUM_THREADS) {
        vector<int> rowDistReduced; vector<int> rowDistSumReduced;
        vector<int> rowDist;        vector<int> rowDistSum;

        reduceNumOfRows(flopsLevel, rowDistReduced, rowDistSumReduced);
        double totalCV = calculateCV(rowDistReduced.size(), 0, rowDistSumReduced);
        double totalCVAfter = balanceRowsCV(rowDistSumReduced, rowDist, rowDistSum, totalCV);
//        cout << "totalCVAfter: " << totalCVAfter << "\n";

        // we reduced number of rows to NUM_THREADS initially
        // output of balanceLevel contains count of rows assigned to each thread &
        // sum of FLOPS assigned to each thread
        // hence apply output of balanceLevel to reduced data structures
        if(totalCV == totalCVAfter)
          // algorithm failed, revert back
          flopsLevelDist = rowDistReduced;
        else {
          int numOfElements = 0; int idx = 0;
          for(int i = 0; i < rowDist.size(); i++) {
            numOfElements = rowDist[i];
            flopsLevelDist.push_back(accumulate(rowDistReduced.begin()+idx, rowDistReduced.begin()+idx+numOfElements, 0));
            idx += numOfElements;
          }
        }
      } else {
        // number of rows < NUM_THREADS
        vector<int> rowDist;        vector<int> rowDistSum;
        double totalCV = calculateCV(flopsLevel.size(), 0, flopsLevel);
        double totalCVAfter = balanceRowsCV(flopsLevel, rowDist, rowDistSum, totalCV);
 //       cout << "totalCVAfter: " << totalCVAfter << "\n";
        if(totalCV == totalCVAfter) {
          // algorithm failed, revert back
          vector<int> byRow(flopsLevel.size(),1);
          flopsLevelDist = byRow;
        } else
          flopsLevelDist = rowDist;
      }
        
      numThreads = flopsLevelDist.size();
  /*    cout << "num threads set: " << numThreads << "\n" << "flopsLevelDist:\n";
      for(auto& i : flopsLevelDist)
        cout << i << ", ";
      cout << "\n";*/

      //vector<int>& flopsLevelDist = partitionList[j];
      //numThreads = flopsLevelDist.size();
    #endif
    if(numThreads > maxNumOfThreads)
      maxNumOfThreads = numThreads;

    int workloadCounter = 0;
//    cout << "level: " << j << "\n";
//    cout << "num threads: " << numThreads << "\n";
    for(int i = 0; i < numThreads-1; i++) {
      #ifdef BALANCE_ROWS
        int flopsBefore = flops;
        flops += writePart(j, i*workloadPerThread, i*workloadPerThread+workloadPerThread, partCounter++, b);
 //       cout << flops-flopsBefore << ",";
      #elif defined BALANCE_ROWS_CV
 /*       cout << "writePart\n";
      for(auto& i : flopsLevelDist)
        cout << i << ", ";
      cout << "\n";*/
        flops += writePart(j, workloadCounter, workloadCounter+flopsLevelDist[i], partCounter++, b);
//        cout << workloadCounter << ", " << workloadCounter+flopsLevelDist[i] << " :: ";
        workloadCounter += flopsLevelDist[i];
//        cout << "workload counter: " << workloadCounter << "\n";
      #else
//        vector<int>& rows = rowDist[i];
//        flops += writePart(rows, 0, rows.size(), partCounter++, b);
      #endif
    }
    
    #ifdef BALANCE_ROWS
      int flopsBefore = flops;
      flops += writePart(j, (numThreads-1)*workloadPerThread, level.size(), partCounter++, b);
 //     cout << flops-flopsBefore << "\n";
    #elif defined BALANCE_ROWS_CV
 //     cout << "writePart last\n";
      flops += writePart(j, workloadCounter, workloadCounter+flopsLevelDist.back(), partCounter++, b);
 /*     cout << workloadCounter << ", " << workloadCounter+flopsLevelDist.back();
      cout << "\n";*/
    #else
//      vector<int>& rows = rowDist.back();
//      flops += writePart(rows, 0, rows.size(), partCounter++, b);
    #endif

    #ifdef REWRITE_ENABLED
      flopsPerLevelRewrite[j] = flops;
    #endif

    vector<int>& partStarts =   tracker[0];
    vector<int>& threadCounts = tracker[1];
    partStarts[(j % TABLE_SIZE)] = partCounter-numThreads;
    threadCounts[(j % TABLE_SIZE)] = numThreads;

    int size = (j == flopsPerLevel.size()-1) ? (flopsPerLevel.size() % TABLE_SIZE) : TABLE_SIZE; 

    if(((j % TABLE_SIZE) == (TABLE_SIZE-1)) || (j == flopsPerLevel.size()-1)) {
      int part = floor((((double)j)/TABLE_SIZE));
   //   cout << "j: " << j << " part: " << part << "\n"; 
      std::ofstream stream(fileName + "/run" + to_string(part) + ".c");
      writeHeader(headerCounter, partCounterStart, partCounter);
      writeFunc(stream, tracker, size, part, maxNumOfThreads, headerCounter);
      headerCounter++;
      partCounterStart = partCounter;
    }
  }


  std::ofstream stream(fileName + "/main.c");
  if(!stream.is_open())
    std::cout << "Cannot open output file!\n";
  //cout << "input to main.c: " << (int)(ceil((((double)flopsPerLevel.size())/TABLE_SIZE))) << "\n";
  allocateMemory(stream, (int)(ceil((((double)flopsPerLevel.size())/TABLE_SIZE))));
    
  writeMakefile(analyzer->getNumOfLevels());

  return 0;
}

// TODO: update code according ALB updates (vector -> int*, validation check code, etc)
int Rewrite::rewrite() {
  int rows = matrixCSR->getNumOfRows();
  int cols = matrixCSR->getNumOfCols();
  int nnzs = matrixCSR->getNumOfVals();

  int *rowPtr = matrixCSR->getRowPtr();
  int *colIdx = matrixCSR->getColIdx();
  double *vals = matrixCSR->getVals();
  struct timeval t1{}, t2{};
  
  if(rows != cols) {
    printf("This is not a square matrix, return.\n");
    return -1;
  }

  vector<double> xRef(cols, 1);
  vector<double> b(cols, 0);
  
  DAG& dag = analyzer->getDAG();
  vector<vector<double>>& values = analyzer->getValues();

  Part* LCSR = matrixCSR->getL();

  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////
  cout << "NUM_THREADS:" << NUM_THREADS << "\n";
  cout << "LOWER_BOUND: " << LOWER_BOUND << "\n";
  #ifdef BALANCE_ROWS
    cout << "BALANCE METHOD: BALANCE_ROWS\n";
  #else
    cout << "BALANCE METHOD: BALANCE_ROWS_CV\n";
  #endif
  cout << "TABLE_SIZE: " << TABLE_SIZE << "\n";
  /////////////////////////////// PRINT SOME BASIC STATS ///////////////////////////////


  // run spmv to get b
  gettimeofday(&t1, nullptr);
  
  std::fill(b.begin(), b.end(), 0);
  for(int i = 0; i < rows; i++) {
    #ifdef REWRITE_ENABLED
     vector<int>& parents = rewritingStrategy->isRewritten(i) ? rewritingStrategy->getInitialParentsOf(i) : dag[i].first;
    #else
     vector<int>& parents = dag[i].first;
    #endif

    vector<double>& rowValues = values[i];

    if(rowValues.empty())
      cout << "empty values for row " << i << "\n"; 

    b[i] += rowValues.back() * xRef[i];
    for(int j = 0 ; j < parents.size() ; j++)
      b[i] += rowValues[j] * xRef[parents[j]];
  }

  if(analyzer->getSingleLoopRows()) {
    string pathB ("/tmp/" + fileName + "_b.bin");
//    cout << "B array:\n";
    dumpToMem(b, pathB.c_str(), LCSR->getRows()-1);
  
    string pathVals ("/tmp/" + fileName + "_vals.bin");
//    cout << "vals array:\n";
    dumpToMem(LCSR->getVals(), pathVals.c_str(), LCSR->getNNZs());
  
    string pathParents ("/tmp/" + fileName + "_parents.bin");
    string pathRowPtr ("/tmp/" + fileName + "_rowPtr.bin");
//    cout << "cols & rows array:\n";
    flattenDumpToMem(LCSR, pathParents.c_str(), pathRowPtr.c_str());

    /*string pathParents ("/tmp/" + fileName + "_parents.bin");
    flattenDumpToMem(LCSR, pathParents.c_str());
    //dumpToMem(LCSR->getColIdx(), pathParents.c_str(), LCSR->getNNZs());

    string pathRowPtr ("/tmp/" + fileName + "_rowPtr.bin");
    cout << "rowPtr array:\n";
    dumpToMem(LCSR->getRowPtr(), pathRowPtr.c_str(), LCSR->getRows()-1);*/
  }

  gettimeofday(&t2, nullptr);
  double timeSpMV = (double) (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

  vector<double> x(cols, 0);

  gettimeofday(&t1, nullptr);
  rewriteExecutor(b, x);

  gettimeofday(&t2, nullptr);
  double chainBreaker =
      (double) (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
  printf("Chainbreaker used %4.2f ms, throughput is %4.2f gflops\n",
         chainBreaker, 2 * LCSR->getNNZs() / (1e6 * chainBreaker));

  // print SpMV throughput
  printf("\nAs a reference, L's SpMV used %4.2f ms, throughput is %4.2f "
         "gflops\n",timeSpMV, 2 * LCSR->getNNZs() / (1e6 * timeSpMV));

  if(analyzer->getSingleLoopRows()) {
    string pathRowIndices ("/tmp/" + fileName + "_rowIndices.bin");
    cout << "rowIndices array:\n";
    dumpToMem(rowIndices, pathRowIndices.c_str(), rowIndices.size());
    cout << "startIndex:\n";
    for(auto& level : startIndex) {
      for(auto& part : level)
        cout << part << ", ";
      cout << "\n";
    }
  }

  return 0;
}
