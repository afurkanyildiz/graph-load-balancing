
rm failed.csv rewritten.csv numOfLevels.csv features.csv

echo  "matrix,failed indegree, failed row cost, failed level size" >> failed.csv
./extract_feature.sh "fail" | cut -d: -f2 | sed -e '/^$/d' | paste - - - - -d,  >> failed.csv

echo -e "\nmatrix,# of rewritten rows"  >> rewritten.csv
./extract_feature.sh "rewritten rows" | cut -d: -f2 | sed -e '/^$/d' | paste - - -d,  >> rewritten.csv

echo  -e "\nmatrix, # of levels after rewrite"  >> numOfLevels.csv
matrix=`./extract_feature.sh "num. of levels"  | head -1 | cut -d: -f2`
levels=`./extract_feature.sh "num. of levels"  |  sed -e '/^$/d' | tail -1 | cut -d, -f2`
echo -e "\n$matrix,$levels" >> numOfLevels.csv

cat failed.csv rewritten.csv numOfLevels.csv >> features.csv


