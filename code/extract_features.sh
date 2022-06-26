
rm failed.csv rewritten.csv levelSize.csv features.csv

echo  "matrix,failed indegree, failed row cost, failed level size" >> failed.csv
./extract_feature.sh "fail" | cut -d: -f2 | sed -e '/^$/d' | paste - - - - -d,  >> failed.csv

echo -e "\nmatrix,# of rewritten rows"  >> rewritten.csv
./extract_feature.sh "rewritten rows" | cut -d: -f2 | sed -e '/^$/d' | paste - - -d,  >> rewritten.csv

echo  -e "\nmatrix, # of levels after rewrite"  >> numOfLevels.csv
./extract_feature.sh "after" | cut -d: -f2 | sed -e '/^$/d' | paste - - -d, >> numOfLevels.csv

cat failed.csv rewritten.csv numOfLevels.csv >> features.csv


