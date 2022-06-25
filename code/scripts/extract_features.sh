
echo "matrix,failed indegree, failed row cost, failed level size" >> features.csv
./extract_feature.sh "fail" | cut -d: -f2 | sed -e '/^$/d' | paste - - - - -d,  >> features.csv

echo "matrix,# of rewritten rows"  >> features.csv
./extract_feature.sh "rewritten rows" | cut -d: -f2 | sed -e '/^$/d' | paste - - -d,  >> features.csv

echo "matrix, level size after rewrite"  >> features.csv
./extract_feature.sh "after" | cut -d: -f2 | sed -e '/^$/d' | paste - - -d, >> features.csv

