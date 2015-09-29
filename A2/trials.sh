#!/bin/bash
#./binarize_openmp images/bee.png images/bee_pthread.png 4
#./binarize_pthreads images/bee.png images/bee_pthread.png 4
#./binarize_sequential images/bee.png images/bee_pthread.png

#echo "Sequential Executing"
#for i in {1..100}
#do
#  ./sobel_sequential images/bee.png images/bee_sobel_seq.png
#done

#echo "Pthread Executing"
#for i in 2 4 8 16 32
#do
#  echo "doing"
#  echo $i
#  echo " "
#  for j in {1..100}
#  do
#    ./sobel_pthreads images/bee.png images/bee_sobel_pthread.png $i
#  done
#done

echo "OpenMP executing"
for i in 2 4 8 16 32
do
  echo "doing"
  echo $i
  echo ""
  for j in {1..100}
  do
    ./sobel_openmp images/bee.png images/bee_sobel_omp.png $i
  done
done