# set number of processors for each thread
export OMP_NUM_THREADS = 1
# send all jobs to Euler
# -n = number of cores (should be the same as OMP_NUM_THREADS)
# -W = limit time for each thread
# -R rusage[mem=XXX] - operational memory per processor
# -R rusage[scratch=XXX] - disk memory per processor
for((i_0=0;i_0<=1;i_0++))
do
  for((i_1=0;i_1<=1;i_1++))
  do
    for((i_2=0;i_2<=1;i_2++))
    do
      for((i_3=0;i_3<=1;i_3++))
      do
        for((i_4=0;i_4<=2;i_4++))
        do
          bsub -n 1 -W 240 -R "rusage[mem=1024]"  -R "rusage[scratch=0]" cancer_3d.exe ..\test_FS1\driver_adv_$i_0\mutation_rate_$i_1\driver_mutation_rate_$i_2\seed_$i_3\tree_dist_type_$i_4\ init_conf_file.dat
        done
      done
    done
  done
done
