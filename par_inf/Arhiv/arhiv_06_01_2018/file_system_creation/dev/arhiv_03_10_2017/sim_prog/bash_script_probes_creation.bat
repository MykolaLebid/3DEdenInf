# set number of processors for each thread
export OMP_NUM_THREADS = 1
# send all jobs to Euler
# -n = number of cores (should be the same as OMP_NUM_THREADS)
# -W = limit time for each thread
# -R rusage[mem=XXX] - operational memory per processor
# -R rusage[scratch=XXX] - disk memory per processor
for((i_0=0;i_0<=2;i_0++))
do
  for((i_1=0;i_1<=1;i_1++))
  do
    for((i_2=0;i_2<=12;i_2++))
    do
      bsub -n 1 -W 240 -R "rusage[mem=1024]"  -R "rusage[scratch=0]" ./cancer_3d.exe ../FS3/driver_adv_$i_0/driver_advantage_$i_1/dist_type_$i_2/ init_conf_file.dat
    done
  done
done
