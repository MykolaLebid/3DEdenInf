g++ simulation.cpp main.cpp functions.cpp sample_processing.cpp treedist.cpp phylip.cpp cons.cpp -lboost_system -std=c++11 -w -O3 -o cancer_3d.exe
chmod +x cancer_3d.exe
export OMP_NUM_THREADS=24
bsub -n 24 -W 120:00 -R "rusage[mem=30240]" -R "rusage[scratch=20240]" ./script_run_multi_times.bat
