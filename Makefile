countme: methyl_extract.cpp
	g++ -std=c++17 -O2 -o countme methyl_extract.cpp -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include 

