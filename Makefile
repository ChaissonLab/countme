all: countme mecons example

ABPOA_DIR=abPOA-1.4.1
SPOA_DIR=/project2/mchaisso_100/mchaisso/projects/CARD/tr_methylation/spoa/inst/

countme: methyl_extract.cpp
	g++ -std=c++17 -O3 -o countme methyl_extract.cpp -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include 


msa.o: msa.cpp msa.h
	g++ -g -c msa.cpp -I$(ABPOA_DIR)/include

mecons: methyl_consensus.cpp
	g++ -O2 -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib64/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons 
#	g++ -fsanitize=address -g -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib64/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons 


example: example.c
	g++ -g -I$(ABPOA_DIR)/include example.c  -L$(ABPOA_DIR)/lib -labpoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o example
