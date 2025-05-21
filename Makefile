all: countme mecons 
SPOA_DIR=spoa/install/usr/local


countme: methyl_extract.cpp
	g++ -std=c++17 -O3 -o countme methyl_extract.cpp -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include 

$(SPOA_DIR)/lib64/libspoa.so:
	cd spoa && mkdir build
	cd spoa/build && meson 
	cd spoa/build && ninja
	cd spoa/build && meson install --destdir=../install/

mecons: methyl_consensus.cpp $(SPOA_DIR)/lib64/libspoa.so
	g++ -O2 -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib64/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons  -Wl,-rpath-link=$(SPOA_DIR)/lib64

#	g++ -fsanitize=address -g -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib64/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons 


