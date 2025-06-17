all: countme mecons 
SPOA_DIR=spoa/install/usr/local


countme: methyl_extract.cpp
	g++ -std=c++17 -g -o countme methyl_extract.cpp -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include 

$(SPOA_DIR)/lib/libspoa.so:
	cd spoa && mkdir -p build
	cd spoa/build && meson
	cd spoa/build && meson setup --libdir=lib
	cd spoa/build && ninja
	cd spoa/build && meson install --destdir=../install/

mecons: methyl_consensus.cpp $(SPOA_DIR)/lib/libspoa.so
	g++ -O2 -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons  -Wl,-rpath=$(PWD)/$(SPOA_DIR)/lib


reextract: region_extract.cpp 
	g++ -O2 region_extract.cpp -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o reextract  

#	g++ -fsanitize=address -g -I$(SPOA_DIR)/include  methyl_consensus.cpp -L$(SPOA_DIR)/lib/ -lspoa -L $(CONDA_PREFIX)/lib -lhts -I $(CONDA_PREFIX)/include -lz  -o mecons 


