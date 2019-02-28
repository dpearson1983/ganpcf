OS = $(shell lsb_release -si)
FC = gfortran
ifeq ($(OS),Fedora)
	CXX = cuda-g++
else
	CXX = g++
endif
ARCHS = -gencode arch=compute_50,code=sm_50 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75
VXX = nvcc $(ARCHS) -ccbin=$(CXX)
FCFLAGS = -march=native -mtune=native -O3 -cpp
CXXFLAGS = -march=native -mtune=native -O3 -fPIC
VXXFLAGS = -Xptxas -dlcm=ca -lineinfo --compiler-options "$(CXXFLAGS)" -O3
LDFLAGS = -lstdc++ -lcudart -lcuda

emulator: library obj/emulator.o
	$(FC) $(LDFLAGS) -lganpcf obj/emulator.o -o emulator
	
library: obj/ganpcf_mod.o obj/ganpcf_capi.o obj/ganpcf.o
	$(FC) $(LDFLAGS) $^ -fPIC -shared -o libganpcf.so
	
obj/emulator.o: source/emulator.f90
	$(FC) $(FCFLAGS) -c source/emulator.f90 -o obj/emulator.o
	
obj/ganpcf_mod.o: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -fPIC -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.o
	
obj/ganpcf_capi.o: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) -c source/ganpcf_capi.cpp -o obj/ganpcf_capi.o
	
obj/ganpcf.o: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.o

clean:
	rm obj/*.o
	rm libnpcf.mod
	rm libganpcf.so
	rm emulator
