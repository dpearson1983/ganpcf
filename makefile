FC = gfortran
CXX = g++
ARCHS = -gencode arch=compute_50,code=sm_50 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75
VXX = nvcc $(ARCHS) -ccbin=cuda-g++
FCFLAGS = -march=native -mtune=native -O3
CXXFLAGS = -march=native -mtune=native -O3
VXXFLAGS = -Xptxas -dlcm=ca -lineinfo --compiler-options "$(CXXFLAGS)" -O3
LDFLAGS = -lstdc++ -lgfortran -lcuda -lcudart

emulator: obj/ganpcf_mod.o obj/ganpcf_capi.o obj/ganpcf.o obj/emulator.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o emulator
	
obj/emulator.o: source/emulator.f90
	$(FC) $(FCFLAGS) -c source/emulator.f90 -o obj/emulator.o
	
obj/ganpcf_mod.o: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.o
	
obj/ganpcf_capi.o: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) -c source/ganpcf_capi.cpp -o obj/ganpcf_capi.o
	
obj/ganpcf.o: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.o

clean:
	rm obj/*.o
	rm libnpcf.mod
	rm emulator
