include makefile.cudacompute
OS = $(shell lsb_release -si)
FC = gfortran
ifeq ($(OS),Fedora)
	CXX = cuda-g++
# 	FC = cuda-gfortran
else
	CXX = g++
endif


VXX = nvcc $(ARCHS) -ccbin=$(CXX)
FCFLAGS = -march=native -mtune=native -O3 -cpp
CXXFLAGS = -march=native -mtune=native -O3 -fPIC
VXXFLAGS = -Xptxas -dlcm=ca -lineinfo --compiler-options "$(CXXFLAGS)" -O3
LDFLAGS = -lstdc++ -lcudart -lcuda

library: obj/ganpcf_mod.o obj/ganpcf_capi.o obj/ganpcf.o
	$(FC) $(LDFLAGS) $^ -fPIC -shared -o libganpcf.so
	$(MAKE) test/emulator
	
obj/ganpcf_mod.o: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -fPIC -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.o
	
obj/ganpcf_capi.o: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) -c source/ganpcf_capi.cpp -o obj/ganpcf_capi.o
	
obj/ganpcf.o: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.o
	
test/emulator: source/emulator.f90
	$(shell cp libganpcf.so test)
	$(FC) $(FCFLAGS) -lganpcf source/emulator.f90 -o test/emulator

clean:
	rm obj/*.o
	rm libnpcf.mod
	rm libganpcf.so
	rm test/libganpcf.so
	rm test/emulator
