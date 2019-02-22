ARCHS = -gencode arch=compute_50,code=sm_50 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75
VXX = nvcc $(ARCHS)
ifeq ($(OS),Windows_NT)
	FC = gfortran
	CXX = cl
	CXXFLAGS = -O2 -D WIN32 
	CXXINCLUDE = -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include"
	ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
		CXXFLAGS += -D AMD64
	else
		ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
			CXXFLAGS += -D AMD64
		endif
		ifeq ($(PROCESSOR_ARCHITECTURE),x86)
			CXXFLAGS += -D IA32
		endif
	endif
else
	VXX += -ccbin=cuda-g++
	FC = gfortran
	CXX = g++
	CXXFLAGS = -march=native -mtune=native -O3
endif
FCFLAGS = -march=native -mtune=native -O3
VXXFLAGS = -Xptxas -dlcm=ca -lineinfo --compiler-options "$(CXXFLAGS)" -O3
LDFLAGS = -lstdc++ -lgfortran -lcuda -lcudart

emulator: obj/ganpcf_mod.o obj/ganpcf_capi.o obj/ganpcf.o obj/emulator.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o emulator
	
obj/emulator.o: source/emulator.f90
	$(FC) $(FCFLAGS) -c source/emulator.f90 -o obj/emulator.o
	
obj/ganpcf_mod.o: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.o

if ($(OS),Windows_NT)
obj/ganpcf_capi.o: souce/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c source/ganpcf_capi.cpp /OUT:"obj/ganpcf_capi.o"
else
obj/ganpcf_capi.o: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) -c source/ganpcf_capi.cpp -o obj/ganpcf_capi.o
endif
	
obj/ganpcf.o: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.o

clean:
	rm obj/*.o
	rm libnpcf.mod
	rm emulator
