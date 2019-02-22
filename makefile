ARCHS = -gencode arch=compute_50,code=sm_50 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75
VXX = nvcc $(ARCHS)
ifeq ($(OS),Windows_NT)
	FC = gfortran
	CXX = cl
	FCFLAGS = -march=native -mtune=native -O3
	CXXFLAGS = -O2 -D WIN32 -D_USE_MATH_DEFINES
	LDFLAGS = -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\lib\x64" -lcuda -lcudart -lgfortran
	CXXINCLUDE = -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include"
	OBJECTS = obj/ganpcf_mod.obj obj/ganpcf_capi.obj obj/ganpcf.obj obj/emulator.obj
	ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
		CXXFLAGS += -D AMD64
		FCFLAGS += -D AMD64
	else
		ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
			CXXFLAGS += -D AMD64
			FCFLAGS += -D AMD64
		endif
		ifeq ($(PROCESSOR_ARCHITECTURE),x86)
			CXXFLAGS += -D IA32
			FCFLAGS += -D IA32
		endif
	endif
else
	VXX += -ccbin=cuda-g++
	FC = gfortran
	CXX = g++
	FCFLAGS = -march=native -mtune=native -O3
	CXXFLAGS = -march=native -mtune=native -O3
	LDFLAGS = -lstdc++ -lgfortran -lcuda -lcudart
	OBJECTS = obj/ganpcf_mod.o obj/ganpcf_capi.o obj/ganpcf.o obj/emulator.o
endif
VXXFLAGS = -Xptxas -dlcm=ca -lineinfo --compiler-options "$(CXXFLAGS)" -O3

emulator: $(OBJECTS)
ifeq ($(OS),Windows_NT)
	$(VXX) $(CXXFLAGS) $(LDFLAGS) $^ -o emulator.exe
else
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o emulator
endif
	
obj/emulator.obj: source/emulator.f90
	$(FC) $(FCFLAGS) -c source/emulator.f90 -o obj/emulator.obj
	
obj/emulator.o: source/emulator.f90
	$(FC) $(FCFLAGS) -c source/emulator.f90 -o obj/emulator.o
	
obj/ganpcf_mod.obj: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.obj
	
obj/ganpcf_mod.o: source/ganpcf_mod.f90
	$(FC) $(FCFLAGS) -c source/ganpcf_mod.f90 -o obj/ganpcf_mod.o

obj/ganpcf_capi.obj: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c source/ganpcf_capi.cpp -Fo"obj\ganpcf_capi.obj"

obj/ganpcf_capi.o: source/ganpcf_capi.cpp
	$(CXX) $(CXXFLAGS) -c source/ganpcf_capi.cpp -o obj/ganpcf_capi.o
	
obj/ganpcf.obj: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.obj
	
obj/ganpcf.o: source/ganpcf.cu
	$(VXX) $(VXXFLAGS) -dw source/ganpcf.cu -o obj/ganpcf.o

clean:
ifeq ($(OS),Windows_NT)
	rm obj/*.obj
	rm emulator.exp
	rm emulator.lib
else
	rm obj/*.o
endif
	rm libnpcf.mod
	rm emulator
