# Copyright (C) 2012-2017 IBM Corp.
#
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.
# 
CC = g++ 
LD = g++
AR = ar

$(info Required NTL version 10.4.0 or later. You should generate dependency)
$(info file makedep.d before compiling by running 'make depend'.)
$(info )

CFLAGS = -g -O2 -march=native -Wfatal-errors -Wshadow -Wall -I. -pthread -std=c++11 
ARFLAGS =ruv
LDFLAGS = -L/usr/local/lib 
LDLIBS  = -lntl -lgmp -lm

OBJ = obfus.o GGH15.o TDMatrixParams.o TDMatrix.o CRTmatrix.o DGaussSampler.o Gaussian1Dsampler.o sampleG.o mat_l.o vec_l.o utils/argmap.o utils/tools.o utils/timing.o GGHnodes.o

PROGS = initialize_x obfuscate_x evaluate_x

TESTPROGS = eval_timing_x init_timing_x test_GenMatPair_x testOBF_x test_threads_x makeRandomBP_x test_CRTmatrix_x test_GGH15_x test_sampleG_x test_timing_x test_DGSamp_x test_IO_x test_TDmatrix_x

all: obf.a $(PROGS)

obf.a: $(OBJ)
	$(AR) $(ARFLAGS) obf.a $(OBJ)

# old-style make depend (FIXME: change to modern style?)
depend:
	g++ -std=gnu++11 -MM *.cpp > makedep.d

# if makedep.d exists in the current directory, include it
include $(wildcard makedep.d)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

./%_x: programs/%.cpp obf.a
	$(CC) $(CFLAGS) -o $@ $< obf.a $(LDFLAGS) $(LDLIBS)

test: $(TESTPROGS)
	./test_CRTmatrix_x threads=4
	./test_DGSamp_x
	./test_sampleG_x
	./test_TDmatrix_x
	./test_GGH15_x threads=4
	./test_IO_x
	./testOBF_x pipe=1 threads=4

clean:
	rm -f *.o *_x *_x.exe *.a core.* utils/*.o programs/*_x programs/*_x.exe
	rm -rf *.dSYM
