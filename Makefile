################################
# Makefile
#
# author: Sizhen Li
# edited by: 08/2021
################################

CC=g++
SANKOFFDEPS=src/LinearSankoffInterface.h src/LinearSankoff.h src/LinearSankoff_array.h src/HMMAlign.h src/LinearFold.h src/LinearPartition.h src/bpp.cpp src/check_mem.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
MULTILIGNDEPS=src/LinearMultilign.h src/LinearSankoff.h src/HMMAlign.h src/LinearFold.h src/check_mem.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
ALIGNDEPS=src/HMMAlignInterface.h src/HMMAlign.h src/check_mem.h src/Utils/utility.h
CFLAGS=-std=c++11 -O3
.PHONY : clean bin/*
objects=bin/linearmultilign bin/linearsankoff bin/linearsankoff_dynalign bin/linearalignment

linearsankoff: src/LinearSankoffInterface.cpp $(SANKOFFDEPS)
		chmod +x linearsankoff
		mkdir -p bin
		$(CC) src/LinearSankoffInterface.cpp src/LinearSankoff.cpp src/LinearSankoff_array.cpp src/HMMAlign.cpp src/LinearFold.cpp src/LinearPartition.cpp src/check_mem.cpp src/Utils/utility.cpp src/Utils/energy_parameter.cpp src/Utils/feature_weight.cpp src/Utils/intl11.cpp src/Utils/intl21.cpp src/Utils/intl22.cpp $(CFLAGS) -Dlv -o bin/linearsankoff
		# $(CC) src/LinearSankoffInterface.cpp src/LinearSankoff.cpp src/HMMAlign.cpp src/LinearFold.cpp src/LinearPartition.cpp src/check_mem.cpp src/Utils/utility.cpp src/Utils/energy_parameter.cpp src/Utils/feature_weight.cpp src/Utils/intl11.cpp src/Utils/intl21.cpp src/Utils/intl22.cpp $(CFLAGS) -Dlv -Ddynalign -o bin/linearsankoff_dynalign

linearmultilign: src/LinearMultilign.cpp $(MULTILIGNDEPS)
		chmod +x linearmultilign
		mkdir -p bin
		$(CC) src/LinearMultilign.cpp src/LinearSankoff.cpp src/HMMAlign.cpp src/LinearFold.cpp src/check_mem.cpp src/Utils/utility.cpp src/Utils/energy_parameter.cpp src/Utils/feature_weight.cpp src/Utils/intl11.cpp src/Utils/intl21.cpp src/Utils/intl22.cpp $(CFLAGS) -Dlv -Dmultilign -o bin/linearmultilign

linearalignment: src/HMMAlignInterface.cpp $(ALIGNDEPS)
		chmod +x linearalignment
		mkdir -p bin
		$(CC) src/HMMAlignInterface.cpp src/HMMAlign.cpp src/Utils/utility.cpp $(CFLAGS) -o bin/linearalignment

clean:
	-rm $(objects)