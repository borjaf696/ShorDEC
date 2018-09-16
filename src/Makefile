#Primero irian las librerias

#Si quiero extender el Makefile poner el +=
CXXFLAGS = -Wall -std=c++11 -g

#Poner luego un bin
BIN = test

#El resto de la ejecucion
release: CXXFLAGS += -O3 -DNDEBUG
release: test

#Read, seq, kmer
obj := ${patsubst %.cpp,%.o,${wildcard ReadData/*.cpp}}
ReadData/%.o: ReadData/%.cpp ReadData/*.h 
	${CXX} -c ${CXXFLAGS} $< -o $@

#Extender
extender_obj := ${patsubst %.cpp,%.o,${wildcard Extender/*.cpp}}
Extender/%.o: Extender/%.cpp Extender/*.h Utils/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Dbg
dbg_obj := ${patsubst %.cpp,%.o,${wildcard DBG/*.cpp}}
DBG/%.o: DBG/%.cpp DBG/*.h Utils/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

test: ${obj} ${dbg_obj} ${extender_obj}
	${CXX} ${obj} ${extender_obj} ${dbg_obj}  -o ${BIN} ${LDFLAGS}

#Clean
clean:
	-rm ${obj}
	-rm ${dbg_obj}
	-rm ${BIN}
