#Primero irian las librerias

#Si quiero extender el Makefile poner el +=
CXXFLAGS = -Wall -std=c++11 -g

#Poner luego un bin
BIN = test

#El resto de la ejecucion
release: CXXFLAGS += -O3 -DNDEBUG
release: test

#test

obj := ${patsubst %.cpp,%.o,${wildcard ReadData/*.cpp}}
ReadData/%.o: ReadData/%.cpp ReadData/*.h Utils/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

dbg_obj := ${patsubst %.cpp,%.o,${wildcard DBG/*.cpp}}
DBG/%.o: DBG/%.cpp DBG/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

test: ${obj} ${dbg_obj}
	${CXX} ${obj} ${dbg_obj} -o ${BIN} ${LDFLAGS}

#Clean
clean:
	-rm ${obj}
	-rm ${dbg_obj}
	-rm ${BIN}
