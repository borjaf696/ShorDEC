CXXFLAGS = -std=c++17 -Wall -Wreorder -g

BIN_release = ../bin/viaDBG
BIN_test = ../bin/viaDBG_test

release: CXXFLAGS += -O3 -DNDEBUG
release: test

tests: CXXFLGS += -O3 -DNDEBUG
tests: test_environment

#Print
p_obj := ${patsubst %.cpp,%.o,${wildcard Utils/*.cpp}}
Utils/%.o: Utils/%.cpp Utils/*h
	${CXX} -c ${CXXFLAGS} ${LDFLAGS} $< -o $@

#Read, seq, kmer
obj := ${patsubst %.cpp,%.o,${wildcard ReadData/*.cpp}}
ReadData/%.o: ReadData/%.cpp ReadData/*.h Utils/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Extender
extender_obj := ${patsubst %.cpp,%.o,${wildcard Extender/*.cpp}}
Extender/%.o: Extender/%.cpp Extender/*.h 
	${CXX} -c ${CXXFLAGS} $< -o $@

#Dbg
dbg_obj := ${patsubst %.cpp,%.o,${wildcard DBG/*.cpp}}
DBG/%.o: DBG/%.cpp DBG/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@ ${OPENMPFLAG}

#Main_code
m_main_obj := ${patsubst %.cpp,%.o,${wildcard Main/*.cpp}}
Main/%.o: Main/%.cpp 
	${CXX} -c ${CXXFLAGS} $< -o $@
m_code_release := Main/main.o
m_code_tests := Main/test_environment.o


test: ${p_obj} ${obj} ${dbg_obj} ${extender_obj} ${m_main_obj}
	${CXX} ${p_obj} ${obj} ${extender_obj} ${dbg_obj} ${m_code_release} -o ${BIN_release} ${LDFLAGS}

test_environment: ${obj} ${dbg_obj} ${extender_obj} ${m_main_obj}
	${CXX} ${obj} ${extender_obj} ${dbg_obj} ${m_code_tests} -o ${BIN_test} ${LDFLAGS}

#Clean
clean:
	-rm ${obj}
	-rm ${p_obj}
	-rm ${dbg_obj}
	-rm ${extender_obj}
	-rm ${m_main_obj}
