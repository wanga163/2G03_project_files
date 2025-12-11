CXX  = c++
LIBS = -lcpgplot -lpgplot -lX11 -lgfortran -lm

OBJS = functions.o pk_model.o pk_simulation.o

all: pk_simulation

pk_simulation: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o pk_simulation

functions.o: functions.cpp pk_model.h
	$(CXX) -c functions.cpp

pk_model.o: pk_model.cpp pk_model.h
	$(CXX) -c pk_model.cpp

pk_simulation.o: pk_simulation.cpp pk_model.h
	$(CXX) -c pk_simulation.cpp

clean:
	rm -f $(OBJS) pk_simulation
