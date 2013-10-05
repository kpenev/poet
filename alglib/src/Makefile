CPP = g++
FLAGS = -g -Wall
ALL = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o

default: all

alglibinternal.o: alglibinternal.h alglibinternal.cpp ap.h stdafx.h
	$(CPP) $(FLAGS) -c -o alglibinternal.o alglibinternal.cpp

alglibmisc.o: alglibmisc.h alglibmisc.cpp ap.h alglibinternal.h stdafx.h
	$(CPP) $(FLAGS) -c -o alglibmisc.o alglibmisc.cpp

ap.o: ap.h ap.cpp stdafx.h
	$(CPP) $(FLAGS) -c -o ap.o ap.cpp

dataanalysis.o: dataanalysis.h dataanalysis.cpp ap.h alglibinternal.h linalg.h statistics.h alglibmisc.h specialfunctions.h solvers.h optimization.h stdafx.h
	$(CPP) $(FLAGS) -c -o dataanalysis.o dataanalysis.cpp

diffequations.o: diffequations.h diffequations.cpp ap.h alglibinternal.h stdafx.h
	$(CPP) $(FLAGS) -c -o diffequations.o diffequations.cpp

fasttransforms.o: fasttransforms.h fasttransforms.cpp alglibinternal.h ap.h stdafx.h
	$(CPP) $(FLAGS) -c -o fasttransforms.o fasttransforms.cpp

integration.o: integration.h integration.cpp ap.h alglibinternal.h linalg.h specialfunctions.h stdafx.h
	$(CPP) $(FLAGS) -c -o integration.o integration.cpp

interpolation.o: interpolation.h interpolation.cpp ap.h alglibinternal.h alglibmisc.h linalg.h solvers.h optimization.h specialfunctions.h integration.h stdafx.h
	$(CPP) $(FLAGS) -c -o interpolation.o interpolation.cpp

linalg.o: linalg.h linalg.cpp ap.h alglibinternal.h alglibmisc.h stdafx.h
	$(CPP) $(FLAGS) -c -o linalg.o linalg.cpp

optimization.o: optimization.h optimization.cpp ap.h alglibinternal.h alglibmisc.h linalg.h solvers.h stdafx.h
	$(CPP) $(FLAGS) -c -o optimization.o optimization.cpp
	
solvers.o: solvers.h solvers.cpp ap.h alglibinternal.h linalg.h alglibmisc.h stdafx.h
	$(CPP) $(FLAGS) -c -o solvers.o solvers.cpp
	
specialfunctions.o: specialfunctions.h specialfunctions.cpp ap.h alglibinternal.h stdafx.h
	$(CPP) $(FLAGS) -c -o specialfunctions.o specialfunctions.cpp

statistics.o: statistics.h statistics.cpp ap.h alglibinternal.h linalg.h specialfunctions.h stdafx.h
	$(CPP) $(FLAGS) -c -o statistics.o statistics.cpp

all: $(ALL)

clean:
	rm -f *.o
