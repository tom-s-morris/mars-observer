all:	MarsOpp

OBJS 		= MarsOpp.o LibEphcom4.o MarsObserver.o MarsDataLogger.o LReader.o \
	EventsFinder.o
CXXFLAGS	= -g -std=c++11 -I../ -I../include -fPIC
#JPLEPH_VER1 = -L../lib -ljpleph
JPLEPH_VER1 =

MarsOpp: $(OBJS)
	g++ -g -std=c++11 -o $@ $(JPLEPH_VER1) -ldl -llua $(OBJS)

clean:
	rm -f $(OBJS) MarsOpp

