all:  example

example: mt2_bisect.cpp example.cpp mt2_bisect.h 
	g++ -o example mt2_bisect.cpp example.cpp 
ttbar:
	g++ -o CallSimulatorttbar mt2_bisect.cpp CallSimulatorttbar.C `root-config --cflags --libs` 
clean:
	rm -f example CallSimulatorttbar
tar:
	tar -cvzf mt2-1.01a.tar.gz mt2_bisect.cpp mt2_bisect.h Makefile example.cpp README update_notes
