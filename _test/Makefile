CXX=g++
CXXFLAGS=-Wall -Wextra --std=c++14
COMPILE=$(CXX) $(CXXFLAGS)

.PHONY: build clean

build: test1 test2 test3 test4 test5 test6 test7

test1: test1.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test2: test2.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test3: test3.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test4: test4.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test5: test5.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test6: test6.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

test7: test7.cpp
	$(COMPILE) -I.. -L.. $< -o $@ -laegraph

clean:
	rm test1 test2 test3 test4 test5 test6 test7
