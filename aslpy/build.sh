mkdir lib 
g++ -std=c++0x -fpic -shared -o lib/libaslpyc.so src/aslpy_c.cpp asl/build/lib/*.o -Iasl/build/include -Iasl/src/solvers -Iinclude -Llib