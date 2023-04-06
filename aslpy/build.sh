mkdir lib 
clang++ -fpic -shared -o lib/libaslpyc.so  src/aslpy_c.cpp -Iasl/build/include -Iasl/src/solvers -Iinclude -Llib -Lasl/build/lib -lasl