git clone https://github.com/ampl/asl.git asl
cd asl
mkdir build
cd build
cmake .. -DUSE_LTO=OFF
make
cd lib
ar -x libasl.a