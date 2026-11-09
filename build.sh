mkdir -p build && cd build
cmake -DCUES_BUILD_TESTS=1 -DCMAKE_BUILD_TYPE=RELEASE ../cues
make
