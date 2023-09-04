# Characteristic function vector
Vector index 0, position 1 -> binary (...00001) -> players (1)
Vector index 1, position 2 -> binary (...00010) -> players (2)
Vector index 2, position 3 -> binary (...00011) -> players (1, 2)
Vector index 3, position 4 -> binary (...00100) -> players (3)
Vector index 4, position 5 -> binary (...00101) -> players (1, 3)
Vector index 5, position 6 -> binary (...00110) -> players (2, 3)
Vector index 6, position 7 -> binary (...00111) -> players (1, 2, 3)
Vector index 7, position 8 -> binary (...01000) -> players (4)
Vector index 8, position 9 -> binary (...01001) -> players (1, 4)

# Build

- install cmake
- install emscripten sdk
- set toolchain to PATH_TO_EMSCRIPTEN_SDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake
- build with cmake

```
cmake --build .
```
