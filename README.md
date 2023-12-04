# Fresh-AIR
FreshAIR: Adaptive Spatial Indexing in Dynamic Environments 

This is an implementation of the proposed update-able spatial adaptive index, FreshAIR.

## Dependencies

- `gcc`

## Usage

The index is defined and implemented in the `fresh_air.h` file. The runner file, `general_testing.cpp` includes the header, creates an instance of the tree, and tests its use. This can be tuned using many flags at compile time, and input values at runtime. An example of the compilation command is as follows:  

```
g++ -D data_size6M -D versionGSR -D fo16 -D maxt64 -D mint32 -D SSL8 -D dh0 -D mt128 -mcmodel=large general_testing_fresh_air_deletion.cpp -O3 -o fresh
```

Then the running code will be:

```
./fresh data.txt 100000 query.txt times.txt 100 2000 
```

## Compile time configurations

The flags that you can use at compile time to set parameters or perform other tasks, are as follows:

- `-D datasizeX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this determines the size of the pre-loaded data.
