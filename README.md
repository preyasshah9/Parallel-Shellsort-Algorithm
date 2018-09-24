# Parallel-Shellsort-Algorithm
Code: proj_Preyas.c
Job Script File: jobscript
Make File: Makefile

These files are in the Shellsort folder.

Execution Instructions:
Place the following files into your home directory:
1. jobscript
2. Makefile
3. parallelShellsort.c

To generate the object file
Run Command: make

Execute the job using:
qsub -q wsuq jobscript

Generated output files:
output.txt (Used for debugging purposes only)
jobscript.o* (stdout)
jobscript.e* (stderr)

Output is stored in the jobscript.o** file.
output is of format: 
Execution Time: xxx
Line 1: output Array row 1
Line 2: output Array row 2
.
.
.
Line 8: output Array row 8
 
To change the array size: 
-> Go to the parallelShellsort.c code file,
-> Modify the constant definition: #define ARRAY_SIZE
   Set LOCAL_ARRAY_SIZE to be ARRAY_SIZE/PROC_NUM

To execute the other sorting algorithm, the corresponding changes needs to be done in the jobscript and makefile.
Code files are Quicksort_serial.c and odd_even_sort.c

NOTE: Execution time is printed in seconds. 
I made some changes to print the execution time in milliseconds for the smaller arrays.
However, I did not include these changes in the final code.
