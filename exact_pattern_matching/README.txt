General:
	To run the program, execute the following commands:
	$ make
	$ ./main text.fa pattern.fa

Notes:
	Best algorithm is printed in the last line of output.
	Runtimes and comparisons don't include preprocessing.
	g++ compiler is used and support for ISO C++ 2011 standard.

For Rabin-Karp algorithm:
	- q = 997, change it to a larger prime if length of pattern >= 997
	- alphabet is assumed "ACGT"

For Boyer-Moore:
	- alphabet is assumed "ACGT"