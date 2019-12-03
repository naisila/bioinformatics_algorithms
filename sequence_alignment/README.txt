Invoke with:
	./allalign --mode global --input sequences.fasta --gapopen -5
	./allalign --mode aglobal --input sequences.fasta --gapopen -5 --gapext -2 allalign --mode local --input sequences.fasta --gapopen -5
	./allalign --mode alocal --input sequences.fasta --gapopen -5 --gapext -2

NOTE: there is no default state of arguments,
	  they have to always be provided

NOTE: arguments parsed in an "ugly" way.
	  Started parsing them with gflags
	  then I gave up since you might not
	  have gflags installed.