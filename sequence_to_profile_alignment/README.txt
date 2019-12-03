Set environment variables, for example:
	gap_penalty=-2
	match_score=1
	mismatch_penalty=-1

Run program with:
	./alignSeqToProfile --fasta seq.fasta --aln aligned_sequences.aln --out seq.aln --gap ${gap_penalty} --match ${match_score} --mismatch ${mismatch_penalty}

OR simply
	./alignSeqToProfile --fasta seq.fasta --aln aligned_sequences.aln --out seq.aln --gap -2 --match 1 --mismatch -1