First, run
		make

For index generation, run
		./bwtfm index text.fa
This will automatically generate text.fa.bwt and text.fa.fm.

For search, run
		./bwtfm search text.fa pattern.fa

Alphabet supported is {A, C, G, T} (case sensitive)

Suffix Array generation is naive O(m^2logm)