/**
 * @file alignSeqToProfile.cpp
 * @brief Sequence to Profile Alignment
 * @author Naisila Puka
 * @version 30/11/2019
*/

#include <iostream>
#include <unistd.h>
#include <map>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include <getopt.h>
#include <chrono>

int main(int argc, char** argv) {
	struct option opts[] =
	{
		{ "fasta", required_argument, 0, 0},
		{ "aln", required_argument, 0, 1},
		{ "out", required_argument, 0, 2},
		{ "gap", required_argument, 0, 3},
		{ "match", required_argument, 0, 4},
		{ "mismatch", required_argument, 0, 5},
		{ 0, 0, 0, 0}
	};

	int c;
	int option_index = 0, all = 0;
	char* fastaFileName = NULL, *alnFileName = NULL, *outFileName = NULL;
	int match, mismatch, gap;

	while (1) {
	    c = getopt_long(argc, argv, "", opts, &option_index);
	    if (c == -1)
	    	break;
	    switch (c) {
	    	case 0:
	    		fastaFileName = optarg;
	    		all++;
	    		break;
	    	case 1:
	    		alnFileName = optarg;
	    		all++;
	    		break;
	    	case 2:
	    		outFileName = optarg;
	    		all++;
	    		break;
	    	case 3:
	    		gap = atoi(optarg);
	    		all++;
	    		break;
	    	case 4:
	    		match = atoi(optarg);
	    		all++;
	    		break;
	    	case 5:
	    		mismatch = atoi(optarg);
	    		all++;
	    		break;
	    }
	}

	if (all < 5) {
		printf("Error! Not all arguments are provided!\n");
		return -1;
	}

	std::string line;

	std::string seqs[10];
	std::string names[10];
	//n: number of sequences, 2<=n<=10, l: profile length, m: newseq length
	int n = 0, l = 0, m = 0;

	std::string fasta, seq, fastaname;
	fasta += fastaFileName;
	std::string aln;
	aln += alnFileName;

    //get sequence
	std::ifstream input_fasta(fasta.c_str());
	while (std::getline(input_fasta, line)) {
		if (line.empty() || line[0] == '>') {
			if (line[0] == '>' && line.size() > 1) {
				fastaname += line.substr(1, line.size() - 1);
			}
		}
		else {
			seq += line;
		}
	}

	//get alignment sequences
	std::ifstream input_aln(aln.c_str());
	while (std::getline(input_aln, line)) {
		if (line.empty()) {
			continue;
		}
		else {
			seqs[n] += line;
			n++;
		}
	}

	int first, last;
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < seqs[i].size(); j++) {
			if ((seqs[i][j] == ' ' ||  seqs[i][j] == '\t') && j != seqs[i].size() - 1) {
				first = j;
				break;
			}
		}
		for (int j = 0; j < seqs[i].size(); j++) {
			if ((seqs[i][j] == ' ' ||  seqs[i][j] == '\t') && j != seqs[i].size() - 1)
				last = j;
		}
		names[i] += seqs[i].substr(0, first);
		seqs[i] = seqs[i].substr(last + 1, seqs[i].size() - last - 1);
	}

	// just to make names with the same length
	int * name_lengths = new int[n + 1];
	for(int g = 0; g < n; g++)
		name_lengths[g] = names[g].size();
	name_lengths[n] = fastaname.size();
	int* max_length;
    max_length = std::max_element(name_lengths, name_lengths + n + 1);
    int max = *max_length;
    for(int g = 0; g < n; g++) {
    	int current = names[g].size();
    	for(int i = 0; i < (max - current); i++)
			names[g] += " ";
    }
    int current = fastaname.size();
    for(int i = 0; i < (max - current); i++)
		fastaname += " ";

	l = seqs[0].size();
	m = seq.size();

    //profile matrix
	std::map<char, float*> P;
	P['A'] = new float[l];
	P['C'] = new float[l];
	P['G'] = new float[l];
	P['T'] = new float[l];
	P['-'] = new float[l];
	for(int i = 0; i < l; i++) {
		P['A'][i] = 0;
		P['C'][i] = 0;
		P['G'][i] = 0;
		P['T'][i] = 0;
		P['-'][i] = 0;
	}

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < n; j++) {
			P[seqs[j][i]][i] += 1/((float)(n));
		}
	}

    // score for aligning char with column int
	std::map<char, float*> S;
	S['A'] = new float[l];
	S['C'] = new float[l];
	S['G'] = new float[l];
	S['T'] = new float[l];
	S['-'] = new float[l];

	for(int i = 0; i < l; i++) {
		S['A'][i] = match * P['A'][i] + mismatch * P['C'][i] + mismatch * P['G'][i] + mismatch * P['T'][i] + gap * P['-'][i];
		S['C'][i] = mismatch * P['A'][i] + match * P['C'][i] + mismatch * P['G'][i] + mismatch * P['T'][i] + gap * P['-'][i];
		S['G'][i] = mismatch * P['A'][i] + mismatch * P['C'][i] + match * P['G'][i] + mismatch * P['T'][i] + gap * P['-'][i];
		S['T'][i] = mismatch * P['A'][i] + mismatch * P['C'][i] + mismatch * P['G'][i] + match * P['T'][i] + gap * P['-'][i];
		S['-'][i] = gap;
	}

	// v(m, l) is the sequence to profile alignment score
	float ** v = new float*[m + 1];
	for(int i = 0; i <= m; i++) {
		v[i] = new float[l + 1];
	}

	// backtrack
	int ** b = new int*[m + 1];
	for(int i = 0; i <= m; i++) {
		b[i] = new int[l + 1];
	}

	v[0][0] = 0;
	for(int j = 1; j <= l; j++) {
		v[0][j] = v[0][j - 1] + S['-'][j - 1];
	}
	for(int i = 1; i <= m; i++) {
		v[i][0] = v[i - 1][0] + gap;
	}

	// B(i,j) is the arrow showing the way diag: 0, up: 1, left: -1
	b[0][0] = -1;
	for(int j = 1; j <= l; j++) {
		b[0][j] = -1;
	}
	for(int i = 1; i <= m; i++) {
		b[i][0] = 1;
	}

	float diag, up, left;
	auto t1 = std::chrono::high_resolution_clock::now();
	for(int i = 1; i <= m; i++) {
		for(int j = 1; j <= l; j++) {
			diag = v[i - 1][j - 1] + S[seq[i - 1]][j - 1];
			up = v[i - 1][j] + gap;
			left = v[i][j - 1] + S['-'][j - 1];

			if(diag >= up && diag >= left) {
				v[i][j] = diag;
				b[i][j] = 0;
			}
			else if(up >= diag && up >= left) {
				v[i][j] = up;
				b[i][j] = 1;
			}
			else {
				v[i][j] = left;
				b[i][j] = -1;
			}
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	float score = v[m][l];
	auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "\nSequence to profile alignment score calculated in " << runtime << " microseconds" << std::endl;
    printf("Score is %.2f\n", score);

	// reconstructing the solution
	int k = m + l;
	int i = m;
	int j = l;

	int xpos = k;
	int ypos[10];
	for(int g = 0; g < n; g++)
		ypos[g] = k;

	char ** old_seqs = new char*[n];
	for(int g = 0; g < n; g++)
		old_seqs[g] = new char[k + 1];

	char * new_seq = new char[k + 1];

	while (!(i == 0 || j == 0)) {
		if (b[i][j] == 0) {
			new_seq[xpos--] = seq[i - 1];
			for(int g = 0; g < n; g++) {
				old_seqs[g][ypos[g]--] = seqs[g][j - 1];
			}
			i--; j--;
		}
		else if (b[i][j] == 1) {
			new_seq[xpos--] = seq[i - 1];
			for(int g = 0; g < n; g++)
				old_seqs[g][ypos[g]--] = '-';
			i--;
		}
		else if (b[i][j] == -1) {
			new_seq[xpos--] = '-';
			for(int g = 0; g < n; g++)
				old_seqs[g][ypos[g]--] = seqs[g][j - 1];
			j--;
		}
	}

	while (xpos > 0) {
		if (i > 0) new_seq[xpos--] = seq[--i];
		else new_seq[xpos--] = '-';
	}
	while (ypos[0] > 0) {
		if (j > 0) {
			for(int g = 0; g < n; g++)
				old_seqs[g][ypos[g]--] = seqs[g][--j];
		}
		else {
			for(int g = 0; g < n; g++)
				old_seqs[g][ypos[g]--] = '-';
		}
	}

	int skip = 1;
	int check = 1;
	for(int h = 1; h <= k; h++) {
		for(int g = 0; g < n; g++) {
			if(old_seqs[g][h] != '-') {
				check = 0;
				break;
			}
		}
		if(new_seq[h] == '-' && check == 1)
			skip++;
		else
			break;
	}

	std::string out;
	out += outFileName;

	std::ofstream myfile(out.c_str());
	if (myfile.is_open()) {
		for(int g = 0; g < n; g++) {
			myfile << names[g];
			myfile << " ";
			for(int i = skip; i <= k; i++)
				myfile << old_seqs[g][i];
			myfile << "\n";
		}
		myfile << fastaname;
		myfile << " ";
		for(int i = skip; i <= k; i++)
			myfile << new_seq[i];
		myfile << "\n";
		myfile.close();
	}
	else
		std::cout << "Unable to open file";
	return 0;
}