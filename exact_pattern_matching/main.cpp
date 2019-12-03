/**
 * @file main.cpp
 * @brief This program compares exact string matching algorithms
 * @author Naisila Puka
 * @version 20/10/2019
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

# define ASCII 128

// preprocessing failure function
int * failure(const std::string& pattern) {
	// like in lecture slides
    int m = pattern.size();
    int * failure = new int[m];
    failure[0] = 0;
    int i = 1;
    int j = 0;

    while(i < m) {
        if(pattern[i] == pattern[j]) {
            failure[i] = j + 1;
            i = i + 1;
            j = j + 1;
        }
        else if(j > 0)
            j = failure[j - 1];
        else {
            failure[i] = 0;
            i = i + 1;
        }
    }
    return failure;
}

// preprocessing for bad character rule
int ** badChar(const std::string& pattern, const std::string& alph) {
	int** bc = new int*[ASCII];
	int a = alph.size();
	int m = pattern.size();
	for(int i = 0; i < a; i++) {
		bc[int(alph[i])] = new int[m + 1];
		bc[int(alph[i])][0] = 0;
	}

	for(int i = 1; i <= m; i++) {
		for(int j = 0; j < a; j++) {
			if(int(alph[j]) == int(pattern[i - 1]))
				bc[int(alph[j])][i] = i;
			else
				bc[int(alph[j])][i] = bc[int(alph[j])][i - 1];
		}
	}
	return bc;
}

// preprocessing for good suffix rule 1
void gs1(const std::string& pattern, int *gs, int *failure) {
	int m = pattern.size();
	failure[m] = m + 1; //failure is nothing

	//1-based coordinate
    int j = m + 1;

    for (int i = m; i > 0; i--){
        // keep going, no failures updated
        while(j <= m && pattern[i - 1] != pattern[j - 1]) {
            // since z != y (from slide 32 in 02-exact-kmp-boyer_moore.pdf)
            if(gs[j] == 0)
                gs[j] = j - i;
            j = failure[j]; //position of next failure
        }

        // pattern[i - 1] == pattern[j - 1], store new failure
        failure[i - 1] = j - 1;
        j--;
    }
}

//Preprocessing for good suffix rule 2 combined with 1
void gs2(const std::string& pattern, int *gs, int *failure) {
	int m = pattern.size();
    int gs2 = failure[0];
    for(int i = 0; i <= m; i++) {
        // use good suffix rule 2 (gsr2) if no appropriate substring is obtained from gsr1
        // since we can use failure function for a prefix of Pattern, which
        // matches with suffix of text in current position
        if(gs[i] == 0)
            gs[i] = gs2;

        // suffix shorter than current failure, make j next failure
        if (i == gs2)
            gs2 = failure[gs2];
    }
}

int main(int argc, char **argv) {
	printf("----------------------\n");
	printf("Exact Pattern Matching\n");
	printf("----------------------\n");

    if (argc < 3) {
        printf("You need to provide both text and pattern files in FASTA format!\n");
        return -1;
    }

    std::ifstream input_text(argv[1]);
    std::ifstream input_pattern(argv[2]);

    if (!input_text.good() || !input_pattern.good()) {
        printf("Error opening text and pattern files!\n");
        return -1;
    }

    std::string line, text, pattern, best;

    while (std::getline(input_text, line)) {
        if(line.empty() || line[0] == '>')
            continue;
        else
            text += line;
    }

    while (std::getline(input_pattern, line)) {
        if(line.empty() || line[0] == '>')
            continue;
        else
            pattern += line;
    }

    // printf("\nText: %s\n\nPattern: %s\n", text.c_str(), pattern.c_str());

    int n = text.size();
    int m = pattern.size();
    if(n < m) {
        printf("P is not in T\n");
        return 0;
    }

    printf("\nNotes: if nothing is printed under each");
    printf("\nalgorithm's name then P is not in T\n");
    printf("Runtimes and comparisons don't include preprocessing.\n");

    int min;

    // brute force exact string matching
    printf("\n>Brute Force\n");

    int i = 1;
    int check = n - m + 2; // since we are using 1-based coordinate
    int j;
    int cmp = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    while(i < check) {
        for(j = 0; j < m; j++) {
            cmp++;
            if(pattern[j] != text[i + j - 1]) {
                break;
            }
        }
        if(j == m) {
            printf("P is in T at position %d. (%d comparisons)\n", i, cmp);
            break;
        }
        i++;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << runtime << " microseconds" << std::endl;
    min = int(runtime);
    best = "Brute Force";

    // Knuth-Morris-Pratt
    printf("\n>Knuth-Morris-Pratt\n");

    int * f = failure(pattern);
    i = 0;
    j = 0;
    cmp = 0;
    t1 = std::chrono::high_resolution_clock::now();
    while(i < n) {
        cmp++;
        if(text[i] == pattern[j]) {
            if(j == m - 1) {
                printf("P is in T at position %d. (%d comparisons)\n", i - j + 1, cmp);
                break;
            }
            else {
                i++;
                j++;
            }
        }
        else {
            if(j > 0)
                j = f[j - 1];
            else {
                i++;
                j = 0;
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << runtime << " microseconds" << std::endl;
    if((int)(runtime) < min) {
    	min = (int)runtime;
    	best = "Knuth-Morris-Pratt";
    }

    // Rabin-Karp
    printf("\n>Rabin-Karp\n");

    int * text_int = new int[n];
    int * pattern_int = new int[m];

    for(int k = 0; k < n; k++) {
        if(text[k] == 'A')
            text_int[k] = 0;
        else if(text[k] == 'C')
            text_int[k] = 1;
        else if(text[k] == 'G')
            text_int[k] = 2;
        else if(text[k] == 'T')
            text_int[k] = 3;
    }

    for(int k = 0; k < m; k++) {
        if(pattern[k] == 'A')
            pattern_int[k] = 0;
        else if(pattern[k] == 'C')
            pattern_int[k] = 1;
        else if(pattern[k] == 'G')
            pattern_int[k] = 2;
        else if(pattern[k] == 'T')
            pattern_int[k] = 3;
    }

    int q = 997; //fix this, always pick q > m
    int c = 1;
    for(int h = 0; h < m - 1; h++)
    	c = (c * 4) % q;
    int fp = 0;
    int ft = 0;
    cmp = 0;
    for(i = 0; i < m; i++) {
        fp = (4 * fp + pattern_int[i]) % q;
        ft = (4 * ft + text_int[i]) % q;
    }

    t1 = std::chrono::high_resolution_clock::now();
    // for loop as in slides
    for(int s = 0; s < n - m + 1; s++) {
        if(fp == ft) {
        	int g;
            for(g = 0; g < m; g++) {
                cmp++;
                if(pattern_int[g] != text_int[s + g]) {
                    break;
                }
            }
            if(g == m) {
                printf("P is in T at position %d. (%d comparisons)\n", s + 1, cmp);
                break;
            }
        }
        ft = ((ft - text_int[s] * c) * 4 + text_int[s + m]) % q;
        if (ft < 0)
            ft = ft + q;
    }
    t2 = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << runtime << " microseconds" << std::endl;
    if((int)(runtime) < min) {
    	min = (int)runtime;
    	best = "Rabin-Karp";
    }

    // Boyer-Moore
    printf("\n>Boyer-Moore\n");

    // algorithm similar to failure function
    // failure will be calculated inside gs1 and used in gs2
    int failure[m + 1], gs[m + 1];

    //bad character preprocessing
    int ** bc = badChar(pattern, "ACTG");

    //good suffix preprocessing
    for(int i = 0; i < m + 1; i++)
    	gs[i] = 0;
    gs1(pattern, gs, failure);
    gs2(pattern, gs, failure);

    int s = 0; //position in text
    cmp = 0;

    t1 = std::chrono::high_resolution_clock::now();
    while(s < n - m + 1) {
  		j = m - 1; // left to right comparison, position in pattern

        // check until mismatch
        while(j >= 0 && (pattern[j] == text[s + j])) { // s + j in text since we have left to right comparison
            cmp++;
            j--;
        }

        if (j < 0) {
            printf("P is in T at position %d. (%d comparisons)\n", s + 1, cmp);
            break;
        }

        // mismatch, pick max(bad character, min(gs1, gs2))
        else {
        	cmp++; // count the mismatch!
        	// (+1)'s below since in preprocessing algorithms, 1-based coordinates are used
        	int gsu = gs[j + 1]; // amount of shift we get from good suffix rules
        	int bch = j + 1 - bc[int(text[s+j])][j + 1]; // from bad character rule
            s += fmax(gsu, bch);
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << runtime << " microseconds" << std::endl;
    if((int)(runtime) < min) {
    	min = (int)runtime;
    	best = "Boyer-Moore";
    }

    printf("\nBest runtime %d microseconds: %s\n", min, best.c_str());
    return 0;
}
