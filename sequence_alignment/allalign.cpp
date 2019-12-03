/**
 * @file allalign.cpp
 * @brief Sequence Alignment
 * @author Naisila Puka
 * @version 15/11/2019
*/

#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>

int main(int argc, char** argv) {
	if(argc < 7) {
		printf("ERROR! Not enough arguments provided!");
		return -1;
	}
	std::string task;
	task += argv[2];

	std::string s1, s2;
	std::string s1name, s2name;
	std::string line;
	int count = 0;
	int gap = atoi(argv[6]);

    //get sequences and their names
	std::ifstream input_text(argv[4]);
	while (std::getline(input_text, line)) {
		if (line.empty() || line[0] == '>') {
			if (line[0] == '>' && line.size() > 1) {
				if(count == 1)
					s2name += line.substr(1, line.size() - 1);
				else 
					s1name += line.substr(1, line.size() - 1);
				count++;
			}
			continue;
		}
		else {
			if(count == 1)
				s1 += line;
			else
				s2 += line;
		}
	}
	if(s1name.size() > s2name.size()) {
		for(int i = 0; i < (s1name.size() - s2name.size()); i++)
			s2name += " ";
	}
	else if(s1name.size() < s2name.size()) {
		for(int i = 0; i < (-(s1name.size() - s2name.size())); i++)
			s1name += " ";
	}

    // printf("s1name: %s\n", s1name.c_str());
    // printf("s1: %s\n", s1.c_str());
    // printf("s2name: %s\n", s2name.c_str());
    // printf("s2: %s\n", s2.c_str());

	int m = s1.size();
	int n = s2.size();

	if (task.compare("global") == 0) {
		int ** score = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			score[i] = new int[n + 1];
		}

		for(int i = 0; i <= m; i++) {
			score[i][0] = i * gap;
		}

		for(int i = 0; i <= n; i++) {
			score[0][i] = i * gap;
		}

		int match;

		for(int i = 1; i <= m; i++) {
			for(int j = 1; j <= n; j++) {
				if(s1[i - 1] == s2[j - 1]) {
					match = 2;
				}
				else
					match = -3;
				score[i][j] = fmax(score[i - 1][j - 1] + match, score[i][j - 1] + gap);
				score[i][j] = fmax(score[i][j], score[i - 1][j] + gap);
			}
		}
		// printf("score is %d\n", score[m][n]);

		// reconstructing the solution
		int l = m + n;
		int i = m;
		int j = n;

		int xpos = l; 
    	int ypos = l; 
  
		char * s1ans = new char[l + 1];
		char * s2ans = new char[l + 1];

		while (!(i == 0 || j == 0)) { 
			if (s1[i - 1] == s2[j - 1])
				match = 2;
			else
				match = -3;
			if (score[i - 1][j - 1] + match == score[i][j]) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = s2[j - 1]; 
				i--; j--; 
			} 
			else if (score[i - 1][j] + gap == score[i][j]) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = '-'; 
				i--; 
			} 
			else if (score[i][j - 1] + gap == score[i][j]) { 
				s1ans[xpos--] = '-'; 
				s2ans[ypos--] = s2[j - 1]; 
				j--; 
			} 
		}
		while (xpos > 0) { 
			if (i > 0) s1ans[xpos--] = s1[--i]; 
			else s1ans[xpos--] = '-'; 
		} 
		while (ypos > 0) { 
			if (j > 0) s2ans[ypos--] = s2[--j]; 
			else s2ans[ypos--] = '-'; 
		}

		int skip = 1;
		for(int k = 1; k <= l; k++) {
			if(s1ans[k] == '-' && s2ans[k] == '-')
				skip++;
			else
				break;
		}

		// for(int k = skip; k <= l; k++)
		// 	printf("%c", s1ans[k]);
		// printf("\n");
		// for(int k = skip; k <= l; k++)
		// 	printf("%c", s2ans[k]);
		// printf("\n");

		std::string outbwt;
        outbwt = "global-naiveGap.aln";

        std::ofstream myfile(outbwt);
        if (myfile.is_open()) {
        	myfile << "Score = ";
        	myfile << score[m][n];
        	myfile << "\n\n";

        	int i = skip;
        	int j = skip;
        	int count = 0;
        	while(true) {
        		int now = i;
        		myfile << s1name;
        		myfile << " ";
        		for(; i <= l && i < (now + 60); i++)
        			myfile << s1ans[i];
        		myfile << "\n";
        		myfile << s2name;
        		myfile << " ";
        		for(; j <= l && j < (now + 60); j++)
        			myfile << s2ans[j];
        		myfile << "\n\n";
        		if(i > l)
        			break;
        	}
            myfile.close();
        }
        else
            std::cout << "Unable to open file";

		return 0;
	}
	else if (task.compare("local") == 0){
		int ** score = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			score[i] = new int[n + 1];
		}

		for(int i = 0; i <= m; i++) {
			score[i][0] = 0;
		}

		for(int i = 0; i <= n; i++) {
			score[0][i] = 0;
		}

		int match;
		int max = 0;
		int indexi = 0;
		int indexj = 0;

		for(int i = 1; i <= m; i++) {
			for(int j = 1; j <= n; j++) {
				if(s1[i - 1] == s2[j - 1]) {
					match = 2;
				}
				else
					match = -3;
				score[i][j] = fmax(score[i - 1][j - 1] + match, score[i][j - 1] + gap);
				score[i][j] = fmax(score[i][j], score[i - 1][j] + gap);
				score[i][j] = fmax(score[i][j], 0);
				if(score[i][j] > max) {
					max = score[i][j];
					indexi = i;
					indexj = j;
				}
			}
		}
		// printf("score is %d\n", score[indexi][indexj]);

		// reconstructing the solution
		int l = indexi + indexj;
		int i = indexi;
		int j = indexj;

		int xpos = l; 
    	int ypos = l; 
  
		char * s1ans = new char[l + 1];
		char * s2ans = new char[l + 1];

		while (!(i == 0 || j == 0 || score[i][j] == 0)) { 
			if (s1[i - 1] == s2[j - 1])
				match = 2;
			else
				match = -3;
			if (score[i - 1][j - 1] + match == score[i][j]) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = s2[j - 1]; 
				i--; j--; 
			} 
			else if (score[i - 1][j] + gap == score[i][j]) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = '-'; 
				i--; 
			} 
			else if (score[i][j - 1] + gap == score[i][j]) { 
				s1ans[xpos--] = '-'; 
				s2ans[ypos--] = s2[j - 1]; 
				j--; 
			} 
		}
		while (xpos > 0) { 
			// if (i > 0) s1ans[xpos--] = s1[--i]; 
			/*else*/ s1ans[xpos--] = '-'; 
		} 
		while (ypos > 0) { 
			// if (j > 0) s2ans[ypos--] = s2[--j]; 
			/*else*/ s2ans[ypos--] = '-'; 
		}

		int skip = 1;
		for(int k = 1; k <= l; k++) {
			if(s1ans[k] == '-' && s2ans[k] == '-')
				skip++;
			else
				break;
		}

		// for(int k = skip; k <= l; k++)
		// 	printf("%c", s1ans[k]);
		// printf("\n");
		// for(int k = skip; k <= l; k++)
		// 	printf("%c", s2ans[k]);
		// printf("\n");

		std::string outbwt;
        outbwt = "local-naiveGap.aln";

        std::ofstream myfile(outbwt);
        if (myfile.is_open()) {
        	myfile << "Score = ";
        	myfile << score[indexi][indexj];
        	myfile << "\n\n";

        	int i = skip;
        	int j = skip;
        	int count = 0;
        	while(true) {
        		int now = i;
        		myfile << s1name;
        		myfile << " ";
        		for(; i <= l && i < (now + 60); i++)
        			myfile << s1ans[i];
        		myfile << "\n";
        		myfile << s2name;
        		myfile << " ";
        		for(; j <= l && j < (now + 60); j++)
        			myfile << s2ans[j];
        		myfile << "\n\n";
        		if(i > l)
        			break;
        	}
            myfile.close();
        }
        else
            std::cout << "Unable to open file";
		return 0;
	}
	else if (task.compare("aglobal") == 0){
		int gape = atoi(argv[8]);

		int ** e = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			e[i] = new int[n + 1];
		}

		int ** f = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			f[i] = new int[n + 1];
		}

		int ** g = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			g[i] = new int[n + 1];
		}

		int ** v = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			v[i] = new int[n + 1];
		}

		// backtrack
		int ** b = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			b[i] = new int[n + 1];
		}

		// E(i,j) is the max value of any alignment where s2j matches a space
		e[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			e[0][j] = 0;
		}
		for(int i = 1; i <= m; i++) {
			e[i][0] = i * gape + gap;
		}

		// F(i,j) is the max value of any alignment where s1i matches a space
		f[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			f[0][j] = j * gape + gap;
		}
		for(int i = 1; i <= m; i++) {
			f[i][0] = 0;
		}

		// G(i,j) is the max value of any alignment where s1i and s2j match (or mismatch)
		g[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			g[0][j] = 0;
		}
		for(int i = 1; i <= m; i++) {
			g[i][0] = 0;
		}

		// V(i,j) is the result
		v[0][0] = 0;
		for(int j = 1; j <= n; j++) {
			v[0][j] = j * gape + gap;
		}
		for(int i = 1; i <= m; i++) {
			v[i][0] = i * gape + gap;
		}

		// B(i,j) is the arrow showing the way diag: 0, up: 1, left: -1
		b[0][0] = -1;
		for(int j = 1; j <= n; j++) {
			b[0][j] = -1;
		}
		for(int i = 1; i <= m; i++) {
			b[i][0] = 1;
		}
		
		int match;

		for(int i = 1; i <= m; i++) {
			for(int j = 1; j <= n; j++) {
				if(s1[i - 1] == s2[j - 1]) {
					match = 2;
				}
				else
					match = -3;
				// E
				e[i][j] = fmax(e[i][j - 1] + gape, v[i][j - 1] + gap + gape);

				// F
				f[i][j] = fmax(f[i - 1][j] + gape, v[i - 1][j] + gap + gape);

				//G
				g[i][j] = v[i - 1][j - 1] + match;

				//V
				if(g[i][j] >= e[i][j] && g[i][j] >= f[i][j]) {
					v[i][j] = g[i][j];
					b[i][j] = 0;
				}
				else if(e[i][j] >= g[i][j] && e[i][j] >= f[i][j]) {
					v[i][j] = e[i][j];
					b[i][j] = -1;
				}
				else {
					v[i][j] = f[i][j];
					b[i][j] = 1;
				}
			}
		}
		// printf("score is %d\n", v[m][n]);

		// reconstructing the solution
		int l = m + n;
		int i = m;
		int j = n;

		int xpos = l; 
    	int ypos = l; 
  
		char * s1ans = new char[l + 1];
		char * s2ans = new char[l + 1];

		while (!(i == 0 || j == 0)) { 
			if (b[i][j] == 0) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = s2[j - 1]; 
				i--; j--; 
			} 
			else if (b[i][j] == 1) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = '-'; 
				i--; 
			} 
			else if (b[i][j] == -1) { 
				s1ans[xpos--] = '-'; 
				s2ans[ypos--] = s2[j - 1]; 
				j--; 
			} 
		}
		while (xpos > 0) { 
			if (i > 0) s1ans[xpos--] = s1[--i]; 
			else s1ans[xpos--] = '-'; 
		} 
		while (ypos > 0) { 
			if (j > 0) s2ans[ypos--] = s2[--j]; 
			else s2ans[ypos--] = '-'; 
		}

		int skip = 1;
		for(int k = 1; k <= l; k++) {
			if(s1ans[k] == '-' && s2ans[k] == '-')
				skip++;
			else
				break;
		}

		std::string outbwt;
        outbwt = "global-affineGap.aln";

        std::ofstream myfile(outbwt);
        if (myfile.is_open()) {
        	myfile << "Score = ";
        	myfile << v[m][n];
        	myfile << "\n\n";

        	int i = skip;
        	int j = skip;
        	int count = 0;
        	while(true) {
        		int now = i;
        		myfile << s1name;
        		myfile << " ";
        		for(; i <= l && i < (now + 60); i++)
        			myfile << s1ans[i];
        		myfile << "\n";
        		myfile << s2name;
        		myfile << " ";
        		for(; j <= l && j < (now + 60); j++)
        			myfile << s2ans[j];
        		myfile << "\n\n";
        		if(i > l)
        			break;
        	}
            myfile.close();
        }
        else
            std::cout << "Unable to open file";
		return 0;
	}
	else if (task.compare("alocal") == 0){
		int gape = atoi(argv[8]);

		int ** e = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			e[i] = new int[n + 1];
		}

		int ** f = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			f[i] = new int[n + 1];
		}

		int ** g = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			g[i] = new int[n + 1];
		}

		int ** v = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			v[i] = new int[n + 1];
		}

		// backtrack
		int ** b = new int*[m + 1];
		for(int i = 0; i <= m; i++) {
			b[i] = new int[n + 1];
		}

		// E(i,j) is the max value of any alignment where s2j matches a space
		e[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			e[0][j] = -10000;
		}
		for(int i = 1; i <= m; i++) {
			e[i][0] = gape + gap;
		}

		// F(i,j) is the max value of any alignment where s1i matches a space
		f[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			f[0][j] = gape + gap;
		}
		for(int i = 1; i <= m; i++) {
			f[i][0] = -10000;
		}

		// G(i,j) is the max value of any alignment where s1i and s2j match (or mismatch)
		g[0][0] = -10000;
		for(int j = 1; j <= n; j++) {
			g[0][j] = -10000;
		}
		for(int i = 1; i <= m; i++) {
			g[i][0] = -10000;
		}

		// V(i,j) is the result
		v[0][0] = 0;
		for(int j = 1; j <= n; j++) {
			v[0][j] = 0;
		}
		for(int i = 1; i <= m; i++) {
			v[i][0] = 0;
		}

		// B(i,j) is the arrow showing the way diag: 0, up: 1, left: -1, stop: -2
		b[0][0] = -2;
		for(int j = 1; j <= n; j++) {
			b[0][j] = -2;
		}
		for(int i = 1; i <= m; i++) {
			b[i][0] = -2;
		}
		
		int match;
		int max = 0;
		int indexi = 0;
		int indexj = 0;

		for(int i = 1; i <= m; i++) {
			for(int j = 1; j <= n; j++) {
				if(s1[i - 1] == s2[j - 1])
					match = 2;
				else
					match = (-3);
				// E
				e[i][j] = fmax(e[i][j - 1] + gape, v[i][j - 1] + gap + gape);

				// F
				f[i][j] = fmax(f[i - 1][j] + gape, v[i - 1][j] + gap + gape);

				//G
				g[i][j] = v[i - 1][j - 1] + match;

				// //V
				// v[i][j] = fmax(g[i][j], e[i][j]);
				// v[i][j] = fmax(v[i][j], f[i][j]);
				// v[i][j] = fmax(v[i][j], 0);

				// if(i == 8 && j == 5)
				// 	printf("8, 5: g %d, e %d, f %d \n", g[i][j], e[i][j], f[i][j]);

				// if(i == 7 && j == 5)
				// 	printf("7, 5: g %d, e %d, f %d \n", g[i][j], e[i][j], f[i][j]);

				//V
				if(0 >= g[i][j] && 0 >= e[i][j] && 0 >= f[i][j]) {
					v[i][j] = 0;
					b[i][j] = -2;
				}
				else if(g[i][j] >= e[i][j] && g[i][j] >= f[i][j] && g[i][j] >= 0) {
					v[i][j] = g[i][j];
					b[i][j] = 0;
				}
				else if(e[i][j] >= g[i][j] && e[i][j] >= f[i][j] && e[i][j] >= 0) {
					v[i][j] = e[i][j];
					b[i][j] = -1;
				}
				else /*if(f[i][j] >= e[i][j] && f[i][j] >= g[i][j] && f[i][j] >= 0)*/{
					v[i][j] = f[i][j];
					b[i][j] = 1;
				}
				// else {
				// 	v[i][j] = 0;
				// 	b[i][j] = -2;
				// }

				//update max
				if(v[i][j] > max) {
					max = v[i][j];
					indexi = i;
					indexj = j;
				}
			}
		}
		// printf("score is %d\n", v[indexi][indexj]);

		// reconstructing the solution
		int l = indexi + indexj;
		int i = indexi;
		int j = indexj;

		int xpos = l; 
    	int ypos = l; 
  
		char * s1ans = new char[l + 1];
		char * s2ans = new char[l + 1];

		while (!(i == 0 || j == 0 || b[i][j] == -2)) {
			// printf("b[%d][%d] is %d\n", i, j, b[i][j]);
			if (b[i][j] == 0) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = s2[j - 1]; 
				i--;
				j--; 
			} 
			else if (b[i][j] == 1) { 
				s1ans[xpos--] = s1[i - 1]; 
				s2ans[ypos--] = '-'; 
				i--; 
			} 
			else if (b[i][j] == -1) { 
				s1ans[xpos--] = '-'; 
				s2ans[ypos--] = s2[j - 1]; 
				j--; 
			}
		}
		while (xpos > 0) { 
			//if (i > 0) s1ans[xpos--] = s1[--i]; 
			/*else*/ s1ans[xpos--] = '-'; 
		} 
		while (ypos > 0) { 
			//if (j > 0) s2ans[ypos--] = s2[--j]; 
			/*else*/ s2ans[ypos--] = '-'; 
		}

		int skip = 1;
		for(int k = 1; k <= l; k++) {
			if(s1ans[k] == '-' && s2ans[k] == '-')
				skip++;
			else
				break;
		}

		std::string outbwt;
        outbwt = "local-affineGap.aln";

        std::ofstream myfile(outbwt);
        if (myfile.is_open()) {
        	myfile << "Score = ";
        	myfile << v[indexi][indexj];
        	myfile << "\n\n";

        	int i = skip;
        	int j = skip;
        	int count = 0;
        	while(true) {
        		int now = i;
        		myfile << s1name;
        		myfile << " ";
        		for(; i <= l && i < (now + 60); i++)
        			myfile << s1ans[i];
        		myfile << "\n";
        		myfile << s2name;
        		myfile << " ";
        		for(; j <= l && j < (now + 60); j++)
        			myfile << s2ans[j];
        		myfile << "\n\n";
        		if(i > l)
        			break;
        	}
            myfile.close();
        }
        else
            std::cout << "Unable to open file";
	}
	else {
		printf("Wrong mode. ERROR\n");
		return -1;
	}
}