/**
 * @file buildUPGMA.cpp
 * @brief Guide Tree Construction
 * @author Naisila Puka
 * @version 21/12/2019
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

// finds minimum distance in the distance matrix and its location
float minDistance(std::vector<std::vector<float> >& matrix, int& x, int& y) {
    float minDistance = matrix[1][0];
    x = 1;
    y = 0;

    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] < minDistance) {
                minDistance = matrix[i][j];
                x = i;
                y = j;
            }
        }
    }

    // Uncomment next line to track the minimum distance
    // printf("min is %f\n", minDistance);
    return minDistance;
}

// updates the netwick format "tree"
void netwickElement(std::vector<std::string>& labels, int x, int y, float distance, std::map<std::string, float>& treeHeights, std::map<std::string, int>& elemsNo) {
    if (y < x) {
        int temp = x;
        x = y;
        y = temp;
    }
    // (A:4.5, B:4.5)
    std::string newlabel = "(";
    newlabel += labels[x];
    newlabel += ":";
    newlabel += std::to_string((distance / 2) - treeHeights[labels[x]]);
    newlabel += ", ";
    newlabel += labels[y];
    newlabel += ":";
    newlabel += std::to_string((distance / 2) - treeHeights[labels[y]]);
    newlabel += ")";

    // updates the number of elements in the current label
    // eg, if we are merging A and B, (A:4.5, B:4.5), elemsNo becomes 2
    elemsNo[newlabel] = elemsNo[labels[x]] + elemsNo[labels[y]];
    labels[x] = newlabel;
    treeHeights[labels[x]] = distance / 2;
    labels.erase(labels.begin() + y);
}

// updates the matrix after merging the elements whose distance is minimum
void updateMatrix(std::vector<std::vector<float> >& matrix, int x, int y, std::vector<std::string>& labels, std::map<std::string, int>& elemsNo) {
    if (y < x) {
        int temp = x;
        x = y;
        y = temp;
    }

    // Uncomment the followinf printf lines to track the matrix update

    // For the lower index x, reconstruct the row (x, i), where i < x
    std::vector<float> newRow;
    for (int i = 0; i < x; i++) {
        // printf("%f, %f \n", matrix[x][i], matrix[y][i]);
        newRow.push_back((matrix[x][i] * elemsNo[labels[x]] + matrix[y][i]* elemsNo[labels[y]]) / (elemsNo[labels[x]] + elemsNo[labels[y]]));
        // printf("%f \n", newRow[i]);
    }
    matrix[x] = newRow;

    // Reconstruct the column (i, x), where y > i > x
    // Since the matrix is lower triangular, row y only contains values for indices < y
    for (int i = x + 1; i < y; i++) {
        // printf("%f, %f \n", matrix[i][x], matrix[y][i]);
        matrix[i][x] = (matrix[i][x] * elemsNo[labels[x]] + matrix[y][i] * elemsNo[labels[y]]) / (elemsNo[labels[x]] + elemsNo[labels[y]]);
        // printf("%f \n", matrix[i][x]);
    }

    for (int i = y + 1; i < matrix.size(); i++) {
        // printf("%f, %f \n", matrix[i][x], matrix[i][y]);
        matrix[i][x] = (matrix[i][x] * elemsNo[labels[x]] + matrix[i][y] * elemsNo[labels[y]]) / (elemsNo[labels[x]] + elemsNo[labels[y]]);
        // printf("%f \n", matrix[i][x]);
        matrix[i].erase(matrix[i].begin() + y);
    }

    matrix.erase(matrix.begin() + y);
}

// overall upgma function which utilizes all of the above
std::string upgma(std::vector<std::vector<float> >& matrix, std::vector<std::string>& labels, std::map<std::string, float>& treeHeights, std::map<std::string, int>& elemsNo) {
    int x, y;
    float distance;
    while (labels.size() > 1) {
        distance = minDistance(matrix, x, y);
        updateMatrix(matrix, x, y, labels, elemsNo);
        netwickElement(labels, x, y, distance, treeHeights, elemsNo);
    }
    // netwick format tree result will be in labels[0] because we always delete the maximum index element when merging
    // so 0 is never deleted because its the minimum index
    return labels[0];
}

// calculates edit distance of two sequences by using the needleman-wunsch alignment
// code utilized from previous homework
float distance(const std::string& s1, const std::string& s2, int match, int mismatch, int gap, int gape, int m, int n) {
    int** e = new int*[m + 1];
    for (int i = 0; i <= m; i++)
        e[i] = new int[n + 1];

    int** f = new int*[m + 1];
    for (int i = 0; i <= m; i++)
        f[i] = new int[n + 1];

    int** g = new int*[m + 1];
    for (int i = 0; i <= m; i++)
        g[i] = new int[n + 1];

    int** v = new int*[m + 1];
    for (int i = 0; i <= m; i++)
        v[i] = new int[n + 1];

    // backtrack
    int** b = new int*[m + 1];
    for (int i = 0; i <= m; i++)
        b[i] = new int[n + 1];

    // E(i,j) is the max value of any alignment where s2j matches a space
    e[0][0] = -10000;
    for (int j = 1; j <= n; j++)
        e[0][j] = 0;
    for (int i = 1; i <= m; i++)
        e[i][0] = i * gape + gap;

    // F(i,j) is the max value of any alignment where s1i matches a space
    f[0][0] = -10000;
    for (int j = 1; j <= n; j++)
        f[0][j] = j * gape + gap;
    for (int i = 1; i <= m; i++)
        f[i][0] = 0;

    // G(i,j) is the max value of any alignment where s1i and s2j match (or mismatch)
    g[0][0] = -10000;
    for (int j = 1; j <= n; j++)
        g[0][j] = 0;
    for (int i = 1; i <= m; i++)
        g[i][0] = 0;

    // V(i,j) is the result
    v[0][0] = 0;
    for (int j = 1; j <= n; j++)
        v[0][j] = j * gape + gap;
    for (int i = 1; i <= m; i++)
        v[i][0] = i * gape + gap;

    // B(i,j) is the arrow showing the way diag: 0, up: 1, left: -1
    b[0][0] = -1;
    for (int j = 1; j <= n; j++)
        b[0][j] = -1;
    for (int i = 1; i <= m; i++)
        b[i][0] = 1;

    int match1;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1])
                match1 = match;
            else
                match1 = mismatch;
            // E
            e[i][j] = fmax(e[i][j - 1] + gape, v[i][j - 1] + gap + gape);

            // F
            f[i][j] = fmax(f[i - 1][j] + gape, v[i - 1][j] + gap + gape);

            //G
            g[i][j] = v[i - 1][j - 1] + match1;

            //V
            if (g[i][j] >= e[i][j] && g[i][j] >= f[i][j]) {
                v[i][j] = g[i][j];
                b[i][j] = 0;
            }
            else if (e[i][j] >= g[i][j] && e[i][j] >= f[i][j]) {
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

    char* s1ans = new char[l + 1];
    char* s2ans = new char[l + 1];

    while (!(i == 0 || j == 0)) {
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
        if (i > 0)
            s1ans[xpos--] = s1[--i];
        else
            s1ans[xpos--] = '-';
    }
    while (ypos > 0) {
        if (j > 0)
            s2ans[ypos--] = s2[--j];
        else
            s2ans[ypos--] = '-';
    }

    int skip = 1;
    for (int k = 1; k <= l; k++) {
        if (s1ans[k] == '-' && s2ans[k] == '-')
            skip++;
        else
            break;
    }

    i = skip;
    int result = 0;
    for (; i <= l; i++) {
        if (s1ans[i] != s2ans[i])
            result++;
    }

    // Uncomment the following lines to see the alignments

    // printf("\n");
    //     i = skip;
    //     for(; i <= l; i++) {
    //     	printf("%c", s1ans[i]);
    //     }
    //     printf("\n");
    //     i = skip;
    //     for(; i <= l; i++) {
    //     	printf("%c", s2ans[i]);
    //     }
    //     printf("\n");

    for (int i = 0; i <= m; i++) {
        delete[] e[i];
        delete[] f[i];
        delete[] g[i];
        delete[] v[i];
        delete[] b[i];
    }

    delete[] e;
    delete[] f;
    delete[] g;
    delete[] v;
    delete[] b;

    delete[] s1ans;
    delete[] s2ans;

    return float(result);
}

int main(int argc, char** argv) {
    struct option opts[] = {
        { "fasta", required_argument, 0, 0 },
        { "match", required_argument, 0, 1 },
        { "mismatch", required_argument, 0, 2 },
        { "gapopen", required_argument, 0, 3 },
        { "gapext", required_argument, 0, 4 },
        { "out", required_argument, 0, 5 },
        { 0, 0, 0, 0 }
    };

    int c;
    int option_index = 0, all = 0;
    char *fastaFileName = NULL, *outFileName = NULL;
    int match, mismatch, gapopen, gapext;

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
            match = atoi(optarg);
            all++;
            break;
        case 2:
            mismatch = atoi(optarg);
            all++;
            break;
        case 3:
            gapopen = atoi(optarg);
            all++;
            break;
        case 4:
            gapext = atoi(optarg);
            all++;
            break;
        case 5:
            outFileName = optarg;
            all++;
            break;
        }
    }

    if (all < 5) {
        printf("Error! Not all arguments are provided!\n");
        return -1;
    }

    std::string line;

    std::string seqs[25];
    std::string names[25];
    //n: number of sequences, 1<=n<=25
    int n = -1;

    //get sequences
    std::ifstream input_fasta(fastaFileName);
    while (std::getline(input_fasta, line)) {
        if (line.empty() || line[0] == '>') {
            if (line[0] == '>' && line.size() > 1) {
                n++;
                names[n] += line.substr(1, line.size() - 1);
            }
        }
        else
            seqs[n] += line;
    }
    n++;

    std::vector<std::vector<float> > matrix(n);
    for (int i = 0; i < n; i++)
        matrix[i].resize(i);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            matrix[i][j] = distance(seqs[i], seqs[j], match, mismatch, gapopen, gapext, seqs[i].size(), seqs[j].size());

    // Uncomment the following lines to try a toy example
    
    // std::vector<std::vector<float> > matrix(6);
    // for (int i = 0; i < 6; i++)
    //     matrix[i].resize(i);

    // matrix[1][0] = 19;
    // matrix[2][0] = 27;
    // matrix[2][1] = 31;
    // matrix[3][0] = 8;
    // matrix[3][1] = 18;
    // matrix[3][2] = 26;
    // matrix[4][0] = 33;
    // matrix[4][1] = 36;
    // matrix[4][2] = 41;
    // matrix[4][3] = 31;
    // matrix[5][0] = 18;
    // matrix[5][1] = 1;
    // matrix[5][2] = 32;
    // matrix[5][3] = 17;
    // matrix[5][4] = 35;


    // int x, y, a;
    // a = minDistance(matrix, x, y);
    // printf("min: %d, x: %d, y: %d\n", a, x, y);

    // for (int i = 0; i < matrix.size(); i++) {
    //     for (int j = 0; j < matrix[i].size(); j++) {
    //         std::cout << matrix[i][j] << " ";
    //     }
    //     printf("\n");
    // }

    std::map<std::string, float> treeHeights;
    std::map<std::string, int> elemsNo;
    std::vector<std::string> labels;
    for (int i = 0; i < n; i++) {
        labels.push_back(names[i]);
        treeHeights[labels[i]] = 0;
        elemsNo[labels[i]] = 1;
    }

    // std::string name[6];
    // name[0]="A";
    // name[1]="B";
    // name[2]="C";
    // name[3]="D";
    // name[4]="E";
    // name[5]="F";

    // for (int i = 0; i < 6; i++) {
    //     labels.push_back(name[i]);
    //     treeHeights[name[i]] = 0;
    //     elemsNo[name[i]] = 1;
    // }

    std::string result = upgma(matrix, labels, treeHeights, elemsNo);

    std::ofstream myfile(outFileName);
	if (myfile.is_open()) {
		myfile << result;
		myfile.close();
	}
	else
		std::cout << "Unable to open file";

    return 0;
}