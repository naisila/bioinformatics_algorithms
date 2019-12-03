/**
 * @file bwtfm.cpp
 * @brief This program implements Burrows-Wheeler Transformation (BWT)
 *        and Ferragina-Manzini (FM) index.
 * @author Naisila Puka
 * @version 30/10/2019
*/

#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include <chrono>
std::string s;

int cmp(int first, int second) {
	return strcmp(s.substr(first - 1, s.size() - first + 1).c_str(), s.substr(second - 1, s.size() - second + 1).c_str()) < 0? 1 : 0; 
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printf("MOT ENOUGH ARGUMENTS PROVIDED\n");
        return 0;
    }
    std::string task;
    task += argv[1];

    if (task.compare("search") == 0) {
        std::string text_name, pattern_name;
        std::string bwtfile, bwt, fmfile, pattern, line;

        //get text name
        std::ifstream input_text(argv[2]);
        while (std::getline(input_text, line)) {
            if (line.empty() || line[0] == '>') {
                if (line[0] == '>' && line.size() > 1) {
                    text_name += line.substr(1, line.size() - 1);
                    break;
                }
                continue;
            }
        }

        //read bwt string
        bwtfile += argv[2];
        bwtfile += ".bwt";
        std::ifstream input_bwt(bwtfile);
        while (std::getline(input_bwt, line)) {
            if (line.empty() || line[0] == '>')
                continue;
            else
                bwt += line;
        }
        int m = bwt.size();

        //read ranks, sa, occa, occc, occg, occt
        fmfile += argv[2];
        fmfile += ".fm";
        std::ifstream input_text1(fmfile);

        //ranks
        std::map<char, int> rank;
        std::getline(input_text1, line);
        std::stringstream lineStream(line);
        int value;
        lineStream >> value;
        rank['A'] = value;
        lineStream >> value;
        rank['C'] = value;
        lineStream >> value;
        rank['G'] = value;
        lineStream >> value;
        rank['T'] = value;
        lineStream >> value;

        //sa
        int* sa = new int[m];
        std::getline(input_text1, line);
        std::stringstream lineStream1(line);
        int index = 0;
        while (lineStream1 >> value) {
            sa[index] = value;
            index++;
        }

        //occ
        std::map<char, int*> occ;
        occ['A'] = new int[m];
        occ['C'] = new int[m];
        occ['G'] = new int[m];
        occ['T'] = new int[m];
        //occa
        std::getline(input_text1, line);
        std::stringstream lineStream2(line);
        index = 0;
        while (lineStream2 >> value) {
            occ['A'][index] = value;
            index++;
        }
        //occc
        std::getline(input_text1, line);
        std::stringstream lineStream3(line);
        index = 0;
        while (lineStream3 >> value) {
            occ['C'][index] = value;
            index++;
        }
        //occg
        std::getline(input_text1, line);
        std::stringstream lineStream4(line);
        index = 0;
        while (lineStream4 >> value) {
            occ['G'][index] = value;
            index++;
        }
        //occt
        std::getline(input_text1, line);
        std::stringstream lineStream5(line);
        index = 0;
        while (lineStream5 >> value) {
            occ['T'][index] = value;
            index++;
        }

        //pattern
        std::ifstream input_pattern(argv[3]);
        while (std::getline(input_pattern, line)) {
            if (line.empty() || line[0] == '>') {
                if (line[0] == '>' && line.size() > 1)
                    pattern_name += line.substr(1, line.size() - 1);
                continue;
            }
            else
                pattern += line;
        }
        int n = pattern.size();

        // bwt searching
        int i = n - 1;
        char next = pattern[i--];
        int start = rank[next] + occ[next][0];
        int end = rank[next] + (occ[next][m - 1] - 1);
        int startrank, endrank;

        auto t1 = std::chrono::high_resolution_clock::now();
        while (start > 0 && end > 0 && i >= 0) {
            next = pattern[i--];
            startrank = rank[next] + occ[next][start - 1];
            endrank = rank[next] + occ[next][end] - 1;
            start = startrank;
            end = endrank;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "\nPattern search done in " << runtime << " microseconds" << std::endl;
        if (start > 0 && end > 0 && end >= start) {
            printf("Pattern %s found in text %s in position(s): ", pattern_name.c_str(), text_name.c_str());
            for (int i = start; i < end + 1; i++) {

                std::cout << sa[i] << " ";
            }
            printf("\n");
        }
        else
            printf("Pattern not found");

        return 0;
    }
    else if (task.compare("index") == 0) {
        std::ifstream input_text(argv[2]);
        std::string line;

        while (std::getline(input_text, line)) {
            if (line.empty() || line[0] == '>')
                continue;
            else
                s += line;
        }
        s += "$";
        int m = s.size();

		int* sa = new int[m];

		auto t1 = std::chrono::high_resolution_clock::now();  
        for (int i = 0; i < m; i++) 
        	sa[i]= i + 1; 
        std::sort(sa, sa + m, cmp);

        char* bwt = new char[m];
        for (int i = 0; i < m; i++) {
            if (sa[i] - 1 > 0)
                bwt[i] = s[sa[i] - 2];
            else
                bwt[i] = '$';
        }

        std::map<char, int*> occ;

        occ['A'] = new int[m];
        occ['A'][0] = 0;

        occ['C'] = new int[m];
        occ['C'][0] = 0;

        occ['G'] = new int[m];
        occ['G'][0] = 0;

        occ['T'] = new int[m];
        occ['T'][0] = 0;

        occ[bwt[0]][0] = 1;

        for (int i = 1; i < m; i++) {
            occ['A'][i] = occ['A'][i - 1];
            occ['C'][i] = occ['C'][i - 1];
            occ['G'][i] = occ['G'][i - 1];
            occ['T'][i] = occ['T'][i - 1];
            if (bwt[i] != '$')
                occ[bwt[i]][i]++;
        }

        std::map<char, int> rank;

        rank['A'] = (occ['A'][m - 1] > 0) ? 1 : -1; //A
        rank['C'] = rank['A'] + occ['A'][m - 1]; //C
        rank['G'] = rank['C'] + occ['C'][m - 1]; //G
        rank['T'] = rank['G'] + occ['G'][m - 1]; //T
        auto t2 = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "Index generation done in " << runtime << " microseconds" << std::endl;

        std::string outbwt;
        outbwt += argv[2];
        outbwt += ".bwt";

        std::ofstream myfile(outbwt);
        if (myfile.is_open()) {
            for (int i = 0; i < m; i++) {
                myfile << bwt[i];
                if ((i + 1) % 60 == 0)
                    myfile << "\n";
            }
            myfile.close();
        }
        else
            std::cout << "Unable to open file";

        std::string outfm;
        outfm += argv[2];
        outfm += ".fm";
        std::ofstream myfile1(outfm);
        if (myfile1.is_open()) {
            //ranks
            myfile1 << rank['A'] << " " << rank['C'] << " " << rank['G'] << " " << rank['T'];
            myfile1 << "\n";
            //sa
            for (int i = 0; i < m; i++) {
                myfile1 << sa[i] << " ";
            }
            myfile1 << "\n";
            //occ[A]
            for (int i = 0; i < m; i++) {
                myfile1 << occ['A'][i] << " ";
            }
            myfile1 << "\n";
            //occ[C]
            for (int i = 0; i < m; i++) {
                myfile1 << occ['C'][i] << " ";
            }
            myfile1 << "\n";
            //occ[G]
            for (int i = 0; i < m; i++) {
                myfile1 << occ['G'][i] << " ";
            }
            myfile1 << "\n";
            //occ[T]
            for (int i = 0; i < m; i++) {
                myfile1 << occ['T'][i] << " ";
            }
        }
        else
            std::cout << "Unable to open file";
        return 0;
    }
    else {
    	printf("Please provide the right parameters. Refer to README.txt");
        return 0;
    }
}