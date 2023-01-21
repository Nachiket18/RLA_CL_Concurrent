#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <fstream>
#include <dirent.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <bits/stdc++.h>
#include <omp.h>
#include <time.h>
#include <algorithm>
#include <execution>

#include <mutex>
#include <shared_mutex>

#include "RLA_parallel.h"


int numThreads = 6;
int lmerUsed = 3;
int lenMax;
static const int ALPHABET_SIZE_LIST[] = {26, 10, 36, 256};
int totalRecords;
int totalUniqueRecords = 0;



std::vector<std::vector<int> > blockArray;
std::vector<std::pair<int, int> > edgeList;
map<int, vector<int> > approxConnectedComponents;
std::vector<std::pair<int, string> > combinedData;
vector<pair<int,string> > uniqueRecords;
std::vector<std::vector<string>>  vec2D;
std::vector<std::vector<string>>  vec2D_original;
std::vector<std::vector<int> > finalConnectedComponents;
std::vector<int> x_prime_index;
std::vector<std::vector<int> > exactmatches;
std::vector<int> x_prime_ref ;	


// helps edit distance calculation in calculateBasicED()

int calculateBasicED2(string& str1, string& str2, int threshRem)
{
	int row, col, i, j;
	vector < vector < int > > matArr;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;

	matArr.resize(row);
	for(i = 0; i < row; ++i)
		matArr[i].resize(col, 0);

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if((int)str1[i-1] == (int)str2[j-1])
					matArr[i][j]	= min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] 	= min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}

	return (matArr[row-1][col-1]);
}

// calculates edit distance between two string
// takes two strings and a threshold value as input
// returns global variable threshold + 1 if distance exceeds theshold param
// returns edit distance
// core mechanism is a DP algo 

int calculateBasicED(string& str1, string& str2, int threshRem)
{
	int dist	= threshRem;

	if(abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if(str1.compare(str2) == 0)
		return 0;
	else if(dist == 0)
		return threshold + 1;
	else if((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{
		string s1, s2;
		int row, col, diagonal;
		int i, j;
		vector<vector<int> > matArr;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row	 		= s1.length() + 1;
		col 		= 2 * dist + 1;
		diagonal 	= dist + s2.length() - s1.length();

		matArr.resize(row);
		for(i = 0; i < row; ++i)
			matArr[i].resize(col, 0);


		for(i = 0; i < dist + 1; i++)
		{
			for(j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j]	= j - dist;
				else if(j == (dist - i))
					matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for(i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for(j = 0; j < col; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for(i = s2.length() - dist + 1; i < row; i++)
		{
			for(j = 0; j < col - i + s2.length() - dist; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		return matArr[row - 1][diagonal];

	}
}

int power_mod(int base,int pow,int p){
	int res = 0;
	for (int i=2; i <= pow;i++){
	
		res += ((base * base) % p);	
	}
	return res % p;
}

std::vector<int> generateIntKmersMap(string& str,int k){

	std::vector<int> k_mer_list;
	int res = 0;
	int p = 97;
	/*
		The Integer value for first K- Mer
	*/
	
	for (int i= 0;i < k;i++){
		res += (int)str[i] * power_mod(26,(k-i+1),p);	
	}
	
	k_mer_list.push_back(res);

	int tmp_pw = pow(26,k-1);	
	if (str.length() >= k+1) {
		
		for (int i=1; i <= (str.length()- k);i++){
			k_mer_list.push_back(((res - ((int)str[i-1] * tmp_pw) * 26 % p) + (int)str[(i-1) + k] ));
		}
	}
	
	return k_mer_list;
}

double calculateBasicQgramIntMap(string& str1, string& str2,int threshRem ){
	
	int k = 3;
	int min_length_str = (str1.size() - str2.size());

	if (min_length_str < 0) {
		min_length_str = - min_length_str;
	}

	if (min_length_str > 1) {
		return 0.0;
	}	
	
	
	if (str1.length() < k) {
		str1 = str1 + str1;
	}
	if (str2.length() < k) {
		str2 = str2 + str2;
	}
	

	
	std::vector< int > s1,s2;
	
	s1 = generateIntKmersMap(str1,k);
	s2 = generateIntKmersMap(str2,k);
	

	//int min_length = min(s1.size(),s2.size());
	int match = 0;

	for(int i = 0; i < s1.size(); ++i) {
    		
    		if ((s1[i] == s2[i]) || (s1[i] == s2[i+1]) || s1[i] == s2[i-1] ){
    			match += 1;
    		}
    
	}
      	
	
	
	int union_set_size = s1.size() + s2.size() - match;

	return (double)match / (double)union_set_size;
	
}


/*
	Reading the data from set of files and merging together into a data vector.
*/

vector<vector<string>> getData(vector<string> files) {
    string line;
    vector<vector<string>> vec2D;
	int lenMax = 0;
    for (int i=0; i < files.size();i++) {
        
        if (files[i] != "") {
            ifstream records(files[i]);
            int max = 0;
			while(getline(records,line)){
                vector<string> result;
                boost::split(result, line, boost::is_any_of("\t"));
                vector<string> vec;
                string name = result[1] + result[2];
                boost::to_lower(result[1]);
				boost::to_lower(result[2]);
				auto last = std::remove_if(result[1].begin(), result[1].end(), [](auto ch) {
        								return ::isdigit(ch) || ::ispunct(ch) || ::iswpunct(ch);
    								});
				result[1].erase(last, result[1].end());
                vec.push_back(result[0]);
                vec.push_back(result[1]);
                vec.push_back(result[2]);
				vec.push_back(result[3]);
				
                vec2D.push_back(vec);

            }

        records.close();
        }   
    }
    
    return vec2D;
}

vector<std::string> getInputFileList(string directory){
    
    DIR *dr;
    struct dirent *en;
    vector<std::string> fileList;

    dr = opendir(directory.c_str()); //open all directory
    if (dr) {
        while ((en = readdir(dr)) != NULL) {
            fileList.push_back(directory + "/" + (string)en->d_name);
        }
        closedir(dr); //close all directory
    }
   return fileList;

}


/*
	Initial Sorting of the Input data.
*/

void sort_data_parallel(vector<vector<string> > &data_vector,int column_1,int column_2){


	auto sortRuleLambda = [column_1,column_2] (std::vector<string> const& v1, std::vector<string> const& v2) -> bool
    {
       if (v1[column_1] != v2[column_1]) {

			return v1[column_1] < v2[column_1];
	   }

	   return v1[column_2] < v2[column_2];
    };

	std::sort(std::execution::par,data_vector.begin(), data_vector.end(),sortRuleLambda);

}




/*
	1. Find the unique elements in the Data Vector Array for specified range (start,limit)
	2. Used by multiple threads (on seperate start,limit values)
	3. Since reading can be done concurrently, we have used shared mutex,
*/


void find_group_representative(vector<vector<string>> &data_vector,int start,int limit,vector<string> &x_prime,vector<string> &temp,vector<int> &indexTemp){
	
	// use for capturing the index of the deduplicate data. Used further while blocking 

	temp.emplace_back(data_vector[start][1]);
	indexTemp.push_back(start);
	
	for(int i = start+1; i < limit; i++) {	

		if (data_vector[i][1] == data_vector[i-1][1]){
			continue;
		}
		else{

			temp.emplace_back(data_vector[i][1]);
			indexTemp.push_back(i);
		}
		
	}
	
} 


void getCombinedData() {
	string strSample(50, '0');
	combinedData.resize(totalRecords);
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		p.second = vec2D[i][1] + vec2D[i][2] + vec2D[i][3];
		combinedData[i]=p;
		if (max < p.second.size()) {
			max = p.second.size();
		}
	}
}

void perform_blocking(std::vector<string>  &x_prime_t,int lmerUsed,int type){

	if (lmerUsed < 3)
		lmerUsed = 3;
	
	int blockTotal 	= pow(26, lmerUsed);
	string strSample;
	blockArray.resize(blockTotal);
	int numElements = x_prime_t.size() / numThreads;
	std::vector<int> threadDivideData;
	int currentThreadData = 0;
	threadDivideData.push_back(0);

	for(int i = 0; i < x_prime_t.size(); i++){

		if (currentThreadData < numElements){
			currentThreadData += 1;
		}
		else{
			threadDivideData.push_back(i);
			currentThreadData = 0;
		}

	}
	threadDivideData.push_back(x_prime_t.size()-1);	
	
	#pragma omp parallel num_threads(numThreads)
	{	
		std::vector< std::vector<int> > temp_blocks_1;
		temp_blocks_1.resize(blockTotal*2);

		#pragma omp for ordered
		for(int i = 1; i < threadDivideData.size(); i++) {
			
			for (int m = threadDivideData[i-1]; m < threadDivideData[i]; m++){
				
				if (x_prime_t[m].size() < lmerUsed){
					int blockID=0;
					blockID += (x_prime_t[m][0] - 97) * (int) pow(ALPHABET_SIZE_LIST[type],lmerUsed - 1);
					if( !( blockID < 0 || blockID >= blockTotal ) ) {
						temp_blocks_1[blockID].push_back(m);	 
					}
				}
				else{
					for (int j = 0; j < (x_prime_t[m].size() - lmerUsed + 1); j++){
						int blockID=0;
						for (int k = 0; k < lmerUsed; ++k)
							if (type == 0)
								blockID	+= (x_prime_t[m][j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type],lmerUsed - k - 1);
							else if(type == 1)
								blockID	+= (x_prime_t[m][j + k] - 48) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
							else if(type == 2)
							{
								if(x_prime_t[m][j + k] >= 97)
									blockID	+= (x_prime_t[m][j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
								else
									blockID	+= (x_prime_t[m][j + k] - 22) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1); // 48 - 26
							}
						if( !( blockID < 0 || blockID >= blockTotal ) ) {
							temp_blocks_1[blockID].push_back(m);	 
						}
					}
				}	
					
			}
			#pragma omp ordered 
			{
				for (int k = 0; k < temp_blocks_1.size();k++){
					for (int j = 0; j < temp_blocks_1[k].size(); j++){
						blockArray[k].push_back(temp_blocks_1[k][j]);
					}
				}
			}			
		}
	}
}


/*
	Performing the de-duplication of records pre blocking serially.
	Was used to check correctness of the parallel implementation and difference in serial and parallel run time. 
*/

void find_group_representative_serial(vector<vector<string> > &data_vector,int start,int limit, vector<string> &x_prime){
	
	
	x_prime.push_back(data_vector[start][1]);
	
	for(int i = start+1; i < limit; i++) {	

		if (data_vector[i] == data_vector[i-1]){
			continue;
		}
		else{
			x_prime.push_back(data_vector[i][1]);
		}
		
	}

	int totalS = 0;	
	cout << "total size" << x_prime.size() << endl;

}


/*
	OpenMP program to perform concurrent de-duplication.

*/

void find_group_representative_openmp(vector<vector<string>> &data_vector,int start,int limit, vector<string> &x_prime, vector<int> &x_prime_index){
	

		
	#pragma omp parallel num_threads(numThreads)
	{

		vector<string> temp;
		vector<int> indexTemp;
		
		temp.push_back(data_vector[start][1]);
		indexTemp.push_back(start);
	
		#pragma omp for ordered
		for(int i = start+1; i < limit; i++) {	

			if (data_vector[i][1] != data_vector[i-1][1]){
				temp.push_back(data_vector[i][1]);
				indexTemp.push_back(i);
			}
			
		}
		#pragma omp critical
		{
			x_prime.insert(x_prime.end(), temp.begin(), temp.end());
			x_prime_index.insert(x_prime_index.end(), indexTemp.begin(), indexTemp.end());
		}
		
	}

	int totalS = 0;	
	cout << "total size of de-duplication vector" << x_prime.size() << endl;

}

void generate_x_prime_reprentative(vector< vector<string> > &data_vector,vector<string> &x_prime, vector<int> &x_prime_index){

	std::vector<std::thread> threads;
	
	int numElements = data_vector.size() / numThreads;
	int remainingElements = data_vector.size();
	int currentRecordPointer = 0;
	std::vector<int> threadDivideData;
	threadDivideData.push_back(0);
	//x_prime_index.resize(data_vector.size()/4);

	for (int i=0; i < numThreads; i++){

		int start = 0;
		int end = 0;
		std::vector<std::string> out;
		if ( (int) data_vector.size() < (currentRecordPointer + numElements) ){
			
			start = currentRecordPointer;
			end = data_vector.size();
			currentRecordPointer = end;
			remainingElements = 0;
			
		}
		else {
			int idealRange = (currentRecordPointer) + numElements;
			int tmp = (currentRecordPointer) + (numElements+1);
			
			while(data_vector[tmp][1] == data_vector[idealRange][1]){
				tmp +=1;
			}

			start = currentRecordPointer;
			end = tmp;
			
			currentRecordPointer = tmp;
			remainingElements = data_vector.size() - currentRecordPointer;			
			
		}
		threadDivideData.push_back(end);
	
	}

	
	#pragma omp parallel num_threads(numThreads)
	{
		vector<string> temp;
		vector<int> indexTemp;	
		vector< vector<int> > localMatches;
		#pragma omp for ordered
		for (int i = 1; i < threadDivideData.size(); i++){
			
			temp.push_back(data_vector[threadDivideData[i-1]][1]);
			indexTemp.push_back(threadDivideData[i-1]);

			for(int j = threadDivideData[i-1]+1; j < threadDivideData[i]; j++) {	

				if (data_vector[j] == data_vector[j-1] ){
					indexTemp.push_back(j);
				}
				else{
					temp.push_back(data_vector[j][1]);
					localMatches.push_back(indexTemp);
					indexTemp.clear();
					indexTemp.push_back(j);
				}
				
			}
			localMatches.push_back(indexTemp);

			#pragma omp ordered
			{
				int size_x_prime = x_prime_index.size();
				x_prime.insert(x_prime.end(), temp.begin(), temp.end());	
				exactmatches.insert(exactmatches.end(),localMatches.begin(),localMatches.end());
				
				/* 
					DEBUGGING PURPOSES
					//cout << "exatcMatches" << exactmatches.size() << endl;
					//cout << "X-prime,Ref" << x_prime.size() << "," << x_prime_ref.size() << endl;
				*/
				
				
			}
		
		}

	}
	
}



void generate_edges_concurrent(std::vector<std::vector<string>> &data_vector , std::vector<string>  &x_prime_t,std::vector<int> &x_prime_index,int dataSize){

	int currentThreadData = 0;

	#pragma omp parallel num_threads(numThreads) 
	{
		std::set<pair<int, int> > set_of_edges;

		#pragma omp for schedule(dynamic,1) nowait
		for( int m = 0; m < blockArray.size(); m++) {
				
				
				std::vector<int> block = blockArray[m];
				std::vector<std::vector<std::string>> dataArr(block.size());
				int l = 0;
				int set = 0;
				BOOST_FOREACH(int p, block) {

					dataArr[l] = vec2D[uniqueRecords[p].first];
					l++;
					
    			};
				
				//std::sort(block.begin(),block.end());
				//block.erase(unique( block.begin(), block.end() ), block.end());

				if (block.size() == 0) {
					continue;
				}

				
				for (int j = 0; j < block.size() - 1; j++){
					for (int k = j+1; k < block.size();k++){
						int first_j = block[j];
						int first_k = block[k];
						std::pair<int, int> edge_pair;
						edge_pair.first = first_j;
						edge_pair.second = first_k;
						// if (set_of_edges.count(edge_pair) > 0){
						// 	continue;
						// }

						if (first_j > first_k){
							
							if (isLinkageOk(dataArr[j],dataArr[k],1)){
								set_of_edges.insert(edge_pair);								
							}
						}
						else{

							if (isLinkageOk(dataArr[j],dataArr[k],1)){
								set_of_edges.insert(edge_pair);
							}
						}

														
					}
			}		
			
		}
		
		#pragma omp critical
		{
				set<std::pair<int,int>>::iterator itr;
				edgeList.insert(edgeList.end(),set_of_edges.begin(), set_of_edges.end());
				/*
					DEBUGGING PURPOSES
					cout << set_of_edges.size() << endl;
				*/
				

		}	

	}
	
}

void generate_edges_serial(std::vector<std::vector<string>> &data_vector , std::vector<string>  &x_prime_t,std::vector<int> &x_prime_index,int dataSize) {

	std::set<pair<int, int> > set_of_edges;
	for (int i = 0; i < blockArray.size(); i++){
		std::vector<int> block = blockArray[i];
		for (int j = 0; j < block.size(); j++){
					for (int k = j+1; k < block.size();k++){
						int first_j = block[j];
						int first_k = block[k];
						int dist_last_name = calculateBasicED(x_prime_t[first_j],x_prime_t[first_k],1);
						int dist_first_name = calculateBasicED(data_vector[x_prime_index[first_j]][2],data_vector[x_prime_index[first_k]][2],1);
						int dist = dist_first_name + dist_last_name;

						if (first_j > first_k){
							if (dist <= 1){
								std::pair<int, int> edge_pair;
								edge_pair.first = first_j;
								edge_pair.second = first_k;
								if (!set_of_edges.count(edge_pair)){
									set_of_edges.insert(edge_pair);
								}
							}
						}
						else{

							if (dist <= 1){
								std::pair<int, int> edge_pair;
								edge_pair.first = first_k;
								edge_pair.second = first_j;
								if (!set_of_edges.count(edge_pair)){
									set_of_edges.insert(edge_pair);
								}								
							}
						}								
					}
				}
	}
	/* 
		cout << "Edges" << set_of_edges.size() << endl;
	*/
	
}

void getExactMatches() {
	vector<int> tempVec;

	tempVec.push_back(0);

	for (int i = 1; i < combinedData.size(); ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(i);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(i);
		}
	}
	exactmatches.push_back(tempVec);
	totalUniqueRecords = exactmatches.size();
	cout << endl << "Total Exact Clusters: " << totalUniqueRecords << endl;
}

void getUniqueEntries() {

	ofstream unique_log_file;
	unique_log_file.open("unique_log");

	uniqueRecords.resize(totalUniqueRecords);

	for (size_t i = 0; i < totalUniqueRecords; i++)
    {
        uniqueRecords[i] = combinedData[exactmatches[i][0]];
		unique_log_file << uniqueRecords[i].first << " , " <<uniqueRecords[i].second;
	}
	unique_log_file.close();

}


int findRoot(int pointID, vector<int> &parentArr){
	if (parentArr[pointID] != pointID)
		parentArr[pointID]	= findRoot(parentArr[pointID], parentArr);

	return parentArr[pointID];
}

void findConnComp(int totalUniqueRecords)
{
	int i, rootU, rootV, edgeTotal;
    vector<int> parentArr;

	for(i = 0; i < totalUniqueRecords; ++i)
	{
		parentArr.push_back(i);
	}

	edgeTotal	= edgeList.size();

    BOOST_FOREACH(pair p, edgeList) {
        rootU	= findRoot(p.first, parentArr);
		rootV	= findRoot(p.second, parentArr);

		if(rootU != rootV)
		{
            parentArr[rootV] = rootU;
		}
    }
    int componentsInClusters = 0;
    for(int i = 0; i < parentArr.size(); i++) {
        int root;
        if (parentArr[i] == i) {
            // cout<< "i: " << i << " Parent: "<< parentArr[i] << " Val: " << uniqueRecords[i].first << " String "<< uniqueRecords[i].second<<endl;
            root = i;
        } else {
            root = findRoot(i, parentArr);
        }
        if (!approxConnectedComponents.count(root)) {
            vector<int> compononents;
            compononents.push_back(root);
            approxConnectedComponents[root] = compononents;
        } else {
            approxConnectedComponents[root].push_back(i);
            componentsInClusters++;
        }
    }
    cout<< "#Connected Components: " << approxConnectedComponents.size() <<endl;
    cout<< "#Total Non Root Nodes in graph: " << componentsInClusters << endl;
}

bool isLinkageOk(vector<string> &a, vector<string> &b, int threshold)
{
    //int dist = 0;
	int name_dist = calculateBasicED(a[1], b[1], 1);
    //dist +=name_dist;
    if (name_dist <= threshold) {
        int dod_dist = calculateBasicED(a[2], b[2], 1);
        //dist+=dod_dist;   
        if (dod_dist <= threshold) {
            int dob_dist = calculateBasicED(a[3], b[3], 1);
            //dist+=dod_dist;
            // if(dist==0) {
            //     //Self edge?
            //     return false;
            // }
            if (dob_dist <= threshold) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}

void findFinalConnectedComp(int intraCompDist, std::vector<string> x_prime_t) {
    
	int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
		totalNodes+=numComponents;
        bool distmat[numComponents][numComponents];
        std::vector<std::vector<string>> dataArr(numComponents); 
		// to make cache-efficient, keep records in a row


        for(int c = 0; c < p.second.size(); c++) {			
            dataArr[c] = vec2D[uniqueRecords[p.second[c]].first];
		};

		for (int i = 0; i < numComponents; i++) {
            distmat[i][i] = true;
            //cout << "Till - Here" << i << endl;
			for (int j = i+1 ; j < numComponents; j++)
            {
                //cout << "Till - Here" << dataArr[i][1] << dataArr[j][1] << endl;
				if (isLinkageOk(dataArr[i], dataArr[j], intraCompDist)) {
                    //cout << dataArr[i][1] << "," << dataArr[j][1] << endl;
					distmat[i][j] = true;
                    distmat[j][i] = true;
                } else {
                    distmat[i][j] = false;
                    distmat[j][i] = false;
                }
            }
        }
        bool nodesConsidered[numComponents];
        for(int i = 0; i < numComponents; i++) {
            nodesConsidered[i] = false;
        }
		

        for(int i=0; i < numComponents; i++) {
            if(nodesConsidered[i] == false) {
                vector<int> connectedComponent;
                connectedComponent.push_back(p.second[i]);
                nodesConsidered[i] = true;
                for(int j=0; j < numComponents; j++) {
                    if ((distmat[i][j] == true) && (nodesConsidered[j]==false)) {
						connectedComponent.push_back(p.second[j]);
                        nodesConsidered[j] = true;
                    }
                }
                finalConnectedComponents.push_back(connectedComponent);
            }
        }
    }
    cout<< "Total Nodes: "<< totalNodes << endl;
    cout<< "Total Components: "<< finalConnectedComponents.size() << endl;
}

void printApproximateCluster( std::vector<string> &x_prime_n) {
    int count = 0; 
    for (auto const& p : approxConnectedComponents) {
		//cout << "Approximate Cluster Size" << p.second.size() << endl;
        if (p.second.size() > 1) {
			for (int i=0; i < p.second.size(); i++) {
            	cout << uniqueRecords[p.second[i]].second << ",";
        	}
        	cout<< endl;
        	count++;
		}
		
        if (count > 50) break;
    }
}

void radixSort(vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string>> tempArr(numRecords);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(256, 0);
		
		for (int j = 0; j < numRecords; ++j) {
			countArr[(strDataArr[j].second)[i]]++;
		}
		
		for (int k = 1; k < 256; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = numRecords - 1; j >= 0; --j)
			tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
		
		for (int j = 0; j < numRecords; ++j)
			strDataArr[j]	= tempArr[j];
	}
}


void printFinalConnectedClusters(std::vector<std::vector<string>> &vec2D) {
	//cout << vec2D.size() << endl;
    for(int i=0; i < 100; i++) {
        for(int j=0; j < finalConnectedComponents[i].size(); j++) {
			 for(int k = 0; k < exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                //cout << "CC:" << finalConnectedComponents[i][j] << endl;
				cout << vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][1] << ",";
             }
    	}
		cout << endl;
	}
}

void writeFinalConnectedComponentToFile(string& result_file_name,std::vector<std::vector<string>> &vec2D) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i < finalConnectedComponents.size(); i++) {
		for(int j = 0; j < finalConnectedComponents[i].size(); j++) {			
			for(int k = 0; k < exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                out_file << vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][0] << ",";
            }
		}
        out_file << "\n";
	}
	out_file.close();
}

void printSortedRecords() {
    for (int i =0; i < 50; i++) {
    	cout<< combinedData[i].first << " " << combinedData[i].second << endl;
    }
}

void printEdges( std::vector<string> &x_prime_t ) {
    ofstream log_file;
	log_file.open("log_file");
    BOOST_FOREACH(pair p, edgeList) {
        string U	= x_prime_t[p.first];
		string V	= x_prime_t[p.second];
        log_file<< "U: "<< p.first << " V: "<< p.second <<endl;
        log_file<< "U: "<< U <<endl;
        log_file<< "V: "<< V << endl; 
        log_file<< endl;
    }
    log_file.close();
}


void buildGraph (vector<vector<string> > vec2D,Graph &graph,long long int recordTotal){

	
	vector<Edge> edges;

	for( auto i = 0; i < (recordTotal-1)/10; i++) {
		for (auto j = i+1; j < (recordTotal/10); j++) {
			int edit_distance = calculateBasicED(vec2D[i][1], vec2D[j][1], 1);
			if (edit_distance <= 2) {
				Edge edge(i, j);
				edges.push_back(edge);
			}
		}
	}

	graph.insertEdges(edges);
	
}



int main(void)
{

  
  vector<std::string> fileList = getInputFileList(directoryPath);
  
  struct timespec start, finish;
  double elapsed;
  
  /*
	cout << fileList.size() << endl;
  	for (int i = 0; i < fileList.size();i++){
		cout << fileList[i] << endl;
  	}
  */
  

  vec2D_original  = getData(fileList);
  vec2D = vec2D_original;

  
  double startTime_sort = omp_get_wtime();

  sort_data_parallel(vec2D,1,2);
  
  double endTime_sort = omp_get_wtime();
  
  double time_sorting = endTime_sort - startTime_sort;

  std::vector<string> x_prime;
 
  int totalsize = 0;
  totalUniqueRecords = x_prime.size();
  totalRecords = vec2D.size();
  
 
  getCombinedData();
  //getExactMatches();
  
  radixSort(combinedData);
  
  //printSortedRecords();
	

//   clock_t start_serial = clock();

	//clock_gettime(CLOCK_MONOTONIC, &start);


	struct timespec start_serial, finish_serial;
 	double elapsed_serial;

	std::vector<string> x_prime_serial;

	// clock_gettime(CLOCK_MONOTONIC, &start_serial);

	double startTime_serial = omp_get_wtime();

	//find_group_representative_serial(vec2D,0,vec2D.size(),x_prime_serial);

	double stopTime_serial = omp_get_wtime();

	// clock_gettime(CLOCK_MONOTONIC, &finish_serial);
	// elapsed_serial = (finish_serial.tv_sec - start_serial.tv_sec);

 	cout << "Elapsed-serial" << (stopTime_serial - startTime_serial) << endl; 
 	cout << "Size" << x_prime_serial.size() << endl;

	std::vector<string> x_prime_n;
	
	std::vector<string> x_prime_t;

	double startTime_p = omp_get_wtime();
	
	generate_x_prime_reprentative(vec2D,x_prime_t,x_prime_index);

	double endTime_p = omp_get_wtime();
    
	double time_generation_x_prime = endTime_p - startTime_p;

	totalUniqueRecords = exactmatches.size();

	getUniqueEntries();
  
  	cout << " The Unique Records" << uniqueRecords.size() << endl;
	
	ofstream log_file;
  	log_file.open("log_file_sorted_data");
  
  	for(int i=0; i < x_prime_t.size(); i++){
  		log_file <<  x_prime_t[i]  << endl;
  	}
  
  	log_file.close();

	cout << "Elapsed - P" << (endTime_p - startTime_p) << endl;
 	cout << "Size" << x_prime_t.size() << endl;

	
	

	double start_time_block = omp_get_wtime();	
	
	perform_blocking(x_prime_t,3,0);

	int blockArraySize  = 0;
	int blockArrsize = blockArray.size();
	for (int i = 0; i < blockArrsize; i++){
		blockArraySize += blockArray[i].size();
	}

	double end_time_block = omp_get_wtime();

	
	
	cout << "BlockArray Size" << blockArraySize << endl;

		

	cout << "Elapsed - Blocking" << (end_time_block - start_time_block) << endl;

	double time_blocking = end_time_block - start_time_block;
	

	double start_time_edge_generation = omp_get_wtime();
	
	generate_edges_concurrent(vec2D,x_prime_t,x_prime_index,blockArraySize);
	
	double end_time_edge_generation = omp_get_wtime();

	//std::sort(edgeList.begin(),edgeList.end());
	//edgeList.erase(unique( edgeList.begin(), edgeList.end() ), edgeList.end());
	
	cout << "Elapsed - EdgeList" << (end_time_edge_generation - start_time_edge_generation) << endl;
	cout << edgeList.size() << endl;
	
	double time_edge_generation = end_time_edge_generation - start_time_edge_generation;

	printEdges(x_prime_t);

	// double start_time_edge_serial = omp_get_wtime();
	// generate_edges_serial(vec2D,x_prime_t,x_prime_index,blockArraySize);
	// double end_time_edge_serial = omp_get_wtime();
	// cout << "Elapsed - EdgeList - Serial" << (end_time_edge_serial - start_time_edge_serial) << endl;

	double start_time_connected_serial = omp_get_wtime();

	findConnComp(x_prime_t.size());

	double end_time_connected_serial = omp_get_wtime();
	cout << "Elapsed - Connected Components - Serial" << (end_time_connected_serial - start_time_connected_serial) << endl;
	
	double time_connected_comp = end_time_connected_serial - start_time_connected_serial;


	//printApproximateCluster(x_prime_t);


	double start_time_final_connected_serial = omp_get_wtime();
	
	findFinalConnectedComp(1,x_prime_t);
	
	double end_time_final_connected_serial = omp_get_wtime();
	
	cout << "Elapsed - Final Connected Components - Serial" << (end_time_final_connected_serial - start_time_final_connected_serial) << endl;
	
	double time_final_connected = end_time_final_connected_serial - start_time_final_connected_serial;

	//printFinalConnectedClusters(vec2D);

	
	string out_name1 = out_file_path + fileName;
	writeFinalConnectedComponentToFile(out_name1,vec2D);
	
	double time_linking = time_sorting + time_generation_x_prime + time_blocking + time_edge_generation + time_connected_comp + time_final_connected;
	cout << "Linking-Time Parallel" << time_linking;
	//printFinalConnectedClusters(vec2D);
	return 0;
}



/*

RESULTS

Elapsed-serial  0.4523
Elapsed - P     0.151133

*/
