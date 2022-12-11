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
static const int ALPHABET_SIZE_LIST[] = {26, 10, 36, 256};

std::shared_mutex read_mutex;
std::mutex mut;


std::vector<std::vector<int> > blockArray;
std::set<pair<int, int> > edgeList;
map<int, vector<int> > approxConnectedComponents;


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

		//if(procID == 1 && checkTemp == 3164)
			//	cout << str1 << " -- " << str2 << " rt " << dist << endl;

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

		//if(procID == 1 && checkTemp == 3164)
			//cout << str1 << " -- " << str2 << " hukjhk " << matArr[row - 1][diagonal] << endl;

		return matArr[row - 1][diagonal];

	}
}




vector<vector<string>> getData(vector<string> files) {
    string line;
    vector<vector<string>> vec2D;

    // Read from the text file
    
    for (int i=0; i < files.size();i++) {
        
        if (files[i] != "") {
            ifstream records(files[i]);
            //cout << "Here:" <<endl;
            while(getline(records,line)){

                vector<string> result;
                boost::split(result, line, boost::is_any_of("\t"));
                vector<string> vec;
                string name = result[1] + result[2];
                boost::to_lower(name);
                //cout << "Name:" << result[0] <<endl;            
                vec.push_back(result[0]);
                vec.push_back(result[1]);
                vec.push_back(result[2]);
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

void sort_data_parallel(vector<vector<string> > &data_vector,int column){

	std::sort(std::execution::par,data_vector.begin(), data_vector.end(),
              [column](std::vector<string> const& v1, std::vector<string> const& v2)
              {
                  return v1[column] < v2[column];
              });
	// tbb::parallel_sort(data_vector.begin(),data_vector.end(),[column](std::vector<string> const& v1, std::vector<string> const& v2)
	//          {
    //             return v1[column] < v2[column];
    //         });
	
	// 12619012
	// 6568083

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
	
	//return temp;
	//read_mutex.unlock_shared();
	
	//cout << "Temp Size" << temp.size() <<  endl;
	//cout << "Thread - id- New" << std::this_thread::get_id() << endl;
	
	//std::copy(temp.begin(), temp.end(),x_prime.grow_by(temp.size()));	

} 


// Not implemented
void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int>>& blockArr)
{
    
	
}

// Not Implemented



void perform_blocking(std::vector<string>  &x_prime_t,int lmerUsed,int type){

	if (lmerUsed < 3)
		lmerUsed = 3;
	
	int blockTotal 	= pow(26, lmerUsed);
	string strSample;
	blockArray.resize(blockTotal);
	int numElements = x_prime_t.size() / numThreads;

	#pragma omp parallel num_threads(6)
	{	
		std::vector< std::vector<int> > temp_blocks;
		temp_blocks.resize(blockTotal);

		#pragma omp for 
		for (int i = 0; i < x_prime_t.size(); i++){
		
		//cout << "x_prime" << x_prime_t[i] << endl;
			for (int j = 0; j < (x_prime_t[i].size() - lmerUsed + 1); j++){
				int blockID=0;
				for (int k = 0; k < lmerUsed; ++k)
					if (type == 0)
						blockID	+= (x_prime_t[i][j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type],lmerUsed - k - 1);
					else if(type == 1)
						blockID	+= (x_prime_t[i][j + k] - 48) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
					else if(type == 2)
					{
						if(x_prime_t[i][j + k] >= 97)
							blockID	+= (x_prime_t[i][j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
						else
							blockID	+= (x_prime_t[i][j + k] - 22) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1); // 48 - 26
					}
				
				if( !( blockID < 0 || blockID >= blockTotal ) ) {
						temp_blocks[blockID].push_back(i);	 
				}
					
			}
			

		}
		//std::cout << "Till Here" << endl;
		#pragma omp critical
		{		
			int thread_data = 0;
			for (int j=0; j < temp_blocks.size(); j++){
				if (temp_blocks[j].size() != 0){
					thread_data += temp_blocks[j].size();
					for (int k = 0; k < temp_blocks[j].size();k++){
						blockArray[j].push_back(temp_blocks[j][k]);
					}
				}
			}
			//blockArray.insert(blockArray.end(),temp_blocks.begin(),temp_blocks.end());		
		}
	}
}






/*
	Performing the de-duplication of records pre blocking serially.
	Was used to check correctness of the parallel implementation and difference in serial and parallel run time. 
*/

void find_group_representative_serial(vector<vector<string> > &data_vector,int start,int limit, vector<string> &x_prime){
	
	//cout << "Thread - id" << std::this_thread::get_id() << endl;
	
	x_prime.push_back(data_vector[start][1]);
	//std::shared_lock<std::shared_mutex> lock(read_mutex);
	
	for(int i = start+1; i < limit; i++) {	

		if (data_vector[i][1] == data_vector[i-1][1]){
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
	
	//cout << "Thread - id" << std::this_thread::get_id() << endl;

		
	#pragma omp parallel num_threads(4)
	{

		vector<string> temp;
		vector<int> indexTemp;
		
		temp.push_back(data_vector[start][1]);
		indexTemp.push_back(start);
	
		#pragma omp for 
		for(int i = start+1; i < limit; i++) {	

			if (data_vector[i][1] != data_vector[i-1][1]){
				temp.push_back(data_vector[i][1]);
				indexTemp.push_back(i);
			}
			
		}
		//cout << "Till Here" << endl;
		#pragma omp critical
		{
			x_prime.insert(x_prime.end(), temp.begin(), temp.end());
			x_prime_index.insert(x_prime_index.end(), indexTemp.begin(), indexTemp.end());
		}
		
	}

	int totalS = 0;	
	cout << "total size" << x_prime.size() << endl;

}




/*
	Adds the data to x_prime - A vector that stores the deduplicated data pre- blocking
	which corresponds to Step 2 of the algorithm.
	Uses mutex as multiple threads are going to write to x_prime which is Critical Section.
*/

void add_representative_x_prime(vector<string> &temp,vector<int> indexTemp,std::vector<string> &x_prime){

	//read_mutex.lock();
	std::lock_guard<mutex> lock(mut);


	x_prime.insert(x_prime.end(), temp.begin(), temp.end());

}





void generate_x_prime_reprentative(vector< vector<string> > &data_vector,vector<string> &x_prime, vector<int> &x_prime_index){

	std::vector<std::thread> threads;
	
	int numElements = data_vector.size() / numThreads;
	//cout << "numElements" << numElements << endl;
	int remainingElements = data_vector.size();
	int currentRecordPointer = 0;
	std::vector<int> threadDivideData;
	threadDivideData.push_back(0);
	x_prime_index.resize(data_vector.size()/4);

	for (int i=0; i < numThreads; i++){

		int start = 0;
		int end = 0;
		std::vector<std::string> out;
		if ( (int)data_vector.size() < (currentRecordPointer + numElements) ){
			
			start = currentRecordPointer;
			end = data_vector.size();
			currentRecordPointer = end;
			remainingElements = 0;
			
		}
		else {
			//cout << "Here - 2" << endl; 
			int idealRange = (currentRecordPointer) + numElements;
			int tmp = (currentRecordPointer) + (numElements+1);
			
			while(data_vector[tmp][1] == data_vector[idealRange][1]){
				tmp +=1;
			}
			//cout << "tmp" << tmp << endl;

			start = currentRecordPointer;
			end = tmp;
			
			currentRecordPointer = tmp;
			remainingElements = data_vector.size() - currentRecordPointer;			
			
		}
		threadDivideData.push_back(end);
	
	
	//threads.push_back(std::thread(&find_group_representative,std::ref(data_vector),start,end,std::ref(x_prime)));
	}

	// int start = 0;
	
	#pragma omp parallel num_threads(6)
	{
		vector<string> temp;
		vector<int> indexTemp;	

		#pragma omp for 
		for (int i=1; i < threadDivideData.size(); i++){
			
			//find_group_representative(data_vector,start,threadDivideData[i],x_prime,temp,indexTemp)
			
			temp.push_back(data_vector[threadDivideData[i-1]][1]);
			indexTemp.push_back(threadDivideData[i-1]);
	
			for(int j = threadDivideData[i-1]+1; j < threadDivideData[i]; j++) {	

				if (data_vector[j][1] != data_vector[j-1][1]){
					temp.push_back(data_vector[j][1]);
					indexTemp.push_back(j);
				}
		
			}
			
		
			//cout << "Here" << temp.size() << endl;
			#pragma omp critical
			{
				int size_x_prime = x_prime.size();
				x_prime.insert(x_prime.end(), temp.begin(), temp.end());
				x_prime_index.insert(x_prime_index.end(),indexTemp.begin(),indexTemp.end());

			}
		
		}

	}
	
}



void generate_edges_concurrent(std::vector<std::vector<string>> &data_vector , std::vector<string>  &x_prime_t,std::vector<int> &x_prime_index,int dataSize){

	int dataPerThread = dataSize / numThreads;
	cout << "Data=perThread" << dataPerThread << endl;
	std::vector<int> threadDivideDataBlocks;
	threadDivideDataBlocks.resize(numThreads+1);

	threadDivideDataBlocks.push_back(0);

	int currentThreadData = 0;

	for(int i = 0; i < blockArray.size(); i++){

		if ( ( currentThreadData + blockArray[i].size() ) < dataPerThread){
			currentThreadData += blockArray[i].size();
		}
		else{
			threadDivideDataBlocks.push_back(i);
			currentThreadData = 0;
		}

	}


	#pragma omp parallel num_threads(6) 
	{
		std::set<pair<int, int> > set_of_edges;

		#pragma omp for
		for(int i = 1; i < threadDivideDataBlocks.size(); i++){

			for( int m = threadDivideDataBlocks[i-1]; m < threadDivideDataBlocks[i]; m++) {
				
				std::vector<int> block = blockArray[m];
				std::sort(block.begin(),block.end());
				block.erase(unique( block.begin(), block.end() ), block.end());

				for (int j = 0; j < block.size(); j++){
					for (int k = j+1; k < block.size();k++){
						int first_j = block[j];
						int first_k = block[k];
						int dist_last_name = calculateBasicED(x_prime_t[first_j],x_prime_t[first_k],1);
						
						if (dist_last_name > 1)
							continue;
						
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

		}

		#pragma omp critical
		{
				set<std::pair<int,int>>::iterator itr;
				edgeList.insert(set_of_edges.begin(), set_of_edges.end());
				cout << set_of_edges.size() << endl;

		}

	}
	
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
	cout << "Number of edges " << edgeTotal  << endl;

    BOOST_FOREACH(pair p, edgeList) {
        rootU	= findRoot(p.first, parentArr);
		rootV	= findRoot(p.second, parentArr);

		if(rootU != rootV)
		{
            parentArr[rootV] = rootU;
		}
    }
    int componentsInClusters = 0;
    for(int i =0; i<parentArr.size(); i++) {
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
    cout<< "#Connected Components: " << approxConnectedComponents.size()<<endl;
    cout<< "#Total Non Root Nodes in graph: " << componentsInClusters << endl;
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


// BENCHMARK(BM_find_group_representative);

// BENCHMARK_MAIN();








int main(void)
{

  const char* directoryPath = "/home/nachiket/RLA_CL_EXTRACT/data/in_1000k";
  vector<std::string> fileList = getInputFileList(directoryPath);
  
  struct timespec start, finish;
  double elapsed;

  cout << fileList.size() << endl;
  

  std::vector<std::vector<string>>  vec2D  = getData(fileList);
  cout << "Size:" << vec2D.size() << endl;
  
  

  sort_data_parallel(vec2D,1);



  std::vector<string> x_prime;
 
  int totalsize = 0;

//   cout << "Total-Size-P" << x_prime.size() << endl;
  
  // 640498,517764
	

//   clock_t start_serial = clock();

	//clock_gettime(CLOCK_MONOTONIC, &start);


	struct timespec start_serial, finish_serial;
 	double elapsed_serial;

	std::vector<string> x_prime_serial;

	// clock_gettime(CLOCK_MONOTONIC, &start_serial);

	double startTime_serial = omp_get_wtime();

	find_group_representative_serial(vec2D,0,vec2D.size(),x_prime_serial);

	double stopTime_serial = omp_get_wtime();

	// clock_gettime(CLOCK_MONOTONIC, &finish_serial);
	// elapsed_serial = (finish_serial.tv_sec - start_serial.tv_sec);

 	cout << "Elapsed-serial" << (stopTime_serial - startTime_serial) << endl;
 
 	cout << "Size" << x_prime_serial.size() << endl;

	std::vector<string> x_prime_n;
	std::vector<int> x_prime_index;
	// clock_t start_mp = clock();
	
	// x_prime_n.push_back("#");

	// double startTime_mp = omp_get_wtime();

	// find_group_representative_openmp(vec2D,0,vec2D.size(),x_prime_n,x_prime_index);
	
	// double endTime_mp = omp_get_wtime();


	// cout << "Elapsed-MP" << (endTime_mp - startTime_mp) << endl;
 	// cout << "Size" << x_prime_n.size() << endl;


	std::vector<string> x_prime_t;

	double startTime_p = omp_get_wtime();
	
	generate_x_prime_reprentative(vec2D,x_prime_t,x_prime_index);

	double endTime_p = omp_get_wtime();

	cout << "Elapsed - P" << (endTime_p - startTime_p) << endl;
 	cout << "Size" << x_prime_t.size() << endl;


	double start_time_block = omp_get_wtime();	
	
	perform_blocking(x_prime_t,3,0);

	int blockArraySize  = 0 ;
	for (int i = 0; i < blockArray.size(); i++){
		blockArraySize += blockArray[i].size();
	}
	
	cout << "BlockArray Size" << blockArraySize << endl;

	double end_time_block = omp_get_wtime();	

	cout << "Elapsed - Blocking" << (end_time_block - start_time_block) << endl;

	double start_time_edge_generation = omp_get_wtime();

	generate_edges_concurrent(vec2D,x_prime_t,x_prime_index,blockArraySize);

	double end_time_edge_generation = omp_get_wtime();

	cout << "Elapsed - EdgeList" << (end_time_edge_generation - start_time_edge_generation) << endl;
	cout << edgeList.size() << endl;
    
	return 0;
}



/*

RESULTS

Elapsed-serial  0.4523
Elapsed - P     0.151133

*/
