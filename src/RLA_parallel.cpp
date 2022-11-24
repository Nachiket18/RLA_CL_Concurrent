#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <fstream>
#include <dirent.h>
#include <boost/algorithm/string.hpp>
#include <execution>
#include <bits/stdc++.h>
#include <omp.h>

#include <algorithm>

#include <mutex>
#include <shared_mutex>
#include <tbb/tbb.h>

#include "RLA_parallel.h"


int numThreads = 6;
int lmerUsed = 3;


std::shared_mutex read_mutex;
std::mutex mut;


tbb::concurrent_vector<int> deDuplicateIndexData;
std::vector<vector<int>> blockArray;



// helps edit distance calculation in calculateBasicED()

int calculateBasicED2(string& str1, string& str2, int threshRem)
{
	int row, col, i, j;
	vector<vector<int> > matArr;

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




vector<vector<string>,tbb::cache_aligned_allocator<string> > getData(vector<string> files) {
    string line;
    vector<vector<string>,tbb::cache_aligned_allocator<string>> vec2D;

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

void sort_data_parallel(vector<vector<string> , tbb::cache_aligned_allocator<string> > &data_vector,int column){

	std::sort(data_vector.begin(), data_vector.end(),
              [column](std::vector<string> const& v1, std::vector<string> const& v2)
              {
                  return v1[column] < v2[column];
              });
	// tbb::parallel_sort(data_vector.begin(),data_vector.end(),[column](std::vector<string> const& v1, std::vector<string> const& v2)
	//          {
    //             return v1[column] < v2[column];
    //         });
	

}




/*
	1. Find the unique elements in the Data Vector Array for specified range (start,limit)
	2. Used by multiple threads (on seperate start,limit values)
	3. Since reading can be done concurrently, we have used shared mutex,
*/


vector<string> find_group_representative(vector<vector<string>> &data_vector,int start,int limit){
	
	vector<string> temp;
	temp.reserve((limit-start)/3);
	vector<int> indexTemp; // use for capturing the index of the deduplicate data. Used further while blocking 

	//cout << "Thread - id" << std::this_thread::get_id() << endl;
	//read_mutex.lock_shared();
	
	//std::shared_lock<std::shared_mutex> lock(read_mutex);
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
	
	return temp;
	//read_mutex.unlock_shared();
	
	//cout << "Temp Size" << temp.size() <<  endl;
	//cout << "Thread - id- New" << std::this_thread::get_id() << endl;
	
	//std::copy(temp.begin(), temp.end(),x_prime.grow_by(temp.size()));

	//add_representative_x_prime(temp,indexTemp,x_prime);

} 


// Not implemented
void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int>>& blockArr)
{
    
	
}

// Not Implemented

/*

	void perform_blocking(vector<vector<string> > &data_vector){

	if (lmerUsed < 3)
		lmerUsed = 3;
	
	int blockTotal 	= pow(26, lmerUsed);
	string strSample;
	blockArray.resize(blockTotal);
	
	for ( int i = 0; i < deDuplicateIndexData.size(); i++){

		vector<int> tmp = deDuplicateIndexData[i];
		for (int j = 0; j < tmp.size(); j++){
			
			for (int k = tmp[j]; k <= tmp[j+1]; k++){
			
			}
		
		}

	}
}


*/



/*
	Performing the de-duplication of records pre blocking serially.
	Was used to check correctness of the parallel implementation and difference in serial and parallel run time. 
*/

void find_group_representative_serial(vector<vector<string>,tbb::cache_aligned_allocator<string> > &data_vector,int start,int limit, vector<string> &x_prime){
	
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

void find_group_representative_openmp(vector<vector<string>,tbb::cache_aligned_allocator<string> > &data_vector,int start,int limit, vector<string> &x_prime){
	
	//cout << "Thread - id" << std::this_thread::get_id() << endl;

		
	#pragma omp parallel num_threads(4)
	{

		vector<string> temp;
		temp.push_back(data_vector[start][1]);
	
		#pragma omp for 
		for(int i = start+1; i < limit; i++) {	

			if (data_vector[i][1] != data_vector[i-1][1]){
				temp.push_back(data_vector[i][1]);
			}
			
		}
		//cout << "Till Here" << endl;
		#pragma omp critical
		
		x_prime.insert(x_prime.end(), temp.begin(), temp.end());

	}

	int totalS = 0;	
	cout << "total size" << x_prime.size() << endl;

}




/*
	Adds the data to x_prime - A vector that stores the deduplicated data pre- blocking
	which corresponds to Step 2 of the algorithm.
	Uses mutex as multiple threads are going to write to x_prime which is Critical Section.
*/

void add_representative_x_prime(vector<string> &temp,vector<int> indexTemp,tbb::concurrent_vector<string> &x_prime){

	//read_mutex.lock();
	//std::lock_guard<mutex> lock(mut);

	//x_prime.grow_by(temp.size());
	std::copy(temp.begin(), temp.end(),x_prime.grow_by(temp.size()));

	//temp.shrink_to_fit();
	//x_prime.push_back(temp);

	//deDuplicateIndexData.push_back(indexTemp);

	//read_mutex.unlock();

}


void generate_x_prime_reprentative(vector<vector<string> , tbb::cache_aligned_allocator<string> > &data_vector,vector<string> &x_prime){

	std::vector<std::future<vector<string>>> threads;

	int numElements = data_vector.size() / numThreads;
	//cout << "numElements" << numElements << endl;
	int remainingElements = data_vector.size();
	int currentRecordPointer = 0;
	


	for (int i=0; i < numThreads; i++ ){
		//cout << "i:" << i << endl;
		//cout << "Record Pointer" << currentRecordPointer << endl;
		//cout << "Remaining Elements" << remainingElements << endl;
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
			
			//cout << "Till here" << endl;
			//cout << "Remaining elements" << remainingElements << endl;
		}
		
		auto start_t = data_vector.begin() + start;
        auto end_t = data_vector.begin() + end + 1;

		vector<vector<string>> result_tmp(end - start + 1);
		copy(start_t, end_t, result_tmp.begin());
		
		threads.push_back(std::async(std::launch::async,&find_group_representative,std::ref(result_tmp),0,result_tmp.size()));
	}
	
	// for (auto i = 0; i < threads.size(); i++){
	// 	if (threads[i].joinable()) {
	// 		threads[i].join();
	// 	}	
    // }

	for (auto i = 0; i < threads.size(); i++){
		//x_prime.reserve(x_prime.size() + threads[i].get().size());
		vector<string> tmp = threads[i].get();
		x_prime.insert(x_prime.end(),tmp.begin(),tmp.end());
    }
	

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
  
  cout << fileList.size() << endl;
  

  vector<vector<string>,tbb::cache_aligned_allocator<string> >  vec2D  = getData(fileList);
  cout << "Size:" << vec2D.size() << endl;
  
  sort_data_parallel(vec2D,1);

  std::vector<string> x_prime;
 
  


  
  int totalsize = 0;

//   cout << "Total-Size-P" << x_prime.size() << endl;
  
  // 640498,517764
	

//   clock_t start_serial = clock();
    	
//   find_group_representative_serial(vec2D,0,vec2D.size(),x_prime); //
		
//   clock_t duration_serial = clock() - start_serial;	

//   cout << "duration_serial" << duration_serial << endl;
	
	
	// clock_t start_uniq = clock();

	// for(int i = 0; i < vec2D.size(); i++) {	

	// 	x_prime_n.push_back(vec2D[i][1]);
	// }

	// de_duplication_parallel_uniq(x_prime_n);

	// clock_t duration_uniq = clock() - start_uniq;


	// cout << " X prime Size " << x_prime_n.size() << endl; 
	// cout << "Duration Uniq" << duration_uniq << endl;

	// Duration_ Parallel --> 5724118 
	                          
	// Duration_ Serial:      723952


		// Duration_ Parallel5916083

        // Performance counter stats for './RLA_parallel':

        // 348,001,239,875      L1-dcache-loads                                             
        // 986,871,106      L1-dcache-load-misses     #    0.28% of all L1-dcache accesses
        // 1,019      context-switches                                            

        // 140.917704162 seconds time elapsed

        // 145.433407000 seconds user
        // 0.219965000 seconds sys

	clock_t start_serial = clock();
	
	find_group_representative_serial(vec2D,0,vec2D.size(),x_prime);
	
	clock_t duration_serial = clock() - start_serial;

	cout << "\nDurationSerial:" << duration_serial << endl;


	
	std::vector<string> x_prime_n;

	clock_t start_mp = clock();
	
	x_prime_n.push_back("#");

	find_group_representative_openmp(vec2D,0,vec2D.size(),x_prime_n);
	
	clock_t duration_mp = clock() - start_mp;

	cout << "\n Duration OpenMP:" << duration_mp << endl;


  




	
	
	// cout << "totalS" << totalsize;

	// clock_t duration_serial = clock() - start_serial;
	// cout << "Duration, DurationSerial:" <<duration  << "," << duration_serial << endl;


//   for (int i =0;i<4;i++){
// 	for (int j=0;j<2;j++){
// 		cout << vec2D[i][j] << "\t";
// 	}
	
//   }

  //Graph graph(recordTotal);
  //buildGraph(vec2D,graph,recordTotal);
  
   return 0;
}






