#include <string>
#include <vector>

using namespace std;

int threshold = 1;

struct Edge {
    int src, dest;
    Edge(int my_src, int my_dest) {
        src = my_src;
        dest = my_dest;
    }
};

class Graph
{
public:

    vector<int>  adjList;
    Graph(int n)
    {
        adjList.resize(n);
    }

    /*
        Default Constructor 
    */
    Graph(){
        
    }

    void insertEdges(vector<Edge> const &edges){

        for (auto &edge: edges)
        {
            adjList[edge.src] = edge.dest;
        }
    }


};


/*
 --------------- END USERS DO NOT MODIFY-------------------------------------------
*/

int calculateBasicED2(string& str1, string& str2, int threshRem);
int calculateBasicED(string& str1, string& str2, int threshRem);
vector<vector<string> > getData(vector<string> files);
vector<std::string> getInputFileList(std::string directory);
void sort_data_parallel(vector<vector<string>> data_vector);
void buildGraph (vector<vector<string> > vec2D,Graph &graph,long long int recordTotal);
//void find_group_representative(vector<vector<string>,tbb::cache_aligned_allocator<string>> &data_vector,int start,int limit,vector<string> &x_prime);
void find_group_representative(vector<vector<string>> &data_vector,int start,int limit,vector<string> &x_prime,vector<string> &temp,vector<int> &indexTemp);
void add_representative_x_prime(vector<string> &temp,vector<int>indexTemp,vector<string> &x_prime);
void generate_x_prime_reprentative(vector<vector<string> > &data_vector,vector<string> &x_prime);
void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int> >& blockArr);
bool isLinkageOk(vector<string> &a, vector<string> &b, int threshold);
int power_mod(int base,int pow,int p);
std::vector<int> generateIntKmersMap(string& str,int k);
double calculateBasicQgramIntMap(string& str1, string& str2,int threshRem);


/*
    ---------------------- USER INPUTS --------------------------------
*/ 

const char* directoryPath = "/home/nachiket/RLA_CL_EXTRACT/data/1m_ds/"; // The path to input data files
string out_file_path = "/home/nachiket/RLA_CL_parallel/"; // path to output files
string fileName = "output_parallel_rla";  // output file name
