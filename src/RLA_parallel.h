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




int calculateBasicED2(string& str1, string& str2, int threshRem);
int calculateBasicED(string& str1, string& str2, int threshRem);
vector<vector<string>,tbb::cache_aligned_allocator<string> > getData(vector<string> files);
vector<std::string> getInputFileList(std::string directory);
void sort_data_parallel(vector<vector<string>,tbb::cache_aligned_allocator<string> > data_vector);
void buildGraph (vector<vector<string> > vec2D,Graph &graph,long long int recordTotal);
vector<string> find_group_representative(vector<vector<string>> &data_vector,int start,int limit);
void add_representative_x_prime(vector<string> &temp,vector<int>indexTemp,tbb::concurrent_vector<string> &x_prime);
void generate_x_prime_reprentative(vector<vector<string>,tbb::cache_aligned_allocator<string> > &data_vector,vector<string> &x_prime);
void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int> >& blockArr);