//=======================================================================
//
// Inspired and based on 
// https://www.boost.org/doc/libs/1_43_0/libs/graph/example/astar-cities.cpp
//
//=======================================================================

#include <boost/graph/astar_search.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/graph/graphviz.hpp>
#include <time.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <unordered_map> 
#include <math.h>    // for sqrt

using namespace boost;
using namespace std;

#define EXCHANGE_FILE_FULL_NAME "C:\\SNU\\MSData\\exchange.txt" 

typedef std::pair<int, int>                                     Edge;
typedef double                                                  Cost;
typedef adjacency_list< listS, vecS, undirectedS, no_property,
                        property<edge_weight_t, Cost> >         Graph;
typedef property_map<Graph, edge_weight_t>::type                WeightMap;
typedef Graph::edge_descriptor                                  EdgeDescr;
typedef Graph::vertex_descriptor                                VertexDescr;

double EPS = 0.00001;

//-----------------------------------------------------------------------------
// auxiliary types
struct Vertex
{
    Vertex(double x, double y, double z = 0.) :x(x), y(y), z(z) {}
    Cost dist(const Vertex& other) const
    {
        Cost dx = x - other.x;
        Cost dy = y - other.y;
        Cost dz = z - other.z;
        return ::sqrt(dx * dx + dy * dy + dz * dz);

    }
    double x, y, z;
};
class VertexLessFn
{
public:
    bool operator() (const Vertex& v0, const Vertex& v1) const
    {
        if (abs(v0.x - v1.x) < EPS)
            if (abs(v0.y - v1.y) < EPS)
                if (abs(v0.z - v1.z) < EPS)
                    return false;
                else
                    return v0.z < v1.z;
            else
                return v0.y < v1.y;
        else
            return v0.x < v1.x;
    }
};
//-----------------------------------------------------------------------------
// euclidean distance heuristic
template <class Graph, class CostType, class LocMap>
class distance_heuristic : public astar_heuristic<Graph, CostType>
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    distance_heuristic(LocMap l, Vertex goal)
        : m_location(l), m_goal(goal) {}
    CostType operator()(Vertex u)
    {
        return m_location[m_goal].dist(m_location[u]);
    }
private:
    LocMap m_location;
    Vertex m_goal;
};


//-----------------------------------------------------------------------------
// exception for termination
struct found_goal {}; 

//-----------------------------------------------------------------------------
// visitor that terminates when we find the goal
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
    astar_goal_visitor(Vertex goal) : m_goal(goal) {}
    template <class Graph>
    void examine_vertex(Vertex u, Graph& g) {
        if (u == m_goal)
            throw found_goal();
    }
private:
    Vertex m_goal;
};

//-----------------------------------------------------------------------------
int register_vertex( const Vertex&                   v,
                     vector<Vertex>&                 vecVertices,
                     map<Vertex, int, VertexLessFn>& mapLookup )
{
    auto vi = mapLookup.find(v);
    int idx = -1;
    if( vi == mapLookup.end() )
    {
        vecVertices.push_back(v);
        idx = (int) vecVertices.size() - 1;
        mapLookup.insert({ v, idx });
    }
    else 
    {
        idx = vi->second;
    }
    return idx;
}

bool read_file( const string&   strFileName,
                vector<Vertex>& vecVertices,
                vector<Edge>&   vecEdges,
                vector<Cost>&   vecWeights )
{
    map<Vertex, int, VertexLessFn> mapLookup;
    ifstream data(strFileName);
    if (!data.is_open())
        return false;

    bool bFound = false;
    for( string strLine; getline(data, strLine); )
        if (strLine == "voronoi")
        {
            bFound = true;
            break;
        }
    if (!bFound)
        return false;

    int nOfEdges = -1;
    data >> nOfEdges;
    for (int i = 0; i < nOfEdges; ++i)
    {
        double px, py, qx, qy;
        int idMSumObj1, idxArc1, idMSumObj2, idxArc2;
        data >> px >> py >> qx >> qy >> idMSumObj1 >> idxArc1 >> idMSumObj2 >> idxArc2;
        Vertex vrtxP(px, py);
        Vertex vrtxQ(qx, qy);
        int idxP = register_vertex(vrtxP, vecVertices, mapLookup);
        int idxQ = register_vertex(vrtxQ, vecVertices, mapLookup);
        vecEdges.push_back({ idxP, idxQ });
        vecWeights.push_back(vrtxP.dist(vrtxQ));
    }
    return true;
}
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    vector<Vertex> vecVertices;
    vector<Edge> vecEdges;
    vector<Cost> vecWeights;
    if (!read_file(EXCHANGE_FILE_FULL_NAME, vecVertices, vecEdges, vecWeights))
        return -1;

    // create graph
    Graph theGraph(vecVertices.size());
    WeightMap weightmap = get(edge_weight, theGraph);
    for (std::size_t j = 0; j < vecEdges.size(); ++j)
    {
        EdgeDescr e;
        bool inserted;
        boost::tie(e, inserted) = add_edge(vecEdges[j].first, vecEdges[j].second, theGraph);
        weightmap[e] = vecWeights[j];
    }


    VertexDescr start = 0;
    VertexDescr goal = 18;


    cout << "Start vertex: " << start << endl;
    cout << "Goal vertex: " << goal << endl;

    vector<Graph::vertex_descriptor> p(num_vertices(theGraph));
    vector<Cost> d(num_vertices(theGraph));
    try {
        // call astar named parameter interface
        astar_search
        (theGraph, start,
            distance_heuristic<Graph, Cost, vector<Vertex>>
            (vecVertices, goal),
            predecessor_map(&p[0]).distance_map(&d[0]).
            visitor(astar_goal_visitor<VertexDescr>(goal)));
    }
    catch (found_goal)
    {
        // found a path to the goal
        list<VertexDescr> shortest_path;
        for (VertexDescr v = goal;; v = p[v]) {
            shortest_path.push_front(v);
            if (p[v] == v)
                break;
        }
        cout << "Shortest path from " << start << " to "
            << goal << ": ";
        list<VertexDescr>::iterator spi = shortest_path.begin();
        cout << start;
        for (++spi; spi != shortest_path.end(); ++spi)
            cout << " -> " << *spi;
        cout << endl << "Total travel distance: " << d[goal] << endl;
        return 0;
    }

    cout << "Didn't find a path from " << start << "to"
        << goal << "!" << endl;
    return 0;
}
//============================ END OF FILE ===================================