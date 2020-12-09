#include <Python.h>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <iterator>
#include <utility>
#include <string>
#include <numeric>

//(nodeCount, edgeCount, edges)
static PyObject *buildNodeEdgeVector(PyObject *self, PyObject *args){
    int nodeCount, edgeCount;
    Py_buffer edgeBuf;
    if (!PyArg_ParseTuple(args, "iiy*", &nodeCount, &edgeCount, &edgeBuf)) {
        PyErr_SetString(PyExc_RuntimeError, "buildNodeEdgeVector incorrect arguments");
        return NULL;
    }

    //Initalize nodeEdgeVector
    auto nodeEdgeVector = new std::vector<std::unordered_map<int, std::list<int>>>();
    nodeEdgeVector->reserve(nodeCount);
    for(int i = 0; i < nodeCount; ++i){
        nodeEdgeVector->push_back(std::unordered_map<int, std::list<int>>());
    }

    //For each edge, add the connected node to the other node's list for both nodes
    int* edges = static_cast<int*>(edgeBuf.buf);
    for(int i = 0; i < edgeCount; ++i){
        (*nodeEdgeVector)[edges[i * 3]].try_emplace(edges[i * 3 + 2], std::list<int>());
        (*nodeEdgeVector)[edges[i * 3]][edges[i * 3 + 2]].push_back(edges[i * 3 + 1]);

        (*nodeEdgeVector)[edges[i * 3 + 1]].try_emplace(edges[i * 3 + 2], std::list<int>());
        (*nodeEdgeVector)[edges[i * 3 + 1]][edges[i * 3 + 2]].push_back(edges[i * 3]);
    }

    PyBuffer_Release(&edgeBuf);

    PyObject* nodeEdgeVectorCapsule = PyCapsule_New((void*)nodeEdgeVector, "nodeEdgeVectorPtr", NULL);
    PyCapsule_SetPointer(nodeEdgeVectorCapsule, (void*)nodeEdgeVector);

    return Py_BuildValue("O", nodeEdgeVectorCapsule);
}

static int nodeEdgeCount(std::unordered_map<int, std::list<int>> &nodeEdges){
    int count = 0;
    for(auto i = nodeEdges.begin(); i != nodeEdges.end(); ++i){
        count += static_cast<int>(i->second.size());
    }
    return count;
}

static double weight(int node1, int node2, std::vector<std::unordered_map<int, std::list<int>>> &nodeEdgeVector){
    auto node1UniqueKeys = std::set<int>();
    auto node2UniqueKeys = std::set<int>();
    std::transform(nodeEdgeVector[node1].begin(), nodeEdgeVector[node1].end(), std::inserter(node1UniqueKeys, node1UniqueKeys.begin()), [](auto pair){return pair.first;});
    std::transform(nodeEdgeVector[node2].begin(), nodeEdgeVector[node2].end(), std::inserter(node2UniqueKeys, node2UniqueKeys.begin()), [](auto pair){return pair.first;});

    auto commonSymptoms = std::vector<int>();
    std::set_intersection(node1UniqueKeys.begin(), node1UniqueKeys.end(), node2UniqueKeys.begin(), node2UniqueKeys.end(), std::inserter(commonSymptoms, commonSymptoms.end()));
        
    int equalEdges = 0;
    for(auto i = commonSymptoms.begin(); i != commonSymptoms.end(); ++i){
        for(auto edge1 = nodeEdgeVector[node1][*i].begin(); edge1 != nodeEdgeVector[node1][*i].end(); ++edge1){
            for(auto edge2 = nodeEdgeVector[node2][*i].begin(); edge2 != nodeEdgeVector[node2][*i].end(); ++edge2){
                if(node1 == *edge2 || node2 == *edge1){
                   equalEdges += 2;
                   break; 
                }else if(*edge1 == *edge2){
                    equalEdges += 1;
                   break; 
                }
            }
        }
    }

    int total = (nodeEdgeCount(nodeEdgeVector[node1]) + nodeEdgeCount(nodeEdgeVector[node2]));
    return (double)equalEdges / (total > 0 ? total : 1);
}

static void fillNodeSymptomCatagory(int symptom, std::unordered_map<int, std::set<int>> &symptomCatagoryMap, std::list<int> &catagories){
    for(auto i = catagories.begin(); i != catagories.end(); ++i){
        symptomCatagoryMap.try_emplace(symptom, std::set<int>());
        symptomCatagoryMap[symptom].insert(*i);
    }
}

struct WeightEdge{
    int n1;
    int n2;
    double weight;
    WeightEdge(int n1, int n2, double weight): n1(n1), n2(n2), weight(weight) {}
};

static PyObject* buildWeightGraph(PyObject *self, PyObject *args){
    int threadNum = 0, threads = 1;
    PyObject* nodeEdgeVectorObj;
    if (!PyArg_ParseTuple(args, "O|ii", &nodeEdgeVectorObj, &threadNum, &threads)) {
        PyErr_SetString(PyExc_RuntimeError, "buildWeightGraph incorrect arguments");
        return NULL;
    }

    auto nodeEdgeVector = static_cast<std::vector<std::unordered_map<int, std::list<int>>>*>(PyCapsule_GetPointer(nodeEdgeVectorObj, "nodeEdgeVectorPtr"));
    int nodeCount = static_cast<int>(nodeEdgeVector->size());
    //nodeCount = 30; //reduce run time for debug
    auto weightGraph = new std::list<WeightEdge>();
    PyObject* weightGraphCapsule = PyCapsule_New((void*)weightGraph, "weightGraphPtr", NULL);
    PyCapsule_SetPointer(weightGraphCapsule, (void*)weightGraph);

    Py_BEGIN_ALLOW_THREADS //Release global interpreter lock
    for(int i = threadNum; i < nodeCount; i += threads){
        for(int j = i + 1; j < nodeCount; ++j){
            weightGraph->push_back(WeightEdge(i, j, weight(i, j, *nodeEdgeVector)));
        }
        if(threadNum == 0){
            std::cout << "\rThread 0 " << ((double)i / nodeCount) * 100 << "% complete" << std::flush;
        }
    }
    if(threadNum == 0){
        std::cout << std::endl;
    }

    Py_END_ALLOW_THREADS //Reaquire global interpreter lock
    return weightGraphCapsule;
}

static PyObject* mergeWeightGraphs(PyObject *self, PyObject *args){
    PyObject* graphs;
    if (!PyArg_ParseTuple(args, "O", &graphs)) {
        PyErr_SetString(PyExc_RuntimeError, "mergeWeightGraphs incorrect arguments");
        return NULL;
    }

    PyObject* sumWeightGraph = PyList_New(0);
    for(int i = 0; i < PyList_GET_SIZE(graphs); ++i){
        auto subgraph = static_cast<std::list<WeightEdge>*>(PyCapsule_GetPointer(PyList_GET_ITEM(graphs, i), "weightGraphPtr"));
        for(auto j = subgraph->begin(); j != subgraph->end(); ++j){
            PyList_Append(sumWeightGraph, Py_BuildValue("(iid)", j->n1, j->n2, j->weight));
        }
    }

    return sumWeightGraph;
}

static PyObject* buildCatagoryGraph(PyObject *self, PyObject* args){
    PyObject *catagoryNodes;
    Py_buffer originalEdgeArray;
    int nodeCount, edgeCount;
    if (!PyArg_ParseTuple(args, "Oy*ii", &catagoryNodes, &originalEdgeArray, &edgeCount, &nodeCount)) {
        PyErr_SetString(PyExc_RuntimeError, "buildCatagoryGraph incorrect arguments");
        return NULL;
    }

    //Create the vector of nodes that each have a list of clusters they are a part of
    auto nodeToCatagory = std::vector<std::list<int>>();
    int catagoryCount = static_cast<int>(PyList_Size(catagoryNodes));
    std::generate_n(std::inserter(nodeToCatagory, nodeToCatagory.end()), nodeCount, [](){
        return std::list<int>();
    });
    for(int i = 0; i < catagoryCount; ++i){
        PyObject *catNodeList = PyList_GetItem(catagoryNodes, i);
        for(int j = 0; j < PyList_Size(catNodeList); ++j){
            nodeToCatagory[PyLong_AsLong(PyList_GetItem(catNodeList, j))].push_back(i);
        }
    }

    //For each node, find clusters it shares a symptom edge with
    //Vector of nodes with map of symptoms with a set of clusteres
    auto nodeSymptomCatagory = std::vector<std::unordered_map<int, std::set<int>>>();
    int* origionalEdges = static_cast<int*>(originalEdgeArray.buf);
    std::generate_n(std::inserter(nodeSymptomCatagory, nodeSymptomCatagory.end()), nodeCount, [](){
        return std::unordered_map<int, std::set<int>>();
    });
    for(int i = 0; i < edgeCount; ++i){
        fillNodeSymptomCatagory(origionalEdges[i * 3 + 2], nodeSymptomCatagory[origionalEdges[i * 3]], nodeToCatagory[origionalEdges[i * 3 + 1]]);
        fillNodeSymptomCatagory(origionalEdges[i * 3 + 2], nodeSymptomCatagory[origionalEdges[i * 3 + 1]], nodeToCatagory[origionalEdges[i * 3]]);
    }

    //Count the edges between clusters as vector of clusters each with a map of symptom to number of occurances, diretional
    //Vector of clusters with map of symptoms with map of clusters with number of (unique nodes with edges with symptom to the cluster)
    auto catagoryEdgeCounts = std::vector<std::unordered_map<int, std::unordered_map<int, int>>>();
    std::generate_n(std::inserter(catagoryEdgeCounts, catagoryEdgeCounts.end()), catagoryCount, [](){
        return std::unordered_map<int, std::unordered_map<int, int>>();
    });
    //For each node
    for(int i = 0; i < nodeCount; ++i){
        //For each cluster that node is in
        for(auto iCat = nodeToCatagory[i].begin(); iCat != nodeToCatagory[i].end(); ++iCat){
            //For each symptom that node has an edge with
            for(auto iSymp = nodeSymptomCatagory[i].begin(); iSymp != nodeSymptomCatagory[i].end(); ++iSymp){
                catagoryEdgeCounts[*iCat].try_emplace(iSymp->first, std::unordered_map<int, int>());
                //For each cluster the node has an edge of said symptom to
                for(auto iDSympCat = iSymp->second.begin(); iDSympCat != iSymp->second.end(); ++iDSympCat){
                    catagoryEdgeCounts[*iCat][iSymp->first].try_emplace(*iDSympCat, 0);
                    catagoryEdgeCounts[*iCat][iSymp->first][*iDSympCat] += 1;
                }
            }
        }
    }

    //Calculate weights and convert back to python friendly format
    PyObject* catagoryEdges = PyList_New(0);
    for(int i = 0; i < catagoryCount; ++i){
        double iNodeCount = static_cast<double>(PyList_Size(PyList_GetItem(catagoryNodes, i)));
        for(auto iSymp = catagoryEdgeCounts[i].begin(); iSymp != catagoryEdgeCounts[i].end(); ++iSymp){
            for(auto desCat = iSymp->second.begin(); desCat != iSymp->second.end(); ++desCat){
                PyList_Append(catagoryEdges, Py_BuildValue("(iiid)", i, desCat->first, iSymp->first, (double)desCat->second / iNodeCount));
            }
        }
    }

    PyBuffer_Release(&originalEdgeArray);
    return catagoryEdges;
}

static PyObject* preprocessCatagories(PyObject* self, PyObject* args){
    PyObject *catagoryNodes, *catagoryEdges;
    if (!PyArg_ParseTuple(args, "OO", &catagoryNodes, &catagoryEdges)) {
        PyErr_SetString(PyExc_RuntimeError, "preprocessCatagories incorrect arguments");
        return NULL;
    }

    //Vector of nodes containing list of catagories
    int catCount = static_cast<int>(PyList_Size(catagoryNodes));
    auto nodes = new std::unordered_map<int, std::list<int>>();
    for(int i = 0; i < catCount; ++i){
        PyObject* catNodes = PyList_GetItem(catagoryNodes, i);
        for(int j = 0; j < PyList_Size(catNodes); ++j){
            int node = static_cast<int>(PyLong_AsLong(PyList_GetItem(catNodes, j)));
            nodes->try_emplace(node, std::list<int>());
            (*nodes)[node].push_back(i);
        }
    }
    PyObject* capCatNodes = PyCapsule_New((void*)nodes, "catNodesPtr", NULL);
    PyCapsule_SetPointer(capCatNodes, (void*)nodes);

    //Vector of catagories containing map of dest catagories containing map of symptoms to weight
    auto edges = new std::vector<std::unordered_map<int, std::unordered_map<int, double>>>();
    std::generate_n(std::inserter(*edges, edges->end()), catCount, [](){
        return std::unordered_map<int, std::unordered_map<int, double>>();
    });
    int edgeCount = static_cast<int>(PyList_Size(catagoryEdges));
    for(int i = 0; i < edgeCount; ++i){
        int source, destination, symptom;
        double weight;
        PyArg_ParseTuple(PyList_GetItem(catagoryEdges, i), "iiid", &source, &destination, &symptom, &weight);
        (*edges)[source].try_emplace(destination, std::unordered_map<int, double>());
        (*edges)[source][destination][symptom] = weight;
    }
    PyObject* capCatEdges = PyCapsule_New((void*)edges, "catEdgesPtr", NULL);
    PyCapsule_SetPointer(capCatEdges, (void*)edges);

    return Py_BuildValue("OO", capCatNodes, capCatEdges);
}

static PyObject* predictSymptoms(PyObject* self, PyObject* args){
    PyObject *newEdges, *preNodes, *preEdges, *catagoryNodes;
    double minWeight;
    if (!PyArg_ParseTuple(args, "OOOOd", &newEdges, &preNodes, &preEdges, &catagoryNodes, &minWeight)) {
        PyErr_SetString(PyExc_RuntimeError, "predictSymptoms incorrect arguments");
        return NULL;
    }
    auto catNodes = static_cast<std::unordered_map<int, std::list<int>>*>(PyCapsule_GetPointer(preNodes, "catNodesPtr"));
    auto catEdges = static_cast<std::vector<std::unordered_map<int, std::unordered_map<int, double>>>*>(PyCapsule_GetPointer(preEdges, "catEdgesPtr"));

    //Convert the drug to drug edges to drug to cluster
    auto clusterEdges = std::set<std::pair<int,int>>();
    int newEdgeCount = static_cast<int>(PyList_Size(newEdges));
    for(int i = 0; i < newEdgeCount; ++i){
        int src, dest, symp;
        PyArg_ParseTuple(PyList_GetItem(newEdges, i), "iii", &src, &dest, &symp);
        for(auto cat = (*catNodes)[dest].begin(); cat != (*catNodes)[dest].end(); ++cat){
            clusterEdges.insert(std::make_pair(*cat, symp));
        }
    }

    //Find potental symptoms
    PyObject* predictedEdges = PyList_New(0);
    int newClusterEdgeCount = static_cast<int>(clusterEdges.size());
    int catCount = static_cast<int>(PyList_Size(catagoryNodes));
    //For each catagory
    for(int cat = 0; cat < catCount; ++cat){
        //Get number of outgoing edges in catagory
        int catEdgeCount = std::accumulate((*catEdges)[cat].begin(), (*catEdges)[cat].end(), 0, [](const int& acc, const auto& dest){
            return acc + static_cast<int>(dest.second.size());
        });
        //Count the similar edges and calculate drug to catagory similarity weight
        int similarEdgeCount = 0;
        for(auto edge = clusterEdges.begin(); edge != clusterEdges.end(); ++edge){
            if(edge->first == cat || (*catEdges)[cat].find(edge->first) != (*catEdges)[cat].end() && (*catEdges)[cat][edge->first].find(edge->second) != (*catEdges)[cat][edge->first].end()){
                similarEdgeCount += 2;
            }
        }
        double weight = (double)similarEdgeCount / (double)(catEdgeCount + newClusterEdgeCount);
        if(weight >= minWeight){
            //For each outgoing edge, see it is aboe the min weight and add it to prediction
            for(auto destCat = (*catEdges)[cat].begin(); destCat != (*catEdges)[cat].end(); ++destCat){
                for(auto catSymp = destCat->second.begin(); catSymp != destCat->second.end(); ++catSymp){
                    if(catSymp->second * weight >= minWeight && clusterEdges.find(std::make_pair(destCat->first, catSymp->first)) == clusterEdges.end()){
                        PyList_Append(predictedEdges, Py_BuildValue("(iid)", destCat->first, catSymp->first, catSymp->second * weight));
                    }
                }
            }
        }
    }

    return predictedEdges;
}

static PyObject* findEdges(PyObject* self, PyObject* args){
    PyObject *edges;
    int node;
    if (!PyArg_ParseTuple(args, "iO", &node, &edges)) {
        PyErr_SetString(PyExc_RuntimeError, "findEdges incorrect arguments");
        return NULL;
    }

    int edgeCount = static_cast<int>(PyList_Size(edges));
    PyObject* nodeEdges = PyList_New(0);
    for(int i = 0; i < edgeCount; ++i){
        int sour, dest, symp;
        if (!PyArg_ParseTuple(PyList_GetItem(edges, i), "iii", &sour, &dest, &symp)){
            PyErr_SetString(PyExc_RuntimeError, "findEdges edge list incorrect");
            return NULL;
        }
        if(sour == node || dest == node){
            PyList_Append(nodeEdges, Py_BuildValue("(iii)", sour, dest, symp));
        }
    }

    return nodeEdges;
}

//Module defs
static PyMethodDef accel_methods[] = {
    {"buildNodeEdgeVector",  buildNodeEdgeVector, METH_VARARGS,
     "Build a vector of dictionaries of edges keyed by symptom \n(nodeCount, edgeCount, edges)"},
    {"buildWeightGraph",  buildWeightGraph, METH_VARARGS,
     "Build the weight graph, supports theading \n(nodeEdgeVector, threadNum=0, threads=1)"},
    {"mergeWeightGraphs",  mergeWeightGraphs, METH_VARARGS,
     "Merge list of weight subgraphs \n(graphs)"},
    {"buildCatagoryGraph",  buildCatagoryGraph, METH_VARARGS,
     "Create the catagory graph \n(catagoryNodes, originalEdgeArray, edgeCount, nodeCount)"},
    {"preprocessCatagories",  preprocessCatagories, METH_VARARGS,
     "Convert the catagory node and edge list into a faster c format \n(catagoryNodes, catagoryEdges)"},
    {"predictSymptoms",  predictSymptoms, METH_VARARGS,
     "Preditct the interdrug interaction symptoms of a partial graph \n(newEdges, preprocessedNodes, preprocessedEdges, catagoryNodes, minWeight)"},
    {"findEdges",  findEdges, METH_VARARGS,
     "Returns list of edges associated with node \n(node, edges)"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef accel_module = {
    PyModuleDef_HEAD_INIT,
    "accel",   /* name of module */
    "Accerate the graph computionations\n", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    accel_methods
};

PyMODINIT_FUNC
PyInit_accel(void) {
   return PyModule_Create(&accel_module);
}