#include <Python.h>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <iterator>

//About 10% faster
static PyObject *edgeComp(PyObject *self, PyObject *args){
    PyObject *edge1, *edge2;
    PyArg_ParseTuple(args, "OO", &edge1, &edge2);

    int e1_node1, e1_node2;
    char* e1_sideEffect, trash;
    PyArg_ParseTuple(edge1, "iiss", &e1_node1, &e1_node2, &e1_sideEffect, &trash);

    int e2_node1, e2_node2;
    char* e2_sideEffect;
    PyArg_ParseTuple(edge2, "iiss", &e2_node1, &e2_node2, &e2_sideEffect, &trash);

    return Py_BuildValue("O", e1_node1 == e2_node1 || e1_node1 == e2_node2 || e1_node2 == e2_node1 || e1_node2 == e2_node2 ? Py_True : Py_False);
}

//Currently using incorrect method
static PyObject *weight(PyObject *self, PyObject *args){
    int node1, node2;
    PyObject* nodeEdgeVector;
    if (!PyArg_ParseTuple(args, "iiO", &node1, &node2, &nodeEdgeVector)) {
      return NULL;
    }

    PyObject* node1Edges = PyList_GetItem(nodeEdgeVector, node1);
    PyObject* node2Edges = PyList_GetItem(nodeEdgeVector, node2);

    int node1Len = static_cast<int>(PyList_Size(node1Edges));
    int node2Len = static_cast<int>(PyList_Size(node2Edges));

    int weight = 0;

    for(int i = 0; i < node1Len; ++i){
        for(int j = 0; j < node2Len; ++j){
            PyObject* edge1 = PyList_GetItem(node1Edges, i);
            PyObject* edge2 = PyList_GetItem(node2Edges, j);

            int e1_node1, e1_node2;
            char* e1_sideEffect, trash;
            if(!PyArg_ParseTuple(edge1, "iiss", &e1_node1, &e1_node2, &e1_sideEffect, &trash)){
                return NULL;
            }

            int e2_node1, e2_node2;
            char* e2_sideEffect;
            if(!PyArg_ParseTuple(edge2, "iiss", &e2_node1, &e2_node2, &e2_sideEffect, &trash)){
                return NULL;
            }
            if(strcmp(e1_sideEffect, e2_sideEffect) == 0 && (e1_node1 == e2_node1 || e1_node1 == e2_node2 || e1_node2 == e2_node1 || e1_node2 == e2_node2)){
                ++weight;
                break;
            }
        }
    }

    double total = node1Len + node2Len;
    if(total == 0){
        total = 1;
    }

    return Py_BuildValue("d", (double)weight / total);
}

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

//Module defs
static PyMethodDef accel_methods[] = {
    {"weight",  weight, METH_VARARGS,
     "Compute weight of simularity between nodes"},
    {"edgeComp",  edgeComp, METH_VARARGS,
     "Compare edges"},
    {"buildNodeEdgeVector",  buildNodeEdgeVector, METH_VARARGS,
     "Build a vector of dictionaries of edges keyed by symptom \n(nodeCount, edgeCount, edges)"},
    {"buildWeightGraph",  buildWeightGraph, METH_VARARGS,
     "Build the weight graph, supports theading \n(nodeEdgeVector, threadNum=0, threads=1)"},
    {"mergeWeightGraphs",  mergeWeightGraphs, METH_VARARGS,
     "Merge list of weight subgraphs \n(graphs)"},
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