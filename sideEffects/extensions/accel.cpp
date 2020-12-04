#include <Python.h>
#include <algorithm>
#include <iostream>
#include <cstring>

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

    int node1Len = PyList_Size(node1Edges);
    int node2Len = PyList_Size(node2Edges);

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

//Module defs
static PyMethodDef accel_methods[] = {
    {"weight",  weight, METH_VARARGS,
     "Compute weight of simularity between nodes"},
    {"edgeComp",  edgeComp, METH_VARARGS,
     "Compare edges"},
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