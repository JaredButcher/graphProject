import snap
import sys
import os
import sideEffects.accel as accel
import time
import concurrent.futures
import threading
import ctypes
import itertools

def loadBaseGraph(rebuild=False):
    '''Loads in the needed files, generatres more useable formats if they don't exist
    '''
    filename = './ChChSe-Decagon_polypharmacy.csv'
    nodeListfilename = './nodelist.csv'
    edgefilename = './edges.csv'
    edgeListfilename = './edgelist.csv'

    nodeList = []
    edgeList = []
    edges = []

    nodeCount = 0
    edgeCount = 0
    edgeValueCount = 0

    if rebuild or not os.path.exists(nodeListfilename) or not os.path.exists(edgeListfilename) or not os.path.exists(edgefilename):
        print("Building lists")
        nodeDict = {}
        edgeDict = {}
        for entry in open(filename, 'r'):
            if entry[0] != '#': #column names
                values = entry.split(',')
                if not values[0] in nodeDict:
                    nodeList.append(values[0])
                    nodeDict[values[0]] = nodeCount
                    nodeCount += 1
                if not values[1] in nodeDict:
                    nodeList.append(values[1])
                    nodeDict[values[1]] = nodeCount
                    nodeCount += 1
                if not values[2] in edgeDict:
                    edgeList.append((values[2], values[3][:-1]))
                    edgeDict[values[2]] = edgeValueCount
                    edgeValueCount += 1
                edges.append((nodeDict[values[0]], nodeDict[values[1]], edgeDict[values[2]]))
                edgeCount += 1
            if not edgeCount % 100000:
                sys.stdout.write(f'\r{edgeCount} edges')

        print("\nSaving node values")
        with open(nodeListfilename, 'w') as nodefile:
            for entry in nodeList:
                nodefile.write(f'{str(entry)}\n')

        print("Saving edge values")
        with open(edgeListfilename, 'w') as edgeValueFile:
            for entry in edgeList:
                edgeValueFile.write(f'{str(entry[0])},{entry[1]}\n')

        print("Saving edges")
        with open(edgefilename, 'w') as edgefile:
            for entry in edges:
                edgefile.write(f'{str(entry[0])},{str(entry[1])},{str(entry[2])}\n')
    else:
        print("Loading edges")
        for entry in open(edgefilename, 'r'):
            values = entry.split(',')
            edges.append((int(values[0]), int(values[1]), int(values[2])))
            edgeCount += 1
            if not edgeCount % 100000:
                sys.stdout.write(f'\r{edgeCount} edges')
        
        print("\nLoading node values")
        for entry in open(nodeListfilename, 'r'):
            nodeList.append(entry)

        print("Loading edge values")
        for entry in open(edgeListfilename, 'r'):
            edgeList.append(entry)

    return (edges, nodeList, edgeList)

#Will proably take about  hours to run
def genWeightedGraph(edges, nodeList, threads = 1):
    '''Generates a complete graph of the drug nodes and the edges storing weights from 0-1 
        representing how similar their common symptom interactions are.
        Will take about 30 min / threads to run with on a desktop processor
    '''
    weightedGraphFilename = "weighted.csv"

    nodeCount = len(nodeList)
    edgeCount = len(edges)

    print("Converting edges to c array")
    arrayEdges = (ctypes.c_int32 * (3 * edgeCount))(*itertools.chain(*edges))
    print("Edge array constructred")

    print("Constructing nodeEdgeVector")
    nodeEdgeVector = accel.buildNodeEdgeVector(nodeCount, edgeCount, arrayEdges)
    print("NodeEdgeVector constructed, building weighted graph")

    subGraphs = []
    subGraphLock = threading.Lock()

    def runBuildWeightGraph(threadNumber):
        print(f'Starting thread {threadNumber}')
        weightedGraph = accel.buildWeightGraph(nodeEdgeVector, threadNumber, threads)
        with subGraphLock:
            print(f'Thread {threadNumber} complete')
            subGraphs.append(weightedGraph)

    with concurrent.futures.ThreadPoolExecutor(threads) as executor:
        print(f'Starting {threads} threads')
        for result in executor.map(runBuildWeightGraph, range(0,threads)):
            pass

    print("Converting into python friendly format")
    weightedGraph = accel.mergeWeightGraphs(subGraphs)
    print("Complete")

    print("Weighted Graph built, writeing to file")
    with open("weight.csv", "w") as weightFile:
        for entry in weightedGraph:
            weightFile.write(f'{str(entry[0])},{str(entry[1])},{str(entry[2])}\n')
    print("Complete")

    return weightedGraph

def clusterWeightedGraph(weightedGraph, edges, nodelist, edgelist):
    '''Use weighted graph to cluster nodes and generate a new catagorization graph.
        Nodes with the weight between them close to 1 should be clustered together. 
        A single node may be part of multiple clusters. 
        The nodes on this catagorization graph represent the clusters in the weighted graph.
        The catagorization graph will have multiple edges between nodes.
        These edges are directed and will store a symptom and a chariteristic weight.
        For an edge AB, the weight is the number of nodes in A that have an edge of the symptom to B divided by the number of nodes in A.
    '''
    #Perge edges in weightedGraph below threashold
    #Cluster weightedGraph
    #Make new graph, nodes corisponding to clusters in weightedGraph
    #Made edges in this new graph containing a symptom and weight
    ##Nodes will have a list of nodes from base graph that are contained within the cluster
    ##for edge AB, the weight is the number of nodes in A that have an edge of the symptom to B divided by the number of nodes in A.
    #Save and return
    return None

def main():
    edges, nodeList, edgeList = loadBaseGraph()
    weightedGraph = genWeightedGraph(edges, nodeList, 10)
    catagorizationGraph = clusterWeightedGraph(weightedGraph, edges, nodeList, edgeList)