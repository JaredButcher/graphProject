import snap
import sys
import os
import sideEffects.accel as accel
import time
import concurrent.futures
import threading
import ctypes
import itertools
import random
import math

def loadBaseGraph(snapFilename='./ChChSe-Decagon_polypharmacy.csv', nodeListfilename = './nodelist.csv', edgefilename = './edges.csv', edgeListfilename = './edgelist.csv', rebuild=False):
    '''Loads in the needed files, generatres more useable formats if they don't exist
    '''

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
        with open(snapFilename, 'r') as snapFile:
            for entry in snapFile:
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

def listOfListsToCArray(iterable, subListSize):
    return (ctypes.c_int32 * (subListSize * len(iterable)))(*itertools.chain(*iterable))

#Will proably take about  hours to run
def genWeightedGraph(edges, nodeList, threads = 1, weightedGraphFilename = "weighted.csv", rebuild=False):
    '''Generates a complete graph of the drug nodes and the edges storing weights from 0-1 
        representing how similar their common symptom interactions are.
        Will take about 30 min / threads to run with on a desktop processor
    '''
    weightedGraph = []

    if rebuild or not os.path.exists(weightedGraphFilename):
        nodeCount = len(nodeList)
        edgeCount = len(edges)

        print("Converting edges to c array")
        arrayEdges = listOfListsToCArray(edges, 3)
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
        with open(weightedGraphFilename, "w") as weightFile:
            for entry in weightedGraph:
                weightFile.write(f'{str(entry[0])},{str(entry[1])},{str(entry[2])}\n')
        print("Complete")
    else:
        with open(weightedGraphFilename, "r") as weightFile:
            for edge in weightFile:
                values = edge.split(',')
                weightedGraph.append((int(values[0]),int(values[1]),float(values[2])))

    return weightedGraph

def clusterWeightedGraph(weightedGraph, edges, nodelist, edgelist, weightThreashold = 0.2, catagoryNodeFilename='catagorynodes.csv', catagoryEdgeFilename='catagoryedges.csv', rebuild=False):
    '''Use weighted graph to cluster nodes and generate a new catagorization graph.
        Nodes with the weight between them close to 1 should be clustered together. 
        A single node may be part of multiple clusters. 
        The nodes on this catagorization graph represent the clusters in the weighted graph.
        The catagorization graph will have multiple edges between nodes.
        These edges are directed and will store a symptom and a chariteristic weight.
        For an edge AB, the weight is the number of nodes in A that have an edge of the symptom to B divided by the number of nodes in A.
    '''
    if rebuild or not os.path.exists(catagoryNodeFilename) or not os.path.exists(catagoryEdgeFilename):
        print("Purgeing edges between weakly related nodes")
        weightedGraph = [x for x in weightedGraph if x[2] > weightThreashold]

        print("Building snap graph")
        snapWeightedGraph = snap.TUNGraph_New(len(nodelist), len(weightedGraph))
        for i in range(0, len(nodelist)):
            snapWeightedGraph.AddNode(i)
        for edge in weightedGraph:
            snapWeightedGraph.AddEdge(edge[0], edge[1])
        print("Clustering weighted graph")
        catagoryNodes = snap.TCnComV()
        print(f'Mod: {snap.CommunityCNM(snapWeightedGraph, catagoryNodes)}\nComponet sizes: ')
        for comp in catagoryNodes:
            if comp.Len() > 2:
                print(f'{len(comp)},', end='')
        print(f'\ntotal {len(catagoryNodes)}')
        
        print("Converting origional edges to c array")
        arrayOrigionalEdges = listOfListsToCArray(edges, 3)
        print("Edge array constructred")

        pyCatagoryNodes = []

        for cat in catagoryNodes:
            temp = []
            for node in cat:
                temp.append(node)
            pyCatagoryNodes.append(temp)

        #Make new graph, nodes corisponding to clusters in weightedGraph
        ##Made edges in this new graph containing a symptom and weight
        ##Nodes will have a list of nodes from base graph that are contained within the cluster
        ##for edge AB, the weight is the number of nodes in A that have an edge of the symptom to B divided by the number of nodes in A.
        print("Buiding catagory edges")
        catagoryEdges = accel.buildCatagoryGraph(pyCatagoryNodes, arrayOrigionalEdges, len(edges), len(nodelist))
        print("Catagory edges built")

        #Save and return
        print("Saving catagory nodes")
        with open(catagoryNodeFilename, 'w') as nodeFile:
            for catagory in catagoryNodes:
                tempcat = [node for node in catagory]
                for node in tempcat[:-1]:
                    nodeFile.write(f'{str(node)},')
                nodeFile.write(f'{str(tempcat[-1])}\n')
        
        print("Saving catagory edges")
        with open(catagoryEdgeFilename, 'w') as edgeFile:
            for edge in catagoryEdges:
                edgeFile.write(f'{str(edge[0])},{str(edge[1])},{str(edge[2])},{str(edge[3])}\n')
    else:
        print("Loading catagory nodes")
        catagoryNodes = []
        with open(catagoryNodeFilename, 'r') as nodeFile:
            for entry in nodeFile:
                nodes = entry.split(',')
                catagoryNodes.append([int(node) for node in nodes])

        print("Loading catagory edges")
        catagoryEdges = []
        with open(catagoryEdgeFilename, 'r') as edgefile:
            for entry in edgefile:
                values = entry.split(',')
                catagoryEdges.append((int(values[0]), int(values[1]), int(values[2]), float(values[3])))

    # Store catagory edges as normal
    return (catagoryNodes, catagoryEdges)

def predictSymptoms(symptoms, catagoryNodes, catagoryEdges, preprocessedNodes=None, preprocessedEdges=None, cutoff=.2):
    '''Symptoms in [(drug index, symptom index)] format
    '''
    if not preprocessedEdges or not preprocessedNodes:
        print("Preprocessing")
        preprocessedNodes, preprocessedEdges = accel.preprocessCatagories(catagoryNodes, catagoryEdges)
    
    print("Predicting symptoms")
    return accel.predictSymptoms(symptoms, preprocessedNodes, preprocessedEdges, cutoff)

def testPrediction(edges, nodeList, catagoryNodes, catagoryEdges, nodesToUse=None, randomNodesToUse=10, edgeFraction=.5, cutoff=.2):
    print("Finding random nodes")
    nodes = []
    trueEdgeSets = []
    while len(nodes) < randomNodesToUse:
        node = random.randrange(0, len(nodeList))
        if not node in nodes:
            nodeEdges = accel.findEdges(node, edges)
            print(len(nodeEdges))
            if(len(nodeEdges) > 20):
                nodes.append(node)
                trueEdgeSets.append(nodeEdges)

    print("Preprocessing")
    preprocessedNodes, preprocessedEdges = accel.preprocessCatagories(catagoryNodes, catagoryEdges)
    
    print("Predicting")
    usedEdgeSets = []
    predictedEdgeSets = []
    for edgeSet in trueEdgeSets:
        reducedSet = edgeSet[:math.floor(len(edgeSet) * edgeFraction)]
        usedEdgeSets.append(reducedSet)
        predictedSet = accel.predictSymptoms(reducedSet, preprocessedNodes, preprocessedEdges, catagoryNodes, cutoff)
        predictedEdgeSets.append(predictedSet)
        print(f'True edge set size: {len(edgeSet)}')
        print(f'Reduced edge set size: {len(reducedSet)}')
        print(f'Predicted edge set size: {len(predictedSet)}')
        print(f'True edge set size divided by the sum of the Reduced and Predicted edge set size: {len(edgeSet) / (len(predictedSet) + len(reducedSet))}')
    print("\nFinished")
    

def main():
    edges, nodeList, edgeList = loadBaseGraph()
    weightedGraph = genWeightedGraph(edges, nodeList, 8)
    catagoryNodes, catagoryEdges = clusterWeightedGraph(weightedGraph, edges, nodeList, edgeList, 0.45)
    testPrediction(edges, nodeList, catagoryNodes, catagoryEdges, randomNodesToUse=3, edgeFraction=.5, cutoff=.2)
    #TODO Create a function that compares a new graph's node to the catagory graph
    ## Similar to current weight function, run between new node and each catagory
    ## Predicted symptoms are weighted by both its weight in catagory graph and the drugs weight to the catagory