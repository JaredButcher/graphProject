import snap
import sys
import os
import sideEffects.accel as accel
import time
import multiprocessing

def loadBaseGraph(rebuild=False):
    filename = './ChChSe-Decagon_polypharmacy.csv'
    nodefilename = './nodes.csv'
    edgefilename = './edges.csv'

    nodes = []
    edges = []

    nodeCount = 0
    edgeCount = 0

    if rebuild or not os.path.exists(nodefilename) or not os.path.exists(edgefilename):
        print("Building lists")
        for entry in open(filename, 'r'):
            if entry[0] != '#': #column names
                values = entry.split(',')
                if not values[0] in nodes:
                    nodes.append(values[0])
                    nodeCount += 1
                if not values[1] in nodes:
                    nodes.append(values[1])
                    nodeCount += 1
                edges.append((nodes[values[0]], nodes[values[1]], values[2], values[3][:-1]))
                edgeCount += 1
            if not edgeCount % 100000:
                sys.stdout.write(f'\r{edgeCount} edges')

        sys.stdout.write(f'\r{edgeCount} edges')
        print("\nSaving nodes")

        with open(nodefilename, 'w') as nodefile:
            for entry in nodes:
                nodefile.write(f'{str(entry)}\n')

        print("Saving edges")

        with open(edgefilename, 'w') as edgefile:
            for entry in edges:
                edgefile.write(f'{str(entry[0])},{str(entry[1])},{str(entry[2])},{str(entry[3])}\n')
    else:
        print("Loading edges")
        for entry in open(edgefilename, 'r'):
            values = entry.split(',')
            edges.append((int(values[0]), int(values[1]), values[2], values[3]))
            edgeCount += 1
            if not edgeCount % 100000:
                sys.stdout.write(f'\r{edgeCount} edges')
        
        print("\nLoading nodes")
        for entry in open(nodefilename, 'r'):
            nodes.append(entry)
    
    return (edges, nodes)

def buildNodeEdgeVector(nodeCount, edges):
    nodeEdgeVector = [{} for i in range(0, nodeCount)]
    for edge in edges:
        if edge[2] in nodeEdgeVector[edge[0]]:
            nodeEdgeVector[edge[0]][edge[2]].append(edge)
        else:
            nodeEdgeVector[edge[0]][edge[2]] = [edge]
        if edge[2] in nodeEdgeVector[edge[1]]:
            nodeEdgeVector[edge[1]][edge[2]].append(edge)
        else:
            nodeEdgeVector[edge[1]][edge[2]] = [edge]
    return nodeEdgeVector

def edgeComp(edge1, edge2):
    return edge1[0] == edge2[0] or edge1[0] == edge2[1] or edge1[1] == edge2[0] or edge1[1] == edge2[1]

def weight(node1, node2, nodeEdgeVector):
    equalEdges = 0

    sympts = set(nodeEdgeVector[node1].keys()).intersection(nodeEdgeVector[node2].keys())

    for sympt in sympts:
        for edge1 in nodeEdgeVector[node1][sympt]:
            for edge2 in nodeEdgeVector[node2][sympt]:
                if edgeComp(edge1, edge2):
                    equalEdges += 1
                    break
    return equalEdges / (len(nodeEdgeVector[node1]) + len(nodeEdgeVector[node2]))

#Will proably take about 2 hours to run
def genWeightedGraph(processes = 1):
    weightedGraphFilename = "weighted.csv"

    edges, nodes = loadBaseGraph()
    #baseGraph = snap.LoadEdgeList_PUndirNet('edges.csv', 0, 1, ',')
    nodeCount = len(nodes)

    print("Building node edge vector")
    nodeEdgeVector = buildNodeEdgeVector(nodeCount, edges)

    if processes > 2:
        print("Multiprocessing still a work in progress")
    else:
        print("Creating weight graph")
        weightedEdges = []
        currentNode = 0
        for i in range(0,nodeCount):
            for j in range(i + 1,nodeCount):
                weightedEdges.append((i, j, weight(i, j, nodeEdgeVector)))
                sys.stdout.write(f'\r{i * nodeCount + j}/{nodeCount**2} nodes weighted')
        print("\n Weight graph created")

        with open(weightedGraphFilename, 'w') as weightedGraphFile:
            for edge in weightedEdges:
                weightedGraphFile.write(f'{str(weightedEdges[0])},{str(weightedEdges[1])},{str(weightedEdges[2])}\n')
        print("File saved")

def main():
    genWeightedGraph()