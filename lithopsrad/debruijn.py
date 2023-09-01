import os
import sys

import lithopsrad.sequence as seq

class deBruijn:
    def __init__(self, reads, k):
        self.nodes = self.construct_graph(reads, k)

    def construct_graph(self, reads, k):
        nodes={}
        for depth, read in reads:
            for left, right in seq.get_kmers_neighbor(read, k):
                if left not in nodes:
                    nodes[left] = Node(left)
                if right not in nodes:
                    nodes[right] = Node(right)

                if right not in nodes[left].edges:
                    nodes[left].edges[right] = Edge(nodes[right])

                nodes[left].edges[right].weight += int(depth)
                nodes[right].nin += 1
                nodes[left].nout += 1
        return(nodes)

    def get_consensus(self):
        nodes = self.nodes
        start = list(nodes.values())[0]
        for k in nodes.keys():
            if nodes[k].nin < start.nin:
                start = nodes[k]

        consensus = start.label
        current = start
        visited=set()


        while len(current.edges) > 0:
            visited.add(current.label)
            next=list(current.edges.keys())[0]
            for n in current.edges.keys():
                if current.edges[n].weight > current.edges[next].weight:
                    next = n
            consensus += current.edges[next].node.label[-1]
            current.edges=dict()
            if next in visited:
                return(None)
            else:
                current = nodes[next]
        return consensus

class Edge:
    def __init__(self, node, weight=0):
        self.node = node
        self.weight = weight

class Node:
    def __init__(self, label):
        self.label = label
        self.nin = 0
        self.nout = 0
        self.edges = dict()
