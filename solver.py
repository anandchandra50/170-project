import networkx as nx
from parse import read_input_file, write_output_file
from utils import is_valid_network, average_pairwise_distance
from os import listdir
from os.path import isfile, join
import sys

# return a list of tuple edges connected to a node
# def get_edgelist(edges_by_node, node):
#     all_edges = []
#     for other_vertex in edges_by_node[node]:
#         all_edges.append((node, other_vertex))
#     return all_edges
#
# def get_connecting_edges(nodes, new_node, edges_by_node):
#     current_nodes_that_connect_to_new = set(edges_by_node[new_node]).intersection(nodes)
#     all_edges = []
#     for current in current_nodes_that_connect_to_new:
#         all_edges.append((new_node, current))
#     return all_edges
def validate_edge(edge, used_nodes):
    return (edge[0] in used_nodes and edge[1] not in used_nodes) or (edge[1] in used_nodes and edge[0] not in used_nodes)

def basic_pairwise_distance(G):
    node_pairwise_distance = {} # node number to total distance from all other nodes
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G))
    for node in G.nodes:
        node_pairwise_distance[node] = sum(path_lengths[node][n] for n in G.nodes)
    node_pairwise_distance = {k: v for k, v in sorted(node_pairwise_distance.items(), key=lambda item: item[1])}
    return node_pairwise_distance

def modified_pairwise_distance(G, T):
    # CHANGE THIS TO INCLUDE T
    node_pairwise_distance = {} # node number to total distance from all other nodes
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G))
    for node in G.nodes:
        node_pairwise_distance[node] = sum(path_lengths[node][n] for n in G.nodes)
    node_pairwise_distance = {k: v for k, v in sorted(node_pairwise_distance.items(), key=lambda item: item[1])}
    return node_pairwise_distance

def update_available_edges(T, available_edges):
    for i in range(len(available_edges)):
        (u, v, w) = available_edges[i]
        new_vertex = u
        if u in list(T.nodes):
            new_vertex = v
        T_copy = nx.Graph()
        T_copy.add_nodes_from(T.nodes)
        T_copy.add_edges_from(T.edges())
        T_copy.add_edge(u, v, weight=w)
        pairwise_distance = basic_pairwise_distance(T_copy)
        new_cost = pairwise_distance[new_vertex]
        available_edges[i] = (u, v, w, new_cost)
    available_edges.sort(key=lambda item: item[3])
    return available_edges

# create MST
def prims(G):
    # init graph
    T = nx.Graph()

    # node distances
    node_pairwise_distance = basic_pairwise_distance(G)

    # edges
    all_edges = list(G.edges.data('weight'))
    # all_edges.sort(key=lambda item: (item[2], node_pairwise_distance[item[0]] + node_pairwise_distance[item[1]]))
    all_edges.sort(key=lambda item: (node_pairwise_distance[item[0]] + node_pairwise_distance[item[1]], item[2]))

    starting_node = list(node_pairwise_distance.keys())[0]
    T.add_node(starting_node)

    # loop
    while len(list(T.nodes)) != len(list(G.nodes)):
        # filter available edges: include only if exactly one of the vertices is in used_nodes
        # if len(list(T.nodes)) == int(0.5 * len(list(G.nodes))): # change threshold
            # all_edges = update_available_edge(G, T, all_edges) # update the sorting of all edges
        available_edges = [all_edges[i] for i in range(len(all_edges)) if validate_edge(all_edges[i], list(T.nodes))]
        if len(list(T.nodes)) >= int(0.5 * len(list(G.nodes))): # change threshold
            available_edges = update_available_edges(T, available_edges)

        # already sorted by weight, so add the next edge
        next_edge = available_edges[0]
        T.add_edge(next_edge[0], next_edge[1], weight=next_edge[2])
    return T


def solve(G):
    """
    Args:
        G: networkx.Graph

    Returns:
        T: networkx.Graph
    """
    T = prims(G)
    # T = nx.minimum_spanning_tree(G)

    if len(T.nodes) <= 2:
        return T

    leaves = []
    # remove leaves first
    for node in T.nodes:
        if T.degree(node) == 1:
            leaves.append(node)

    # remove these leaves
    for node in leaves:
        T.remove_node(node)


    def get_node_distances(G):
        """
        Given a graph G, returns a dictionary that maps each node in G to the sum of
        distances from each node to every other node
        """
        node_pairwise_distance = {} # node number to total distance from all other nodes
        path_lengths = dict(nx.all_pairs_dijkstra_path_length(G))
        for node in G.nodes:
            node_pairwise_distance[node] = sum(path_lengths[node][n] for n in G.nodes)
        node_pairwise_distance = {k: v for k, v in sorted(node_pairwise_distance.items(), key=lambda item: item[1], reverse=True)}
        return node_pairwise_distance

    change = True
    while change:
        # get distances
        node_pairwise_distance = get_node_distances(T)
        change = False
        for node in node_pairwise_distance: # descending order, remove nodes with highest distances when possible
            # makes graph copy
            t_copy = nx.Graph()
            t_copy.add_nodes_from(T.nodes)
            t_copy.add_edges_from(T.edges)
            t_copy.remove_node(node)
            if len(list(t_copy.nodes)) == 0:
                continue
            if is_valid_network(G, t_copy):
                # if we modify T, note that we made a change and recompute distances
                T.remove_node(node)
                change = True
                break

    return T
    # TODO: your code here!
    # create T
    # T = nx.Graph()
    #
    # edges_by_node = {} # node to list of other nodes (goes both ways), represents edge matrix.
    # for node in G.nodes:
    #     edges_by_node[node] = []
    # for edge in G.edges:
    #     edges_by_node[edge[0]].append(edge[1])
    #     edges_by_node[edge[1]].append(edge[0])
    #
    # node_pairwise_distance = {} # node number to total distance from all other nodes
    #
    # path_lengths = dict(nx.all_pairs_dijkstra_path_length(G))
    # for node in G.nodes:
    #     node_pairwise_distance[node] = sum(path_lengths[node][n] for n in G.nodes)
    # node_pairwise_distance = {k: v for k, v in sorted(node_pairwise_distance.items(), key=lambda item: item[1])}
    # relevant_nodes = [] # (node, edge) keep track of neighboring nodes that we can add (ensure we build a tree)
    # print(node_pairwise_distance)
    # # begin the algorithm
    # first_node = list(node_pairwise_distance.keys())[0]
    # T.add_node(first_node)
    # relevant_nodes = relevant_nodes + list(G.neighbors(first_node))
    # relevant_nodes.sort(key=lambda item: node_pairwise_distance[item])
    # print(G.nodes)
    # print(G.edges)
    # print('starting: {}'.format(T.nodes))
    # # now loop
    # while True:
    #     print('running {}'.format(relevant_nodes))
    #     next_node = relevant_nodes.pop(0)
    #     # T.add_node(next_node)
    #     new_edges = get_connecting_edges(T.nodes, next_node, edges_by_node)
    #     print(new_edges[0])
    #     T.add_edge(*new_edges[0])
    #     print('nodes: {}'.format(T.nodes))
    #     if not nx.is_tree(T):
    #         print('no longer tree')
    #         T.remove_node(next_node)
    #         if len(relevant_nodes) == 0:
    #             return T
    #         relevant_nodes.append(next_node)
    #         continue
    #
    #     # if is_valid_network(G, T) and # is valid, but adding the next node increases the numbers
    #     if not is_valid_network(G, T): # still a tree, but doesn't dominate
    #         # keep this next vertex
    #         # add the neighboring nodes
    #         relevant_nodes = relevant_nodes + list(G.neighbors(node))
    #         relevant_nodes = list(set(relevant_nodes))
    #         relevant_nodes.sort(key=lambda item: node_pairwise_distance[item])
    #         continue
    #     else: # is a tree and dominates, optimize this later
    #         return T
    # return T

# Here's an example of how to run your solver.

# Usage: python3 solver.py test.in

# if __name__ == '__main__':
#     assert len(sys.argv) == 2
#     path = sys.argv[1]
#     G = read_input_file(path)
#     T = solve(G)
#     assert is_valid_network(G, T)
#     print("Average  pairwise distance: {}".format(average_pairwise_distance(T)))
#     write_output_file(T, 'outputs/test.out')

if __name__ == '__main__':
    onlyfiles = [f for f in listdir('inputs') if isfile(join('inputs', f))]
    for f in onlyfiles:
        file_name = f.split('.')[0]
        print(file_name)
        G = read_input_file('inputs/' + file_name + '.in')
        T = solve(G)
        assert is_valid_network(G, T)
        write_output_file(T, 'outputs/' + file_name + '.out')
