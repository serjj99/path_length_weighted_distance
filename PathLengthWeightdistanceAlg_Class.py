'''
Algorithm: [Path-length-weights distance]
Authors: [S. Sanjuan, R. Arnau, E. A. Sánchez Pérez, J.M. Calabuig, L.M. García Raffi]
Affiliation: [UPV]
Date: [27/06/2024]
Article: [Título del Artículo, Nombre de la Revista]
Version: 0.0.1

Dependencies:
    - numpy (versión 1.23.5)
    
Usage:
    Example:
    ```python
    import numpy as np
    result = my_algorithm(input_data)
    print(result)
    ```
'''

import numpy as np


'''
This function filters the paths by means of equation 3.3 of the article (cite)
'''
def P_star(P_ab): # P_ab are the paths from a to b
    # List for storing rows that are not dominated
    P_star_list = []
    
    # Traverse each row P in P_ab
    for P in P_ab:
        s_P, l_P = P
        is_dominated = False
        # Verify that there is no P' in P_ab that dominates P
        for P_prime in P_ab:
            if not np.array_equal(P, P_prime):
                s_P_prime, l_P_prime = P_prime
                if s_P_prime <= s_P and l_P_prime >= l_P: # Filter
                    is_dominated = True
                    break
        # If P is not dominated by any P', add it to P_star_list
        if not is_dominated:
            P_star_list.append(P)
    
    # Convert the list of non-dominated paths back to a NumPy array and delete the duplicates
    P_star_array = np.array(P_star_list)
    P_star_array = np.unique(P_star_array, axis=0)
    
    return P_star_array


'''
This function filters the paths by means of equation 6.3 of the article (cite)
'''
def P_star_1(P_ab): # P_ab are the paths from a to b
    # List for storing rows that are not dominated
    P_star_1_list = []
    
    # Traverse each row P in P_ab
    for P in P_ab:
        s_P, l_P = P
        is_dominated = False
        # Verify that there is no P' in P_ab that dominates P
        for P_prime in P_ab:
            if not np.array_equal(P, P_prime):
                s_P_prime, l_P_prime = P_prime
                if s_P_prime <= s_P and (l_P_prime*s_P_prime) >= (l_P*s_P): # Filter
                    is_dominated = True
                    break
        # If P is not dominated by any P', add it to P_star_1_list
        if not is_dominated:
            P_star_1_list.append(P)
    
    # Convert the list of non-dominated paths back to a NumPy array and delete the duplicates
    P_star_1_array = np.array(P_star_1_list)
    P_star_1_array = np.unique(P_star_1_array, axis=0)
    
    return P_star_1_array

'''
This function filters the paths by means of equation 6.5 of the article (cite)
'''
def P_star_2(P_ab): # P_ab are the paths from a to b
    # List for storing rows that are not dominated
    P_star_2_list = []
    
    # Traverse each row P in P_ab
    for P in P_ab:
        s_P, l_P = P
        is_dominated = False
        # Verify that there is no P' in P_ab that dominates P
        for P_prime in P_ab:
            if not np.array_equal(P, P_prime):
                s_P_prime, l_P_prime = P_prime
                if s_P_prime >= s_P and (l_P_prime*s_P_prime) >= (l_P*s_P): # Filter
                    is_dominated = True
                    break
        # If P is not dominated by any P', add it to P_star_2_list
        if not is_dominated:
            P_star_2_list.append(P)
    
    # Convert the list of non-dominated paths back to a NumPy array and delete the duplicates
    P_star_2_array = np.array(P_star_2_list)
    P_star_2_array = np.unique(P_star_2_array, axis=0)
    
    return P_star_2_array


class Graph:
    
    '''
    Function that creates a graph are connections along with the weights
    associated with the length of the paths (for the calculation of the
                                             path lenght weight distance).
    '''
    def __init__(self, vertex, weight=1):
        self.V = vertex
        
        self.edges = np.full((vertex, vertex), np.inf)
        np.fill_diagonal(self.edges, 0)
        
        # weights are the weights associated with the length of the paths
        if weight > 0 and isinstance(weight, int):
            self.weights = [1]
            for i in range(1,vertex+1):
                self.weights.append(1/i**weight)
        else:
            self.weights = weight


    '''
    This function is used to add a weighted edge to the graph.
    '''
    def add_edge(self, u, v, w): # Edge from u to v with weight w
        self.edges[u, v] = w 


    '''
    In case we have the proximity matrix of the graph, with this function
    we can associate this matrix to the new graph.
    '''
    def new_edge_matrix(self, M):
        self.edges = M
        
    
    '''
    Auxiliary function that performs a depth-first search (DFS) to detect cycles in the graph.
    '''
    def isCyclicUtil(self, v, visited, recStack):
        visited[v] = True
        recStack[v] = True
 
        neighbour_v = list(*np.where(self.edges[v] != np.inf))
        neighbour_v.remove(v)
        for neighbour in neighbour_v:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, recStack) == True:
                    return True
            elif recStack[neighbour] == True:
                return True
 
        recStack[v] = False
        return False
 

    '''
    This function uses isCyclicUtil to check if there is a cycle in the whole graph.
    It is the main function that initiates the cycle detection.
    '''
    def isCyclic(self):
        visited = [False] * (self.V + 1)
        recStack = [False] * (self.V + 1)
        for node in range(self.V):
            if visited[node] == False:
                if self.isCyclicUtil(node, visited, recStack) == True:
                    return True
        return False
    
    
    
    def dijkstra(self, vn):
        phi = self.edges
        
        n = len(phi)
        visitados = [False] * n
        distancia = [float('inf')] * n
        distancia[vn] = 0

        for _ in range(n):
            u = -1
            for i in range(n):
                if not visitados[i] and (u == -1 or distancia[i] < distancia[u]):
                    u = i
            visitados[u] = True

            for v in range(n):
                if phi[u][v] != 0:
                    nueva_distancia = distancia[u] + phi[u][v]
                    if nueva_distancia < distancia[v]:
                        distancia[v] = nueva_distancia

        return distancia
    

    
    def path_lenght_weight(self, vn, filt=0): #G = (V, phi) 
        if self.isCyclic():
            print('The graph has cycles')
            return None, None
        
        else:
            M = self.V
            V = {int(i) for i in range(M)}
            phi = self.edges
            W = self.weights
            
            # Select the Pareto frontier filter
            if filt==1:
                pareto_frontier = P_star_1
            elif filt==2:
                pareto_frontier = P_star_2
            else:
                pareto_frontier = P_star
            
            # Initialization
            D = [np.array([[np.inf, 0]]) for _ in range(M)]
            D[vn] = np.array([[0,0]])

            Q = {vn}
            m = 1 
            while (m < M) and (len(Q) != 0): # Each time it loops it adds 1 to the path length (l).
                
                for vi in Q:
                    Ai = {vj for vj in V-{vi} if phi[vj,vi] != np.inf} # Set of all the nodes of the graph that are connected to the nodes of Q (except the same node).
                    
                    for vj in Ai:
                        frontier = False # This variable indicates whether it is necessary to filter the paths (because new ones are added) or not.
                        
                        for s_i,l_i in D[vi]:
                            l_i = int(l_i)
                            s_j = phi[vj,vi] + s_i # Calculation of s
                            
                            if not any(np.array_equal([s_j, l_i+1], fila) for fila in D[vj]): # If the path is new then it is saved, and frontier is activated to subsequently perform the Pareto front.
                                D[vj] = np.vstack([D[vj], [s_j, l_i+1]])
                                frontier = True
                        
                        if frontier:
                            D[vj] = pareto_frontier(D[vj]) # Filters the paths using the frontier pareto set in the start variable 'filt'.
                
                Q = {vi for vi in V-Q if D[vi][0][0] != np.inf} # Set of all NOT VISITED nodes with distance < inf for path length.
                m += 1
            
            # Distance computation
            d = []
            l = []
            for v_i in V:
                min_d, min_l = min((W[int(l_i)] * s_i, int(l_i)) for s_i, l_i in D[v_i])
                d.append(min_d)
                l.append(min_l)
            
        return np.array(d), np.array(l) # The distances together with the lengths of the paths in NumPy array format are returned as output.
    

'''
This function calculates the distance matrices of the graph.
First it calculates the weighted-path distance matrix (using the djikastra algorithm) and,
 if it has no cycles, it calculates the path-lenght-weighted distance matrix (using our algorithm).
'''    
def calc_distance_matrices(graph, filt=0):
    D, t = graph.path_lenght_weight(0)
    
    if D[0] == None: #The graph has cycles: the path-lenght-weighted distance matrix can not be calculated
        print('The graph has cycles')
        distancias_dijkstra = []
        for vn in range(graph.V):
            D_dijkstra = graph.dijkstra(vn)
            distancias_dijkstra.append(D_dijkstra)
        
        distancias_dijkstra = np.vstack(distancias_dijkstra)
        
        return None, distancias_dijkstra
    
    
    else:
        distancias_dijkstra = []
        distancias = []
        for vn in range(graph.V):
            D, t = graph.path_lenght_weight(vn, filt=filt)
            distancias.append(D)
            
            D_dijkstra = graph.dijkstra(vn)
            distancias_dijkstra.append(D_dijkstra)
        
        distancias = np.vstack(distancias)
        distancias = distancias.T
        
        distancias_dijkstra = np.vstack(distancias_dijkstra)

        return distancias, distancias_dijkstra