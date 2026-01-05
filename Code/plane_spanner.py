import math
import matplotlib.pyplot as plt
import networkx as nx
from scipy.spatial import Delaunay
import numpy as np

def generate_squashed_fibonacci_boundary(n_points):
    """
    Generates points based on a Fibonacci spiral (Golden Ratio),
    squashes them into a thin ellipse, and keeps only the outer 'boundary' points.
    This creates a complex zig-zag path along the perimeter.
    """
    P = []
    
    # Generate a dense spiral
    total_generated = n_points * 10 
    
    phi = (1 + math.sqrt(5)) / 2
    golden_angle = 2 * math.pi * (1 - 1/phi)
    
    # Squash factor (0.1 = very thin ellipse)
    aspect_ratio = 0.2
    
    # We only keep the last n_points (the outer shell)
    start_index = total_generated - n_points
    
    for i in range(start_index, total_generated):
        r = math.sqrt(i)
        theta = i * golden_angle
        
        x = r * math.cos(theta)
        y = r * math.sin(theta) * aspect_ratio
        
        P.append((x, y))
        
    return P

def euclidean_dist(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def build_graph_from_delaunay(points):
    """
    Builds a NetworkX graph from the Delaunay triangulation of the points.
    """
    points_array = np.array(points)
    tri = Delaunay(points_array)
    
    G = nx.Graph()
    # Add nodes
    for i, p in enumerate(points):
        G.add_node(i, pos=p)
        
    # Add edges from triangulation
    for simplex in tri.simplices:
        # simplex contains indices of 3 points forming a triangle
        # Add edges (u, v) for each pair in the triangle
        for i in range(3):
            u = simplex[i]
            v = simplex[(i + 1) % 3]
            if not G.has_edge(u, v):
                dist = euclidean_dist(points[u], points[v])
                G.add_edge(u, v, weight=dist)
                
    return G

def find_max_stretch_pair(G, points):
    """
    Finds the pair of nodes (u, v) with the maximum stretch.
    Stretch = (shortest path distance in G) / (Euclidean distance)
    """
    max_stretch = -1.0
    best_pair = None
    path_dist = 0.0
    direct_dist = 0.0
    
    # Compute all pairs shortest paths
    # Note: If the graph is not connected, this will only contain reachable pairs
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))
    
    n = len(points)
    for i in range(n):
        if i not in path_lengths:
            continue
            
        for j in range(i + 1, n):
            # Check reachability
            if j not in path_lengths[i]:
                continue
                
            # Graph distance
            d_G = path_lengths[i][j]
            
            # Euclidean distance
            d_E = euclidean_dist(points[i], points[j])
            
            if d_E > 1e-9: # Avoid division by zero
                stretch = d_G / d_E
                if stretch > max_stretch:
                    max_stretch = stretch
                    best_pair = (i, j)
                    path_dist = d_G
                    direct_dist = d_E
                    
    return best_pair, max_stretch, path_dist, direct_dist

def main():
    # 1. Generate points
    n = 200 # Number of boundary points
    P = generate_squashed_fibonacci_boundary(n)
    print(f"Generated {len(P)} points on Squashed Fibonacci Boundary.")

    # 2. Build Delaunay Graph
    G = build_graph_from_delaunay(P)
    print(f"Built Delaunay graph with {G.number_of_edges()} edges.")

    # 3. Find Max Stretch
    (u, v), max_stretch, d_G, d_E = find_max_stretch_pair(G, P)
    print(f"Max Stretch Pair: ({u}, {v})")
    print(f"  Euclidean Distance: {d_E:.4f}")
    print(f"  Graph Distance:     {d_G:.4f}")
    print(f"  Max Stretch:        {max_stretch:.4f}")

    # 4. Get the path for plotting
    shortest_path = nx.shortest_path(G, source=u, target=v, weight='weight')
    
    # 5. Plotting
    plt.figure(figsize=(10, 6)) # Adjusted aspect ratio
    
    # Plot Edges
    for edge in G.edges():
        p1 = P[edge[0]]
        p2 = P[edge[1]]
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'gray', alpha=0.3, zorder=1)

    # Plot Max Stretch Path
    path_x = [P[node][0] for node in shortest_path]
    path_y = [P[node][1] for node in shortest_path]
    plt.plot(path_x, path_y, 'r-', linewidth=3, label=f'Max Stretch Path (S={max_stretch:.4f})', zorder=2)

    # Plot Points
    x_coords = [p[0] for p in P]
    y_coords = [p[1] for p in P]
    plt.scatter(x_coords, y_coords, c='blue', s=20, zorder=3)
    
    # Highlight start/end of max stretch
    plt.scatter([P[u][0], P[v][0]], [P[u][1], P[v][1]], c='red', s=100, zorder=4)

    plt.title(f'Delaunay Triangulation on Squashed Fibonacci (Stretch ~ {max_stretch:.4f})')
    plt.legend()
    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.3)
    
    plt.show()

if __name__ == "__main__":
    main()
