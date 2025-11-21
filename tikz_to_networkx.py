import re
import math
import os
import random

def parse_tikz_file(file_path):
    """
    Parses a TikZ .tex file to extract a set of vertices and edges.
    
    Args:
        file_path (str): Path to the .tex file.
        
    Returns:
        dict: {'nodes': set(), 'edges': list of tuples}
    """
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} was not found.")

    with open(file_path, 'r', encoding='utf-8') as f:
        tikz_content = f.read()
    
    # Data Structures
    nodes_map = {} 
    edges_list = []
    
    # Regex for Nodes: \filldraw ... (x,y) ... {label};
    node_pattern = re.compile(r'\\filldraw.*?\((\-?[\d\.]+),(\-?[\d\.]+)\).*?node.*?\{(.*?)\};')
    
    # Regex for Coordinates in Draw: (x,y)
    coord_regex = re.compile(r'\((\-?[\d\.]+),(\-?[\d\.]+)\)')
    
    # 1. Parse Nodes
    for line in tikz_content.splitlines():
        line = line.strip()
        if line.startswith(r'\filldraw'):
            match = node_pattern.search(line)
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                raw_label = match.group(3)
                # Clean label: $v_{1}$ -> v1
                clean_label = raw_label.replace('$', '').replace('{', '').replace('}', '').replace('_', '')
                nodes_map[(x, y)] = clean_label

    def get_node_id_at(target_x, target_y, tolerance=0.1):
        closest_label = None
        min_dist = float('inf')
        for (nx, ny), label in nodes_map.items():
            dist = math.hypot(nx - target_x, ny - target_y)
            if dist < min_dist:
                min_dist = dist
                closest_label = label
        if min_dist < tolerance:
            return closest_label
        return None

    # 2. Parse Edges
    for line in tikz_content.splitlines():
        line = line.strip()
        if line.startswith(r'\draw'):
            all_coords = coord_regex.findall(line)
            if len(all_coords) >= 2:
                start_coords = all_coords[0]
                end_coords = all_coords[-1]
                
                s_x, s_y = float(start_coords[0]), float(start_coords[1])
                e_x, e_y = float(end_coords[0]), float(end_coords[1])
                
                start_node = get_node_id_at(s_x, s_y)
                end_node = get_node_id_at(e_x, e_y)
                
                if start_node and end_node:
                    edges_list.append((start_node, end_node))

    # Sorting nodes for consistent output
    sorted_nodes = sorted(list(nodes_map.values()), key=lambda x: int(x.replace('v', '')) if x.replace('v','').isdigit() else x)
    
    return {
        "nodes": sorted_nodes,
        "edges": edges_list
    }

def generate_tikz_file(nodes, edges, output_path="output_graph.tex", iterations=100, spread=6.0):
    """
    Generates a TikZ .tex file using a Force-Directed Graph Layout (Fruchterman-Reingold).
    This simulates physics to untangle the graph and minimize crossings.
    
    Args:
        nodes (list): List of node labels.
        edges (list): List of tuples (start, end).
        output_path (str): Output filename.
        iterations (int): Number of physics simulation steps.
        spread (float): The approximate width/height of the drawing area. 
                        Reduced to 10.0 to minimize whitespace.
    """
    
    # --- 1. Force-Directed Layout Algorithm ---
    
    # Constants
    width = spread
    height = spread
    area = width * height
    k = math.sqrt(area / len(nodes)) if nodes else 1.0 # Optimal distance
    
    # Initial random positions (centered around 0,0)
    # We seed random for reproducibility
    random.seed(42) 
    positions = {node: (random.uniform(-width/2, width/2), random.uniform(-height/2, height/2)) for node in nodes}
    
    # Temperature (controls how much nodes move; cools down over time)
    t = width / 10.0
    dt = t / (iterations + 1)

    for i in range(iterations):
        # Displacements for this iteration
        disp = {node: [0.0, 0.0] for node in nodes}
        
        # A. Repulsive Forces (All pairs repel)
        for v in nodes:
            for u in nodes:
                if v != u:
                    dx = positions[v][0] - positions[u][0]
                    dy = positions[v][1] - positions[u][1]
                    dist = math.hypot(dx, dy)
                    if dist < 0.01: dist = 0.01 # Avoid division by zero
                    
                    # Fr = k^2 / dist
                    repulse = (k * k) / dist
                    disp[v][0] += (dx / dist) * repulse
                    disp[v][1] += (dy / dist) * repulse

        # B. Attractive Forces (Connected nodes attract)
        for (u, v) in edges:
            if u in positions and v in positions:
                dx = positions[v][0] - positions[u][0]
                dy = positions[v][1] - positions[u][1]
                dist = math.hypot(dx, dy)
                if dist < 0.01: dist = 0.01
                
                # Fa = dist^2 / k
                attract = (dist * dist) / k
                
                disp[v][0] -= (dx / dist) * attract
                disp[v][1] -= (dy / dist) * attract
                
                disp[u][0] += (dx / dist) * attract
                disp[u][1] += (dy / dist) * attract

        # C. Apply Displacements & Temperature Cooling
        for v in nodes:
            # Gravity: Add a small pull towards the center to reduce whitespace drift
            # This keeps the graph compact
            disp[v][0] -= positions[v][0] * 0.05  # weak gravity
            disp[v][1] -= positions[v][1] * 0.05

            dx = disp[v][0]
            dy = disp[v][1]
            dist = math.hypot(dx, dy)
            if dist < 0.01: dist = 0.01
            
            # Limit movement by current temperature 't'
            move_dist = min(dist, t)
            
            positions[v] = (
                positions[v][0] + (dx / dist) * move_dist,
                positions[v][1] + (dy / dist) * move_dist
            )
        
        # Cool down
        t -= dt
        if t < 0: t = 0

    # --- 2. Post-Processing: Center the Graph ---
    if positions:
        min_x = min(p[0] for p in positions.values())
        max_x = max(p[0] for p in positions.values())
        min_y = min(p[1] for p in positions.values())
        max_y = max(p[1] for p in positions.values())
        
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        
        for node in positions:
            positions[node] = (positions[node][0] - center_x, positions[node][1] - center_y)

    # --- 3. Generate TikZ Output ---

    tex_header = r"""\documentclass[tikz,border=2pt]{standalone}
    \usepackage{tikz}
    \usetikzlibrary{arrows.meta, positioning, calc}

    \begin{document}
    \begin{tikzpicture}
    """

    tex_footer = r"""
    \end{tikzpicture}
    \end{document}
    """

    content_edges = "% --- Edges ---\n"
    content_nodes = "% --- Nodes ---\n"

    # Generate Edges
    for u, v in edges:
        if u in positions and v in positions:
            ux, uy = positions[u]
            vx, vy = positions[v]
            content_edges += f"\\draw[black, thick] ({ux:.2f},{uy:.2f}) -- ({vx:.2f},{vy:.2f});\n"

    # Generate Nodes
    colors = ['red', 'blue', 'orange', 'magenta']
    for i, node in enumerate(nodes):
        x, y = positions[node]
        color = colors[i % len(colors)]
        
        if node.startswith('v') and node[1:].isdigit():
             display_label = f"$v_{{{node[1:]}}}$"
        else:
             display_label = f"${node}$"
             
        content_nodes += f"\\filldraw[fill={color}, draw=black] ({x:.2f},{y:.2f}) circle (3pt) node[anchor=south] {{{display_label}}};\n"

    full_content = tex_header + content_edges + "\n" + content_nodes + tex_footer
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(full_content)
    
    return os.path.abspath(output_path)