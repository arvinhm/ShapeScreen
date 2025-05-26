#!/usr/bin/env python3
"""
craft_cavity_extractor.py

CRAFT-inspired cavity detection algorithm for extracting binding pockets
from PDB files within a specified grid box.

Based on the CRAFT (Cavity Recognition Assisted by Flow Transfer) algorithm
but simplified for integration with peptide screening.

Usage:
python craft_cavity_extractor.py -p PSMA.pdb -c "23.052,58.975,42.259" -s "30,30,30" \
    --min-cavity-atoms 15 -o psma_craft_cavity --surface-buffer 0.3 --min-burial 1   \
    --susceptible-radius 3.0   --susceptible-volume 2.0
"""

import os
import numpy as np
import argparse
from typing import List, Dict
from scipy.spatial import Delaunay, distance_matrix
from scipy.spatial.distance import cdist
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Try to import dependencies
try:
    from rdkit import Chem
    RDKIT = True
except ImportError:
    RDKIT = False
    print("Warning: RDKit not available")

try:
    import prody as pr
    PRODY = True
except ImportError:
    PRODY = False
    print("Warning: ProDy not available")

class CRAFTCavityExtractor:
    """CRAFT-inspired cavity detection for binding site analysis"""
    
    def __init__(self, pdb_file, center, size, 
                 probe_radius=1.4, 
                 susceptible_radius=4.0,
                 susceptible_volume=4.0,
                 min_cavity_atoms=15,
                 surface_buffer=1.0,
                 min_burial=3,
                 refine_shape=True,
                 connection_distance=3.0):
        
        self.pdb_file = pdb_file
        self.center = np.array(center, dtype=float)
        self.size = np.array(size, dtype=float)
        self.probe_radius = probe_radius
        self.susceptible_radius = susceptible_radius
        self.susceptible_volume = susceptible_volume
        self.min_cavity_atoms = min_cavity_atoms
        self.surface_buffer = surface_buffer
        self.min_burial = min_burial
        self.refine_shape = refine_shape
        self.connection_distance = connection_distance
        
        # Van der Waals radii (Å)
        self.vdw_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
            'P': 1.80, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
            'MG': 1.73, 'CA': 2.31, 'ZN': 1.39, 'FE': 2.00, 'MN': 2.00
        }
        
        # Results
        self.protein_coords = None
        self.protein_atoms = None
        self.atom_radii = None
        self.delaunay = None
        self.tetrahedra = None
        self.susceptible_tetrahedra = []
        self.cavities = []
        self.cavity_points = []
        
    def load_protein(self):
        """Load protein coordinates within grid region"""
        print(f"Loading protein from: {self.pdb_file}")
        
        if PRODY:
            # Use ProDy for better PDB parsing
            pdb = pr.parsePDB(self.pdb_file)
            protein = pdb.select("protein")
            if protein is None:
                raise ValueError("No protein atoms found!")
                
            # Filter to grid region (with buffer for Delaunay)
            buffer = 8.0  # Larger buffer for triangulation
            half = self.size/2 + buffer
            region_sel = protein.select(
                f"x > {self.center[0]-half[0]} and x < {self.center[0]+half[0]} and "
                f"y > {self.center[1]-half[1]} and y < {self.center[1]+half[1]} and "
                f"z > {self.center[2]-half[2]} and z < {self.center[2]+half[2]}"
            )
            
            if region_sel is None:
                raise ValueError("No protein atoms in grid region!")
                
            self.protein_coords = region_sel.getCoords().astype(np.float32)
            
            # Get atom information
            atom_names = region_sel.getNames()
            elements = region_sel.getElements()
            
            # Assign VdW radii
            self.atom_radii = np.array([
                self.vdw_radii.get(elem, 1.7) for elem in elements
            ], dtype=np.float32)
            
            self.protein_atoms = region_sel
            
        elif RDKIT:
            # Fallback to RDKit
            mol = Chem.MolFromPDBFile(self.pdb_file, removeHs=True)
            if mol is None:
                raise ValueError("Could not read PDB file with RDKit")
                
            conf = mol.GetConformer()
            buffer = 8.0
            half = self.size/2 + buffer
            coords = []
            radii = []
            
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                xyz = np.array([pos.x, pos.y, pos.z])
                if all(abs(xyz - self.center) <= half):
                    coords.append(xyz)
                    element = atom.GetSymbol()
                    radii.append(self.vdw_radii.get(element, 1.7))
                    
            if not coords:
                raise ValueError("No protein atoms in grid region!")
                
            self.protein_coords = np.array(coords, dtype=np.float32)
            self.atom_radii = np.array(radii, dtype=np.float32)
            
        else:
            raise ImportError("Neither ProDy nor RDKit available for PDB parsing")
            
        print(f"Loaded {len(self.protein_coords)} protein atoms in extended grid region")
        
    def build_delaunay_triangulation(self):
        """Build Delaunay triangulation of protein atoms"""
        print("Building Delaunay triangulation...")
        
        try:
            self.delaunay = Delaunay(self.protein_coords)
            self.tetrahedra = self.delaunay.simplices
            print(f"Created {len(self.tetrahedra)} tetrahedra")
            
        except Exception as e:
            print(f"Delaunay triangulation failed: {e}")
            print("This can happen with degenerate point configurations")
            raise
            
    def calculate_tetrahedron_properties(self):
        """Calculate volume and circumradius for each tetrahedron"""
        print("Calculating tetrahedron properties...")
        
        volumes = []
        circumradii = []
        circumcenters = []
        
        for tetra in tqdm(self.tetrahedra, desc="Processing tetrahedra"):
            # Get tetrahedron vertices
            vertices = self.protein_coords[tetra]
            
            # Calculate volume using determinant method
            # V = |det(v1-v0, v2-v0, v3-v0)| / 6
            v1, v2, v3, v0 = vertices
            matrix = np.array([v1-v0, v2-v0, v3-v0])
            volume = abs(np.linalg.det(matrix)) / 6.0
            volumes.append(volume)
            
            # Calculate circumcenter and circumradius
            try:
                circumcenter, circumradius = self._calculate_circumsphere(vertices)
                circumradii.append(circumradius)
                circumcenters.append(circumcenter)
            except:
                # Degenerate tetrahedron
                circumradii.append(0.0)
                circumcenters.append(np.mean(vertices, axis=0))
                
        self.volumes = np.array(volumes)
        self.circumradii = np.array(circumradii)
        self.circumcenters = np.array(circumcenters)
        
        print(f"Volume range: {self.volumes.min():.2f} - {self.volumes.max():.2f} Å³")
        print(f"Circumradius range: {self.circumradii.min():.2f} - {self.circumradii.max():.2f} Å")
        
    def _calculate_circumsphere(self, vertices):
        """Calculate circumcenter and circumradius of tetrahedron"""
        # Using the standard circumsphere calculation
        v0, v1, v2, v3 = vertices
        
        # Create matrix for circumcenter calculation
        a = np.array([
            [v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]],
            [v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]],
            [v3[0]-v0[0], v3[1]-v0[1], v3[2]-v0[2]]
        ])
        
        # Right-hand side
        b = np.array([
            np.sum((v1-v0)**2) / 2,
            np.sum((v2-v0)**2) / 2,
            np.sum((v3-v0)**2) / 2
        ])
        
        # Solve for circumcenter offset
        try:
            offset = np.linalg.solve(a, b)
            circumcenter = v0 + offset
            circumradius = np.linalg.norm(offset)
            return circumcenter, circumradius
        except np.linalg.LinAlgError:
            # Degenerate case
            center = np.mean(vertices, axis=0)
            radius = np.max([np.linalg.norm(v - center) for v in vertices])
            return center, radius
            
    def identify_susceptible_tetrahedra(self):
        """Identify tetrahedra that could be part of cavities"""
        print("Identifying susceptible tetrahedra...")
        
        # Find surface tetrahedra (have neighbors = -1)
        surface_tetrahedra = set()
        for i, neighbors in enumerate(self.delaunay.neighbors):
            if -1 in neighbors:
                surface_tetrahedra.add(i)
                
        print(f"Found {len(surface_tetrahedra)} surface tetrahedra")
        
        # Apply peeling process to remove surface tetrahedra
        peeling_radius = 10.0  # Å
        boundary_tetrahedra = set()
        
        # Initial surface tetrahedra with large circumradius
        current_surface = {i for i in surface_tetrahedra 
                          if self.circumradii[i] > peeling_radius}
        boundary_tetrahedra.update(current_surface)
        
        # Iteratively peel layers
        while current_surface:
            next_layer = set()
            for tetra_idx in current_surface:
                neighbors = self.delaunay.neighbors[tetra_idx]
                for neighbor in neighbors:
                    if (neighbor != -1 and 
                        neighbor not in boundary_tetrahedra and
                        self.circumradii[neighbor] > peeling_radius):
                        next_layer.add(neighbor)
                        
            boundary_tetrahedra.update(next_layer)
            current_surface = next_layer
            
        print(f"Identified {len(boundary_tetrahedra)} boundary tetrahedra")
        
        # Find susceptible tetrahedra (large volume and radius, not boundary)
        susceptible = []
        for i, (volume, radius) in enumerate(zip(self.volumes, self.circumradii)):
            if (volume >= self.susceptible_volume and 
                radius >= self.susceptible_radius and
                i not in boundary_tetrahedra):
                susceptible.append(i)
                
        self.susceptible_tetrahedra = susceptible
        self.boundary_tetrahedra = boundary_tetrahedra
        
        print(f"Found {len(self.susceptible_tetrahedra)} susceptible tetrahedra")
        
    def calculate_face_accessibility(self, tetra_idx):
        """Calculate which faces of a tetrahedron are accessible to probe"""
        tetra = self.tetrahedra[tetra_idx]
        accessible_faces = []
        
        # Each tetrahedron has 4 faces (triangles)
        face_indices = [
            [tetra[1], tetra[2], tetra[3]],  # Face opposite to vertex 0
            [tetra[0], tetra[2], tetra[3]],  # Face opposite to vertex 1
            [tetra[0], tetra[1], tetra[3]],  # Face opposite to vertex 2
            [tetra[0], tetra[1], tetra[2]]   # Face opposite to vertex 3
        ]
        
        for face in face_indices:
            # Get face vertices and their radii
            face_coords = self.protein_coords[face]
            face_radii = self.atom_radii[face]
            
            # Calculate maximum probe radius that can pass through this face
            max_probe_radius = self._calculate_max_probe_radius(face_coords, face_radii)
            
            if max_probe_radius >= self.probe_radius:
                accessible_faces.append(face)
                
        return accessible_faces
        
    def _calculate_max_probe_radius(self, face_coords, face_radii):
        """Calculate maximum probe radius that can pass through triangular face"""
        # Simplified calculation - project to 2D and solve
        # This is a simplified version of the CRAFT algorithm
        
        # Calculate triangle center
        center = np.mean(face_coords, axis=0)
        
        # Calculate distances from center to each vertex
        distances = [np.linalg.norm(coord - center) for coord in face_coords]
        
        # Estimate maximum probe radius (simplified)
        # Real CRAFT uses complex geometric calculations
        min_clearance = min(distances) - max(face_radii)
        
        return max(0, min_clearance)
        
    def flow_transfer_algorithm(self):
        """Connect susceptible tetrahedra using flow transfer"""
        print("Running flow transfer algorithm...")
        
        unprocessed = set(self.susceptible_tetrahedra)
        cavities = []
        
        while unprocessed:
            # Start new cavity from first unprocessed tetrahedron
            seed = next(iter(unprocessed))
            cavity_tetrahedra = set()
            queue = [seed]
            
            while queue:
                current = queue.pop(0)
                if current in unprocessed:
                    cavity_tetrahedra.add(current)
                    unprocessed.remove(current)
                    
                    # Check neighbors
                    neighbors = self.delaunay.neighbors[current]
                    for neighbor in neighbors:
                        if (neighbor != -1 and 
                            neighbor in unprocessed):
                            
                            # Check if probe can flow between tetrahedra
                            if self._can_probe_flow(current, neighbor):
                                queue.append(neighbor)
                                
            # Only keep cavities with sufficient atoms
            cavity_atoms = set()
            for tetra_idx in cavity_tetrahedra:
                cavity_atoms.update(self.tetrahedra[tetra_idx])
                
            if len(cavity_atoms) >= self.min_cavity_atoms:
                cavities.append({
                    'tetrahedra': list(cavity_tetrahedra),
                    'atoms': list(cavity_atoms),
                    'volume': sum(self.volumes[list(cavity_tetrahedra)])
                })
                
        self.cavities = cavities
        print(f"Found {len(self.cavities)} cavities")
        
        # Sort by volume (largest first)
        self.cavities.sort(key=lambda x: x['volume'], reverse=True)
        
    def _can_probe_flow(self, tetra1, tetra2):
        """Check if probe can flow between two tetrahedra"""
        # Find shared face between tetrahedra
        atoms1 = set(self.tetrahedra[tetra1])
        atoms2 = set(self.tetrahedra[tetra2])
        shared_atoms = atoms1.intersection(atoms2)
        
        if len(shared_atoms) == 3:  # Shared triangular face
            face_atoms = list(shared_atoms)
            face_coords = self.protein_coords[face_atoms]
            face_radii = self.atom_radii[face_atoms]
            
            max_probe = self._calculate_max_probe_radius(face_coords, face_radii)
            return max_probe >= self.probe_radius
            
        return False
        
    def extract_cavity_points(self, grid_spacing=0.5):
        """Extract refined cavity points with surface overlap removal"""
        print("Extracting refined cavity points...")
        
        if not self.cavities:
            print("No cavities found!")
            self.cavity_points = []
            return
            
        # Generate fine grid within original box (not extended)
        half = self.size / 2
        ranges = [
            np.arange(self.center[i] - half[i], self.center[i] + half[i], 
                     grid_spacing, dtype=np.float32)
            for i in range(3)
        ]
        
        xx, yy, zz = np.meshgrid(*ranges, indexing='xy')
        grid_points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
        
        print(f"Testing {len(grid_points):,} grid points against cavities...")
        
        cavity_points = []
        
        for cavity_idx, cavity in enumerate(self.cavities):
            print(f"Processing cavity {cavity_idx + 1}/{len(self.cavities)} "
                  f"(volume: {cavity['volume']:.1f} Å³)")
            
            cavity_atoms = cavity['atoms']
            cavity_coords = self.protein_coords[cavity_atoms]
            cavity_radii = self.atom_radii[cavity_atoms]
            
            # Refined cavity point extraction with surface overlap removal
            chunk_size = 10000
            valid_points = []
            
            for i in tqdm(range(0, len(grid_points), chunk_size), 
                         desc=f"Cavity {cavity_idx + 1}"):
                chunk = grid_points[i:i+chunk_size]
                
                # Calculate distances to ALL protein atoms (not just cavity atoms)
                all_distances = cdist(chunk, self.protein_coords)
                
                # Calculate distances to cavity atoms specifically
                cavity_distances = cdist(chunk, cavity_coords)
                
                # 1. Must be outside ALL protein atoms (with VdW radii + buffer)
                outside_all_atoms = np.all(
                    all_distances > (self.atom_radii[np.newaxis, :] + self.surface_buffer), axis=1
                )
                
                # 2. Must be reasonably close to cavity atoms (accessibility)
                min_cavity_dist = cavity_distances.min(axis=1)
                close_to_cavity = min_cavity_dist < 5.0  # Within 5Å of cavity atoms
                
                # 3. Must be in "buried" region (not surface exposed)
                # Count nearby protein atoms within 8Å
                nearby_atoms = np.sum(all_distances < 8.0, axis=1)
                buried = nearby_atoms >= self.min_burial  # Configurable burial requirement
                
                # 4. Remove points too close to protein surface
                min_protein_dist = all_distances.min(axis=1)
                not_surface = min_protein_dist > (2.0 + self.surface_buffer)  # Configurable surface distance
                
                # Combine all criteria for clean cavity shape
                valid_mask = outside_all_atoms & close_to_cavity & buried & not_surface
                valid_points.append(chunk[valid_mask])
                
            if valid_points:
                cavity_points_array = np.concatenate(valid_points)
                
                # Additional refinement: remove outliers and smooth the cavity
                if self.refine_shape and len(cavity_points_array) > 100:
                    cavity_points_array = self._refine_cavity_shape(cavity_points_array)
                
                cavity_points.append({
                    'cavity_id': cavity_idx + 1,
                    'points': cavity_points_array,
                    'volume': cavity['volume'],
                    'n_atoms': len(cavity_atoms)
                })
                
        self.cavity_points = cavity_points
        
        # Apply connectivity analysis to merged cavities
        if len(self.cavity_points) > 1:
            self.cavity_points = self._ensure_single_connected_cavity(self.cavity_points)
        
        total_points = sum(len(cp['points']) for cp in self.cavity_points)
        print(f"Extracted {total_points:,} refined cavity points from {len(self.cavity_points)} cavities")
        print("Applied surface overlap removal, shape refinement, and connectivity analysis")
    
    def _refine_cavity_shape(self, points: np.ndarray) -> np.ndarray:
        """Refine cavity shape by removing outliers and smoothing"""
        if len(points) < 100:
            return points
            
        print(f"  Refining cavity shape ({len(points):,} points)...")
        
        # 1. Remove statistical outliers using IQR method
        center = np.mean(points, axis=0)
        distances = np.linalg.norm(points - center, axis=1)
        
        Q1 = np.percentile(distances, 25)
        Q3 = np.percentile(distances, 75)
        IQR = Q3 - Q1
        
        # Keep points within 1.5 * IQR from median
        outlier_threshold = Q3 + 1.5 * IQR
        inlier_mask = distances <= outlier_threshold
        refined_points = points[inlier_mask]
        
        # 2. Density-based filtering: remove isolated points
        if len(refined_points) > 50:
            refined_points = self._remove_isolated_points(refined_points)
        
        # 3. Convex hull filtering: keep points that form a reasonable shape
        if len(refined_points) > 20:
            refined_points = self._convex_hull_filter(refined_points)
        
        removed = len(points) - len(refined_points)
        if removed > 0:
            print(f"    Removed {removed:,} outlier/surface points ({removed/len(points)*100:.1f}%)")
        
        return refined_points
    
    def _remove_isolated_points(self, points: np.ndarray, min_neighbors: int = 3, 
                               neighbor_radius: float = 2.0) -> np.ndarray:
        """Remove points that don't have enough nearby neighbors"""
        from scipy.spatial.distance import cdist
        
        # Calculate pairwise distances
        distances = cdist(points, points)
        
        # Count neighbors within radius for each point
        neighbor_counts = np.sum(distances < neighbor_radius, axis=1) - 1  # Exclude self
        
        # Keep points with sufficient neighbors
        connected_mask = neighbor_counts >= min_neighbors
        return points[connected_mask]
    
    def _convex_hull_filter(self, points: np.ndarray, expansion_factor: float = 1.2) -> np.ndarray:
        """Filter points using expanded convex hull to maintain cavity shape"""
        try:
            from scipy.spatial import ConvexHull
            
            # Calculate convex hull
            hull = ConvexHull(points)
            hull_points = points[hull.vertices]
            
            # Calculate center and expand hull slightly
            center = np.mean(hull_points, axis=0)
            expanded_hull_points = center + expansion_factor * (hull_points - center)
            
            # Create expanded hull
            expanded_hull = ConvexHull(expanded_hull_points)
            
            # Keep points inside the expanded hull
            from scipy.spatial.distance import cdist
            
            # Simple approach: keep points within reasonable distance of hull
            hull_center = np.mean(hull_points, axis=0)
            max_hull_dist = np.max(np.linalg.norm(hull_points - hull_center, axis=1))
            
            point_distances = np.linalg.norm(points - hull_center, axis=1)
            inside_mask = point_distances <= max_hull_dist * expansion_factor
            
            return points[inside_mask]
            
        except Exception as e:
            print(f"    Convex hull filtering failed: {e}, skipping...")
            return points
    
    def _ensure_single_connected_cavity(self, cavity_points_list: List[Dict]) -> List[Dict]:
        """Ensure we have only the largest connected cavity component"""
        print("Analyzing cavity connectivity...")
        
        # Merge all cavity points
        all_points = []
        total_volume = 0
        total_atoms = set()
        
        for cavity_data in cavity_points_list:
            all_points.append(cavity_data['points'])
            total_volume += cavity_data['volume']
            # Note: atoms might overlap between cavities, so we use set union
        
        if not all_points:
            return []
        
        merged_points = np.concatenate(all_points)
        print(f"  Merged {len(cavity_points_list)} cavities into {len(merged_points):,} points")
        
        # Find connected components using spatial clustering
        connected_components = self._find_connected_components(merged_points, self.connection_distance)
        
        if len(connected_components) == 1:
            print("  Single connected cavity found - no separation needed")
            # Return as single unified cavity
            return [{
                'cavity_id': 1,
                'points': merged_points,
                'volume': total_volume,
                'n_atoms': len(total_atoms) if total_atoms else sum(len(c.get('atoms', [])) for c in cavity_points_list)
            }]
        
        # Multiple components found - keep only the largest
        component_sizes = [len(comp) for comp in connected_components]
        largest_idx = np.argmax(component_sizes)
        largest_component = connected_components[largest_idx]
        
        removed_points = len(merged_points) - len(largest_component)
        print(f"  Found {len(connected_components)} disconnected components")
        print(f"  Component sizes: {sorted(component_sizes, reverse=True)}")
        print(f"  Keeping largest component: {len(largest_component):,} points")
        print(f"  Removed {removed_points:,} points from {len(connected_components)-1} smaller components")
        
        # Return only the largest connected component
        return [{
            'cavity_id': 1,
            'points': largest_component,
            'volume': total_volume * (len(largest_component) / len(merged_points)),  # Scale volume proportionally
            'n_atoms': sum(len(c.get('atoms', [])) for c in cavity_points_list)  # Keep total atom count
        }]
    
    def _find_connected_components(self, points: np.ndarray, 
                                  connection_distance: float = 3.0) -> List[np.ndarray]:
        """Find connected components in cavity points using spatial clustering"""
        if len(points) == 0:
            return []
        
        try:
            from sklearn.cluster import DBSCAN
            from scipy.spatial.distance import pdist, squareform
            
            # Use DBSCAN clustering to find connected components
            # Points within connection_distance are considered connected
            clustering = DBSCAN(eps=connection_distance, min_samples=5)
            labels = clustering.fit_predict(points)
            
            # Group points by cluster label
            components = []
            unique_labels = np.unique(labels)
            
            for label in unique_labels:
                if label == -1:  # Noise points in DBSCAN
                    continue
                component_mask = labels == label
                component_points = points[component_mask]
                if len(component_points) > 10:  # Minimum size for a valid component
                    components.append(component_points)
            
            if not components:
                # Fallback: treat all points as one component if clustering fails
                print("    DBSCAN clustering found no valid components, keeping all points")
                return [points]
            
            return components
            
        except ImportError:
            print("    scikit-learn not available, using distance-based clustering")
            return self._distance_based_clustering(points, connection_distance)
        except Exception as e:
            print(f"    Clustering failed: {e}, keeping all points as single component")
            return [points]
    
    def _distance_based_clustering(self, points: np.ndarray, 
                                  connection_distance: float = 3.0) -> List[np.ndarray]:
        """Fallback clustering method using distance-based connectivity"""
        from scipy.spatial.distance import cdist
        
        if len(points) == 0:
            return []
        
        # Build connectivity graph
        distances = cdist(points, points)
        adjacency = distances < connection_distance
        
        # Find connected components using graph traversal
        visited = np.zeros(len(points), dtype=bool)
        components = []
        
        for i in range(len(points)):
            if visited[i]:
                continue
                
            # BFS to find connected component
            component_indices = []
            queue = [i]
            
            while queue:
                current = queue.pop(0)
                if visited[current]:
                    continue
                    
                visited[current] = True
                component_indices.append(current)
                
                # Add unvisited neighbors to queue
                neighbors = np.where(adjacency[current] & ~visited)[0]
                queue.extend(neighbors)
            
            if len(component_indices) > 10:  # Minimum component size
                components.append(points[component_indices])
        
        if not components:
            # Fallback: return all points as single component
            return [points]
            
        return components
        
    def save_results(self, output_prefix="craft_cavity"):
        """Save cavity results to files"""
        print("Saving results...")
        
        if not self.cavity_points:
            print("No cavity points to save!")
            return
        
        # Save individual cavity points
        all_vis_points = []
        for i, cavity_data in enumerate(self.cavity_points):
            points = cavity_data['points']
            if len(points) > 0:
                # Subsample for visualization
                max_points = 5000
                if len(points) > max_points:
                    indices = np.random.choice(len(points), max_points, replace=False)
                    vis_points = points[indices]
                else:
                    vis_points = points
                    
                filename = f"{output_prefix}_cavity_{i+1:02d}.pdb"
                with open(filename, 'w') as f:
                    for j, (x, y, z) in enumerate(vis_points, 1):
                        f.write(f"HETATM{j:5d}  C   CAV A   1    "
                               f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 0.00           C\n")
                        
                print(f"Saved cavity {i+1} to {filename} "
                      f"({len(vis_points):,}/{len(points):,} points)")
                
                # Collect subsampled points for unified cavity (same as individual files)
                all_vis_points.append(vis_points)
        
        # Create unified cavity file by merging all individual cavity subsampled points
        if all_vis_points:
            unified_points = np.concatenate(all_vis_points)
            
            # Save unified cavity using same subsampling as individual files
            unified_filename = f"{output_prefix}_unified_cavity.pdb"
            with open(unified_filename, 'w') as f:
                for j, (x, y, z) in enumerate(unified_points, 1):
                    f.write(f"HETATM{j:5d}  C   CAV A   1    "
                           f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 0.00           C\n")
            print(f"Saved unified cavity to {unified_filename} ({len(unified_points):,} points)")
            print(f"  (This should be similar size to sum of individual cavity files)")
                      
        # Save protein atoms for reference
        if hasattr(self, 'protein_atoms') and PRODY:
            pr.writePDB(f"{output_prefix}_protein.pdb", self.protein_atoms)
            print(f"Saved protein to {output_prefix}_protein.pdb")
            
        # Save cavity analysis
        analysis_file = f"{output_prefix}_analysis.txt"
        with open(analysis_file, 'w') as f:
            f.write("CRAFT Refined Cavity Analysis\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Protein: {self.pdb_file}\n")
            f.write(f"Grid center: {self.center}\n")
            f.write(f"Grid size: {self.size}\n")
            f.write(f"Probe radius: {self.probe_radius} Å\n")
            f.write(f"Surface buffer: {self.surface_buffer} Å\n")
            f.write(f"Min burial atoms: {self.min_burial}\n")
            f.write(f"Shape refinement: {self.refine_shape}\n")
            f.write(f"Connection distance: {self.connection_distance} Å\n")
            f.write(f"Total cavities found: {len(self.cavity_points)}\n")
            if all_vis_points:
                f.write(f"Unified cavity points: {len(unified_points):,}\n\n")
            
            for i, cavity_data in enumerate(self.cavity_points):
                f.write(f"Cavity {i+1}:\n")
                f.write(f"  Volume: {cavity_data['volume']:.1f} Å³\n")
                f.write(f"  Atoms: {cavity_data['n_atoms']}\n")
                f.write(f"  Points: {len(cavity_data['points']):,}\n\n")
                
        print(f"Saved analysis to {analysis_file}")
        
        # PyMOL visualization commands
        print(f"\n🔬 PyMOL Visualization Commands:")
        print(f"   load {output_prefix}_protein.pdb")
        print(f"   load {output_prefix}_unified_cavity.pdb")
        print(f"   show cartoon, {output_prefix}_protein")
        print(f"   show spheres, {output_prefix}_unified_cavity")
        print(f"   color green, {output_prefix}_unified_cavity")
        print(f"   set sphere_scale, 0.3")
        print(f"\n   # Individual cavities:")
        for i in range(len(self.cavity_points)):
            print(f"   load {output_prefix}_cavity_{i+1:02d}.pdb")
            print(f"   color {'blue' if i == 1 else 'red' if i == 2 else 'yellow'}, {output_prefix}_cavity_{i+1:02d}")
        
    def run_analysis(self, output_prefix="craft_cavity"):
        """Run complete CRAFT cavity analysis"""
        print("=" * 60)
        print("CRAFT CAVITY EXTRACTION")
        print("=" * 60)
        
        try:
            self.load_protein()
            self.build_delaunay_triangulation()
            self.calculate_tetrahedron_properties()
            self.identify_susceptible_tetrahedra()
            self.flow_transfer_algorithm()
            self.extract_cavity_points()
            self.save_results(output_prefix)
            
            print("\n" + "=" * 60)
            print("ANALYSIS COMPLETE")
            print("=" * 60)
            
            if self.cavity_points:
                print("✅ Successfully extracted cavities using CRAFT algorithm")
                print(f"Found {len(self.cavity_points)} cavities")
                largest_cavity = max(self.cavity_points, key=lambda x: x['volume'])
                print(f"Largest cavity: {largest_cavity['volume']:.1f} Å³ "
                      f"with {len(largest_cavity['points']):,} points")
                total_points = sum(len(cp['points']) for cp in self.cavity_points)
                print(f"Total unified cavity points: {total_points:,}")
            else:
                print("⚠️  No cavities found - try adjusting parameters")
                
        except Exception as e:
            print(f"❌ Error during analysis: {e}")
            import traceback
            traceback.print_exc()

def main():
    parser = argparse.ArgumentParser(description="CRAFT-inspired cavity extraction")
    parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
    parser.add_argument("-c", "--center", required=True, help="Grid center 'x,y,z'")
    parser.add_argument("-s", "--size", required=True, help="Grid size 'sx,sy,sz'")
    parser.add_argument("--probe-radius", type=float, default=1.4,
                       help="Probe radius in Angstroms (default: 1.4)")
    parser.add_argument("--susceptible-radius", type=float, default=4.0,
                       help="Susceptible tetrahedron radius (default: 4.0)")
    parser.add_argument("--susceptible-volume", type=float, default=4.0,
                       help="Susceptible tetrahedron volume (default: 4.0)")
    parser.add_argument("--min-cavity-atoms", type=int, default=15,
                       help="Minimum atoms for cavity (default: 15)")
    parser.add_argument("--surface-buffer", type=float, default=1.0,
                       help="Buffer distance from protein surface (default: 1.0)")
    parser.add_argument("--min-burial", type=int, default=3,
                       help="Minimum nearby atoms for burial (default: 3)")
    parser.add_argument("--refine-shape", action="store_true",
                       help="Apply shape refinement to remove outliers")
    parser.add_argument("--connection-distance", type=float, default=3.0,
                       help="Maximum distance for cavity connectivity (default: 3.0)")
    parser.add_argument("-o", "--output", default="craft_cavity",
                       help="Output prefix (default: craft_cavity)")
    
    args = parser.parse_args()
    
    # Parse center and size
    center = [float(x) for x in args.center.split(',')]
    size = [float(x) for x in args.size.split(',')]
    
    if len(center) != 3 or len(size) != 3:
        raise ValueError("Center and size must have 3 values each")
    
    # Run CRAFT analysis
    extractor = CRAFTCavityExtractor(
        args.protein, center, size,
        probe_radius=args.probe_radius,
        susceptible_radius=args.susceptible_radius,
        susceptible_volume=args.susceptible_volume,
        min_cavity_atoms=args.min_cavity_atoms,
        surface_buffer=args.surface_buffer,
        min_burial=args.min_burial,
        refine_shape=args.refine_shape,
        connection_distance=args.connection_distance
    )
    
    extractor.run_analysis(args.output)

if __name__ == "__main__":
    main() 