import pygame
import sys
from pygame.locals import *
from math import sin, cos, radians
from pyquaternion import Quaternion
import pandas as pd
import numpy as np

pygame.init()

# Constants
WIDTH, HEIGHT = 800, 600
FPS = 1200

# Sat dimensions
body_height = 20
body_width = 10
body_length = 34.05
reference = max([body_width, body_height, body_length])
body_height /= reference * 2
body_width /= reference * 2
body_length /= reference * 2

# Colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
YELLOW = (255, 255, 0)
VIOLET = (128, 0, 128)  # A common shade of violet
ORANGE = (255, 165, 0)  # A bright orange
PURPLE = (128, 0, 128)  # Same as violet for consistency
BROWN = (139, 69, 19)  # A medium brown
GRAY = (128, 128, 128)

# Initialize Pygame window
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Attitude Visualizer")
clock = pygame.time.Clock()
available_fonts = pygame.font.get_fonts()

font = pygame.font.SysFont("consolas", 18)  # Choose a font and size

# Satellite body vertices
body_vertices = [
    [+ body_height, - body_width, - body_length],
    [+ body_height, + body_width, - body_length],
    [- body_height, + body_width, - body_length],
    [- body_height, - body_width, - body_length],
    [+ body_height, - body_width, + body_length],
    [+ body_height, + body_width, + body_length],
    [- body_height, + body_width, + body_length],
    [- body_height, - body_width, + body_length]
]

# Satellite body edges (connecting vertices)
body_edges = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
]

camera_vertices = []
camera_radius = 2.5
camera_radius /= reference
camera_center = [- body_height / 2, 0, - body_length]
camera_poly_number = 20
for num in range(camera_poly_number + 1):
    x = camera_center[0] + camera_radius * cos(2 * np.pi * num / camera_poly_number)
    y = camera_center[1] + camera_radius * sin(2 * np.pi * num / camera_poly_number)
    camera_vertices.append([x, y, camera_center[2]])
camera_edges = [(i, i + 1) for i in range(camera_poly_number)]

ring_vertices = []
ring_radius = 2 * 16
ring_thickness = 2 * 3.5
antenna_radius = 2 * 40
ring_radius /= reference
ring_thickness /= reference
antenna_radius /= reference
antenna_phase = - 0.25
ring_center = [+ body_height / 2, 0, + body_length * 2.5]
ring_poly_number = 20
antenna_poly_number = 3
for num in range(ring_poly_number + 1):
    x = ring_center[0] + ring_radius * cos(2 * np.pi * num / ring_poly_number)
    y = ring_center[1] + ring_radius * sin(2 * np.pi * num / ring_poly_number)
    ring_vertices.append([x, y, ring_center[2]])
for num in range(ring_poly_number + 1):
    x = ring_center[0] + (ring_radius + ring_thickness) * cos(2 * np.pi * num / ring_poly_number)
    y = ring_center[1] + (ring_radius + ring_thickness) * sin(2 * np.pi * num / ring_poly_number)
    ring_vertices.append([x, y, ring_center[2]])
for num in range(antenna_poly_number + 1):
    x = ring_center[0] + antenna_radius * cos(2 * np.pi * (num / antenna_poly_number + antenna_phase))
    y = ring_center[1] + antenna_radius * sin(2 * np.pi * (num / antenna_poly_number + antenna_phase))
    ring_vertices.append([x, y, ring_center[2]])
ring_vertices.append([ring_center[0], ring_center[1], + body_length])
ring_edges = [(i, i + 1) for i in range(ring_poly_number)]
ring_edges += [(i + 1 + ring_poly_number, i + 2 + ring_poly_number) for i in range(ring_poly_number)]
ring_edges += [(i + 2 + 2 * ring_poly_number, i + 3 + 2 * ring_poly_number) for i in range(antenna_poly_number)]
ring_edges += [(i + 2 + 2 * ring_poly_number, len(ring_vertices) - 1) for i in range(antenna_poly_number)]

# Axes vertices
axes_vertices = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]

# Axes edges
axes_edges = [
    (0, 1),  # x-axis
    (0, 2),  # y-axis
    (0, 3)   # z-axis
]

# Arrow vertices
arrow_vertices = [
    [0, 0, 0],
    [1, 0, 0]
]

# Arrow edges
arrow_edges = [
    (0, 1)  # x-axis
]

# Camera parameters
fov = 800  # Field of view
viewer_distance = 7

def project_vertices(vertices, viewer_distance, rotation_matrix):
    projected_vertices = []
    for vertex in vertices:
        # Apply rotation
        rotated_vertex = np.dot(rotation_matrix, vertex)
        x, y, z = rotated_vertex
        scale = fov / (-z + viewer_distance)
        screen_x = int(WIDTH / 2 + x * scale)
        screen_y = int(HEIGHT / 2 - y * scale)
        projected_vertices.append((screen_x, screen_y))
    return projected_vertices

def rotation_matrix(axis, angle_deg):
    angle_rad = np.radians(angle_deg)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    if axis == 'x':
        return np.array([[1, 0, 0],
                         [0, c, -s],
                         [0, s, c]])
    elif axis == 'y':
        return np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])
    elif axis == 'z':
        return np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")
    
def draw_arrow(arrow_vertices, edge, color):
    arrow_scale = 0.1  # Adjust arrow size as needed
    arrow_angle = 35
    start_vertex = arrow_vertices[0]
    end_vertex = arrow_vertices[1]
    delta_size = [start_vertex[0] - end_vertex[0], start_vertex[1] - end_vertex[1]]
    arrow_size = np.linalg.norm(delta_size) * arrow_scale
    # arrow_size = np.linalg.norm(end_vertex - start_vertex)
    pygame.draw.line(screen, color, start_vertex, end_vertex, 2)
    
    
    dx = end_vertex[0] - start_vertex[0]
    dy = end_vertex[1] - start_vertex[1]
    angle = np.arctan2(dy, dx)
    # Calculate arrowhead points
    arrow_point1 = (end_vertex[0] + arrow_size * cos(angle - radians(180 - arrow_angle)), 
                    end_vertex[1] + arrow_size * sin(angle - radians(180 - arrow_angle)))
    arrow_point2 = (end_vertex[0] + arrow_size * cos(angle + radians(180 - arrow_angle)), 
                    end_vertex[1] + arrow_size * sin(angle + radians(180 - arrow_angle)))
    # Draw arrowhead
    pygame.draw.line(screen, color, end_vertex, arrow_point1, 2)
    pygame.draw.line(screen, color, end_vertex, arrow_point2, 2)

def draw_axes(axes_vertices, axes_edges, X_color, Y_color, Z_color):
    arrow_scale = 0.1  # Adjust arrow size as needed
    arrow_angle = 35
    for edge in axes_edges:
        start_vertex = axes_vertices[edge[0]]
        end_vertex = axes_vertices[edge[1]]
        delta_size = [start_vertex[0] - end_vertex[0], start_vertex[1] - end_vertex[1]]
        arrow_size = np.linalg.norm(delta_size) * arrow_scale
        # arrow_size = np.linalg.norm(end_vertex - start_vertex)
        if edge[1] == 1:  # x-axis
            color = X_color
        elif edge[1] == 2:  # y-axis
            color = Y_color
        else:  # z-axis
            color = Z_color
        pygame.draw.line(screen, color, start_vertex, end_vertex, 2)
        # Draw arrowhead
        if edge[1] != 0:  # Omit arrowhead for the origin
            dx = end_vertex[0] - start_vertex[0]
            dy = end_vertex[1] - start_vertex[1]
            angle = np.arctan2(dy, dx)
            # Calculate arrowhead points
            arrow_point1 = (end_vertex[0] + arrow_size * cos(angle - radians(180 - arrow_angle)), 
                            end_vertex[1] + arrow_size * sin(angle - radians(180 - arrow_angle)))
            arrow_point2 = (end_vertex[0] + arrow_size * cos(angle + radians(180 - arrow_angle)), 
                            end_vertex[1] + arrow_size * sin(angle + radians(180 - arrow_angle)))
            # Draw arrowhead
            pygame.draw.line(screen, color, end_vertex, arrow_point1, 2)
            pygame.draw.line(screen, color, end_vertex, arrow_point2, 2)

# Create rotation matrix for immproved visualization
center_matrix_a = rotation_matrix('x', 270)
center_matrix_b = rotation_matrix('z', 225)
center_matrix_c = rotation_matrix('x', 30)
#center_matrix_c = rotation_matrix
center_matrix = np.dot(center_matrix_c, np.dot(center_matrix_a, center_matrix_b))

# Read quaternion data from the CSV file
file_path = '../../data/simulation_output.csv'
orbit_path = '../../data/orbit_file.csv'
df = pd.read_csv(file_path)
times = df['Time'].values
angular_rates = df[['AngularRateX', 'AngularRateY', 'AngularRateZ']].values
quaternions_e2b = df[['Quaternion_e2b_W', 'Quaternion_e2b_X', 'Quaternion_e2b_Y', 'Quaternion_e2b_Z']].values
quaternions_e2o = df[['Quaternion_e2o_W', 'Quaternion_e2o_X', 'Quaternion_e2o_Y', 'Quaternion_e2o_Z']].values
dfo = pd.read_csv(orbit_path)
magnetic_field = dfo[['MagX', 'MagY', 'MagZ']].values

# Main game loop
for i in range(len(quaternions_e2b)):
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()

    # Get quaternion for the current frame
    current_quaternion = Quaternion(quaternions_e2b[i])
    orbit_quaternion = Quaternion(quaternions_e2o[i])
    
    # Get magnetic field for the current frame
    current_magnetic_field = magnetic_field[i].tolist()
    current_magnetic_unit = current_magnetic_field / np.linalg.norm(current_magnetic_field)

    # Rotate the cube using quaternion rotation
    rotated_body_vertices = [current_quaternion.rotate(v) for v in body_vertices]
    rotated_camera_vertices = [current_quaternion.rotate(v) for v in camera_vertices]
    rotated_ring_vertices = [current_quaternion.rotate(v) for v in ring_vertices]
    rotated_body_axes_vertices = [current_quaternion.rotate(v) for v in axes_vertices]
    rotated_orbit_axes_vertices = [orbit_quaternion.rotate(v) for v in axes_vertices]

    # Project Satellite vertices to 2D screen coordinates
    projected_body_vertices = project_vertices(rotated_body_vertices, viewer_distance, center_matrix)
    projected_camera_vertices = project_vertices(rotated_camera_vertices, viewer_distance, center_matrix)
    projected_ring_vertices = project_vertices(rotated_ring_vertices, viewer_distance, center_matrix)

    # Project axes vertices to 2D screen coordinates
    projected_body_axes_vertices = project_vertices(rotated_body_axes_vertices, viewer_distance, center_matrix)
    projected_orbit_axes_vertices = project_vertices(rotated_orbit_axes_vertices, viewer_distance, center_matrix)
    projected_axes_vertices = project_vertices(axes_vertices, viewer_distance, center_matrix)

    mag_vertices = [[0, 0, 0], current_magnetic_unit]
    projected_mag_vertices = project_vertices(mag_vertices, viewer_distance, center_matrix)

    # Fill screen black
    screen.fill(BLACK)

    # Draw the body, orbit and fixed axes
    draw_axes(projected_axes_vertices, axes_edges, GRAY, GRAY, GRAY)
    draw_axes(projected_orbit_axes_vertices, axes_edges, BROWN, ORANGE, PURPLE)
    draw_axes(projected_body_axes_vertices, axes_edges, RED, GREEN, BLUE)

    draw_arrow(projected_mag_vertices, arrow_edges, GRAY)
    

    # Draw the cube
    for edge in body_edges:
        start_vertex = projected_body_vertices[edge[0]]
        end_vertex = projected_body_vertices[edge[1]]
        pygame.draw.line(screen, WHITE, start_vertex, end_vertex, 2)

    for edge in camera_edges:
        start_vertex = projected_camera_vertices[edge[0]]
        end_vertex = projected_camera_vertices[edge[1]]
        pygame.draw.line(screen, WHITE, start_vertex, end_vertex, 2)

    for edge in ring_edges:
        start_vertex = projected_ring_vertices[edge[0]]
        end_vertex = projected_ring_vertices[edge[1]]
        pygame.draw.line(screen, WHITE, start_vertex, end_vertex, 2)

    # Display time and angular rates on the screen
    time_text = font.render(f"Time: {times[i]} s", True, WHITE)
    screen.blit(time_text, (10, 10))
    angular_rate_text = font.render(f"Angular Rates (X,Y,Z): {angular_rates[i]}", True, WHITE)
    screen.blit(angular_rate_text, (10, 50))

    pygame.display.flip()
    clock.tick(FPS)