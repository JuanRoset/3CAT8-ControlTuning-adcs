import sys
import pygame
import numpy as np
import pandas as pd
from math import sin, cos, radians, atan2
from pygame.locals import *
from pyquaternion import Quaternion


# ============================================================
# CONFIGURATION
# ============================================================

WIDTH, HEIGHT = 800, 600
FPS = 1200
INITIAL_VIEW_MODE = "orbit"

FOV = 800
VIEWER_DISTANCE = 7

BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
GRAY  = (128, 128, 128)
RED   = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE  = (0, 0, 255)
BROWN = (139, 69, 19)
ORANGE = (255, 165, 0)
PURPLE = (128, 0, 128)


# ============================================================
# INITIALIZATION
# ============================================================

pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Attitude Visualizer")
clock = pygame.time.Clock()
font = pygame.font.SysFont("consolas", 18)


# ============================================================
# GEOMETRY
# ============================================================

def normalize_dimensions():
    h, w, l = 20, 10, 34.05
    ref = max(h, w, l)
    return h/(2*ref), w/(2*ref), l/(2*ref), ref


def build_body(h, w, l):
    vertices = [
        [+h,-w,-l],[+h,+w,-l],[-h,+w,-l],[-h,-w,-l],
        [+h,-w,+l],[+h,+w,+l],[-h,+w,+l],[-h,-w,+l]
    ]
    edges = [
        (0,1),(1,2),(2,3),(3,0),
        (4,5),(5,6),(6,7),(7,4),
        (0,4),(1,5),(2,6),(3,7)
    ]
    return vertices, edges


def build_circle(center, radius, n):
    verts = []
    for i in range(n+1):
        x = center[0] + radius*cos(2*np.pi*i/n)
        y = center[1] + radius*sin(2*np.pi*i/n)
        verts.append([x,y,center[2]])
    edges = [(i,i+1) for i in range(n)]
    return verts, edges


def build_camera(ref, h, l):
    radius = 2.5 / ref
    center = [-h/2, 0, -l]
    return build_circle(center, radius, 20)


def build_ring(ref, h, l):
    ring_radius = 2*16/ref
    ring_thick  = 2*3.5/ref
    ant_radius  = 2*40/ref
    ant_phase   = -0.25
    center = [+h/2, 0, +l*2.5]

    ring_v, ring_e = build_circle(center, ring_radius, 20)

    outer_v, outer_e = build_circle(center, ring_radius+ring_thick, 20)
    offset = len(ring_v)
    outer_e = [(a+offset,b+offset) for a,b in outer_e]

    antenna_v = []
    for i in range(4):
        x = center[0] + ant_radius*cos(2*np.pi*(i/3+ant_phase))
        y = center[1] + ant_radius*sin(2*np.pi*(i/3+ant_phase))
        antenna_v.append([x,y,center[2]])

    ant_offset = len(ring_v)+len(outer_v)
    antenna_e = [(ant_offset+i,ant_offset+i+1) for i in range(3)]
    hub_index = ant_offset+4
    antenna_v.append([center[0],center[1],+l])
    antenna_e += [(ant_offset+i,hub_index) for i in range(3)]

    vertices = ring_v + outer_v + antenna_v
    edges = ring_e + outer_e + antenna_e
    return vertices, edges


AXES_VERTICES = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
AXES_EDGES = [(0,1),(0,2),(0,3)]


# ============================================================
# RENDERING UTILITIES
# ============================================================

def rotation_matrix(axis, angle_deg):
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    if axis=='x': return np.array([[1,0,0],[0,c,-s],[0,s,c]])
    if axis=='y': return np.array([[c,0,s],[0,1,0],[-s,0,c]])
    if axis=='z': return np.array([[c,-s,0],[s,c,0],[0,0,1]])
    raise ValueError("Invalid axis")


def rotate(vertices, q):
    return [q.rotate(v) for v in vertices]


def project(vertices, view_matrix):
    projected = []
    for v in vertices:
        x,y,z = np.dot(view_matrix, v)
        scale = FOV / (-z + VIEWER_DISTANCE)
        px = int(WIDTH/2 + x*scale)
        py = int(HEIGHT/2 - y*scale)
        projected.append((px,py))
    return projected


def draw_edges(vertices, edges, color):
    for a,b in edges:
        pygame.draw.line(screen, color, vertices[a], vertices[b], 2)


def draw_arrow(start, end, color):
    pygame.draw.line(screen, color, start, end, 2)

    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = np.linalg.norm([dx,dy])
    arrow_size = length * 0.1
    angle = atan2(dy, dx)
    arrow_angle = radians(35)

    p1 = (end[0] - arrow_size*cos(angle-arrow_angle),
          end[1] - arrow_size*sin(angle-arrow_angle))
    p2 = (end[0] - arrow_size*cos(angle+arrow_angle),
          end[1] - arrow_size*sin(angle+arrow_angle))

    pygame.draw.line(screen, color, end, p1, 2)
    pygame.draw.line(screen, color, end, p2, 2)


def draw_axes(vertices, colors):
    for i,(a,b) in enumerate(AXES_EDGES):
        draw_arrow(vertices[a], vertices[b], colors[i])

def apply_view(vertices, view_matrix):
    """Apply visualization rotation but keep 3D coordinates."""
    return [np.dot(view_matrix, v) for v in vertices]

def draw_edges_depth_sorted(vertices_3d, edges, color):
    """
    vertices_3d : vertices after quaternion + view rotation (NOT projected)
    edges       : list of index pairs
    """

    # Compute average Z for each edge
    edge_depth = []
    for edge in edges:
        z_avg = (vertices_3d[edge[0]][2] + vertices_3d[edge[1]][2]) / 2
        edge_depth.append((z_avg, edge))

    # Sort back-to-front (largest negative z first)
    edge_depth.sort(reverse=True)

    # Project and draw
    for _, edge in edge_depth:
        v1 = vertices_3d[edge[0]]
        v2 = vertices_3d[edge[1]]

        scale1 = FOV / (-v1[2] + VIEWER_DISTANCE)
        scale2 = FOV / (-v2[2] + VIEWER_DISTANCE)

        p1 = (int(WIDTH/2 + v1[0]*scale1),
              int(HEIGHT/2 - v1[1]*scale1))

        p2 = (int(WIDTH/2 + v2[0]*scale2),
              int(HEIGHT/2 - v2[1]*scale2))

        pygame.draw.line(screen, color, p1, p2, 2)


# ============================================================
# VIEW MATRICES
# ============================================================

center_matrix = (
    rotation_matrix('x',30) @
    rotation_matrix('x',270) @
    rotation_matrix('z',225)
)

orbit_view_pitch = 30
orbit_view_matrix = (
    rotation_matrix('x',270) @
    rotation_matrix('x',-orbit_view_pitch) @
    rotation_matrix('y',-90) @
    rotation_matrix('x',45)
)


# ============================================================
# DATA
# ============================================================

df = pd.read_csv('../data/simulation_output.csv')
dfo = pd.read_csv('../data/orbit_file.csv')

times = df['Time'].values
rates = df[['AngularRateX','AngularRateY','AngularRateZ']].values
q_e2b_all = df[['Quaternion_e2b_W','Quaternion_e2b_X','Quaternion_e2b_Y','Quaternion_e2b_Z']].values
q_e2o_all = df[['Quaternion_e2o_W','Quaternion_e2o_X','Quaternion_e2o_Y','Quaternion_e2o_Z']].values
mag_field = dfo[['MagX','MagY','MagZ']].values


# ============================================================
# BUILD OBJECTS
# ============================================================

h,w,l,ref = normalize_dimensions()
body_v, body_e = build_body(h,w,l)
camera_v, camera_e = build_camera(ref,h,l)
ring_v, ring_e = build_ring(ref,h,l)


# ============================================================
# MAIN LOOP
# ============================================================

view_mode = INITIAL_VIEW_MODE

for i in range(len(q_e2b_all)):

    # ---------------- Event Handling ----------------
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()
        if event.type == KEYDOWN and event.key == K_v:
            view_mode = "orbit" if view_mode == "inertial" else "inertial"

    # ---------------- Quaternions ----------------
    q_e2b = Quaternion(q_e2b_all[i])
    q_e2o = Quaternion(q_e2o_all[i])

    if view_mode == "inertial":
        q_body = q_e2b
        q_orbit = q_e2o
        view_matrix = center_matrix
    else:
        q_body = q_e2o.inverse * q_e2b
        q_orbit = Quaternion()
        view_matrix = orbit_view_matrix

    # ---------------- Clear Screen ----------------
    screen.fill(BLACK)

    # ---------------- Rotate Objects ----------------
    body_rot = rotate(body_v, q_body)
    camera_rot = rotate(camera_v, q_body)
    ring_rot = rotate(ring_v, q_body)

    # ---------------- Apply View Transform (3D) ----------------
    body_view = apply_view(body_rot, view_matrix)
    camera_view = apply_view(camera_rot, view_matrix)
    ring_view = apply_view(ring_rot, view_matrix)

    # ---------------- Depth Sorted Drawing ----------------
    draw_edges_depth_sorted(body_view, body_e, WHITE)
    draw_edges_depth_sorted(camera_view, camera_e, WHITE)
    draw_edges_depth_sorted(ring_view, ring_e, WHITE)

    # ---------------- Axes (no depth sorting) ----------------
    body_axes_p  = project(rotate(AXES_VERTICES, q_body), view_matrix)
    orbit_axes_p = project(rotate(AXES_VERTICES, q_orbit), view_matrix)
    fixed_axes_p = project(AXES_VERTICES, view_matrix)

    draw_axes(fixed_axes_p, [GRAY, GRAY, GRAY])
    draw_axes(orbit_axes_p, [BROWN, ORANGE, PURPLE])
    draw_axes(body_axes_p, [RED, GREEN, BLUE])

    # ---------------- Magnetic Vector ----------------
    mag_vec = mag_field[i]
    mag_unit = mag_vec / np.linalg.norm(mag_vec)

    if view_mode == "orbit":
        mag_unit = q_e2o.inverse.rotate(mag_unit)

    mag_rot = apply_view([[0,0,0], mag_unit], view_matrix)
    mag_p = project(mag_rot, np.eye(3))  # already in view frame

    draw_arrow(mag_p[0], mag_p[1], GRAY)

    # ---------------- HUD ----------------
    screen.blit(font.render(f"Time: {times[i]} s", True, WHITE), (10,10))
    screen.blit(font.render(f"Angular Rates: {rates[i]}", True, WHITE), (10,40))
    screen.blit(font.render(f"View: {view_mode}", True, WHITE), (10,70))

    pygame.display.flip()
    clock.tick(FPS)