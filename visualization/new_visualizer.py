from vpython import *

scene.background = color.black
scene.width = 800
scene.height = 600

earth = sphere(radius=6.371e6, texture=textures.earth, make_trail=False)

# Example satellite model
sat = box(pos=vector(7e6, 0, 0), size=vector(2e5, 1e5, 3e5), color=color.white)

# Axes (body frame)
arrow(pos=sat.pos, axis=vector(1e5, 0, 0), color=color.red)
arrow(pos=sat.pos, axis=vector(0, 1e5, 0), color=color.green)
arrow(pos=sat.pos, axis=vector(0, 0, 1e5), color=color.blue)

# Example rotation using quaternions
from scipy.spatial.transform import Rotation as R
quat = [0.9239, 0.3827, 0, 0]  # example quaternion (w,x,y,z)
r = R.from_quat([quat[1], quat[2], quat[3], quat[0]])
sat.axis = vector(*r.apply([1,0,0]))
sat.up = vector(*r.apply([0,1,0]))