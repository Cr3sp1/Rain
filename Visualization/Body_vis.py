from vpython import *


def plot_geometric_bodies(filename):

    # Read the data from the file and visualize

    with open(filename, 'r') as file:

        for line in file:

            parts = line.strip().split(',')

            body_type = parts[0]

            params = list(map(float, parts[1:]))



            if body_type == 'S':  # Sphere

                center = vector(params[0], params[1], params[2])

                radius = params[3]

                sphere(pos=center, radius=radius, color=color.blue)



            elif body_type == 'P':  # Parallelepiped

                center = vector(params[0], params[1], params[2])

                axis1 = vector(params[3], params[4], params[5])

                axis2 = vector(params[6], params[7], params[8])

                axis3 = vector(params[9], params[10], params[11])

                # VPython does not have a native parallelepiped object, so we approximate it with a box

                size = vector(abs(axis1.x) + abs(axis2.x) + abs(axis3.x),

                              abs(axis1.y) + abs(axis2.y) + abs(axis3.y),

                              abs(axis1.z) + abs(axis2.z) + abs(axis3.z))

                box(pos=center, size=size, color=color.green)



            elif body_type == 'C':  # Capsule

                point1 = vector(params[0], params[1], params[2])

                point2 = vector(params[3], params[4], params[5])

                radius = params[6]

                axis = point2 - point1

                cylinder(pos=point1, axis=axis, radius=radius, color=color.red)

                sphere(pos=point1, radius=radius, color=color.red)

                sphere(pos=point2, radius=radius, color=color.red)

plot_geometric_bodies('../data/Temp/Base.dat')


# Adjust camera position and orientation
scene.camera.pos = vector(5, 1, 0.84)
scene.camera.axis = vector(-5, 0, 0)
scene.camera.up = vector(0, 0, 1)

# Capture the scene as an image and save it to file
# img = scene.get_image()
# img.save("scene_image.png")

