from math import cos, sin, pi
from PIL import Image, ImageDraw

# Image parameters
width, height = 200, 200
cx, cy = width // 2, height // 2

radius = 4        # radius of each cell-circle
offset = 90       # distance of each center from image center
n_cells = 3       # number of circles on the ring

# global shift
shift_x = -5       # 5 right
shift_y = -5       # 5 down

# Create blank grayscale image
img = Image.new("L", (width, height), 0)
draw = ImageDraw.Draw(img)

def draw_disk(center, r, gray):
    x, y = center
    bbox = [x - r, y - r, x + r, y + r]
    draw.ellipse(bbox, fill=gray)

# Place 3 circles equally spaced on a ring, then shift
for i in range(n_cells):
    theta = 2 * pi * i / n_cells
    x = cx + int(offset * cos(theta)) + shift_x
    y = cy + int(offset * sin(theta)) + shift_y
    label = 90 + i*3   # 90, 93, 96
    draw_disk((x, y), radius, gray=label)

# Save as TIFF
img.save("astrocytes.tif", format="TIFF")
print("Saved astrocytes.tif")