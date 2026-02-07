from PIL import Image, ImageDraw

# Image parameters
width, height = 200, 200
cx, cy = width // 2, height // 2

radius = 4          # radius of each cell-circle
offset = 35    # half the distance between the two centers

# Centers of the two circles (horizontally spaced around the image center)
c1 = (cx - offset, cy - offset)

# Create blank grayscale image
img = Image.new("L", (width, height), 0)
draw = ImageDraw.Draw(img)

def draw_disk(center, r, gray):
    x, y = center
    bbox = [x - r, y - r, x + r, y + r]
    draw.ellipse(bbox, fill=gray)

# First circle (cell 1, gray value 1)
draw_disk(c1, radius, gray=90)

# Save as TIFF
img.save("pericyte_cells2.tif", format="TIFF")
print("Saved pericyte_cells2.tif")
