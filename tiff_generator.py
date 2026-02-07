from math import cos, sin, pi
from PIL import Image, ImageDraw

def make_centered_circles_tiff(
    filename: str,
    width: int = 200,
    height: int = 200,
    radius: int = 15,
    offset: int = 20,
    n_cells: int = 2
):
    """
    Create a TIFF with n_cells circular regions around the image center.
    Each circle gets a unique gray label (1, 2, ..., n_cells).
    """
    cx, cy = width // 2, height // 2

    img = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(img)

    def draw_disk(center, r, gray):
        x, y = center
        bbox = [x - r, y - r, x + r, y + r]
        draw.ellipse(bbox, fill=gray)

    # Place n_cells equally spaced on a circle of radius = offset
    for i in range(n_cells):
        theta = 2 * pi * i / n_cells
        x = cx + int(offset * cos(theta))
        y = cy + int(offset * sin(theta))
        label = i + 1  # labels 1..n_cells
        draw_disk((x, y), radius, gray=label)

    img.save(filename, format="TIFF")
    print(f"Saved {filename}")

# Example usage:
# two cells, slightly spaced horizontally around center
make_centered_circles_tiff("two_cells_centered.tif", radius=15, offset=20, n_cells=2)