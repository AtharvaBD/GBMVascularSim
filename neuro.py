from PIL import Image, ImageDraw
import math


def make_neuron_square_with_corners(
    filename: str = "neurons_square_perimeter_plus_inner_corners.tif",
    width: int = 300,
    height: int = 300,
    side: int = 250,          # square side length
    n_perimeter: int = 20,    # neurons along perimeter
    neuron_radius: int = 4,
    corner_offset: int = 30,  # how far inside the square corners
):
    """
    Create a TIFF with:
      - n_perimeter neurons along the perimeter of a centered square
      - 4 extra neurons just inside the square corners.
    Each neuron gets a unique gray label for Morpheus (1..N).
    """

    cx, cy = width // 2, height // 2
    half_side = side // 2

    # square bounds
    x_min = cx - half_side
    x_max = cx + half_side
    y_min = cy - half_side
    y_max = cy + half_side

    img = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(img)

    def draw_disk(center, r, gray):
        x, y = center
        bbox = [x - r, y - r, x + r, y + r]
        draw.ellipse(bbox, fill=gray)

    # place neurons along the perimeter
    label = 1
    for i in range(n_perimeter):
        t = i / n_perimeter  # goes from 0 to <1

        if t < 0.25:
            # top edge: left -> right
            alpha = t / 0.25
            x = x_min + int(alpha * (x_max - x_min))
            y = y_min
        elif t < 0.5:
            # right edge: top -> bottom
            alpha = (t - 0.25) / 0.25
            x = x_max
            y = y_min + int(alpha * (y_max - y_min))
        elif t < 0.75:
            # bottom edge: right -> left
            alpha = (t - 0.5) / 0.25
            x = x_max - int(alpha * (x_max - x_min))
            y = y_max
        else:
            # left edge: bottom -> top
            alpha = (t - 0.75) / 0.25
            x = x_min
            y = y_max - int(alpha * (y_max - y_min))

        draw_disk((x, y), neuron_radius, gray=label)
        label += 1

    # 4 inner-corner neurons, moved toward center by corner_offset
    c_tl = (x_min + corner_offset, y_min + corner_offset)
    c_tr = (x_max - corner_offset, y_min + corner_offset)
    c_bl = (x_min + corner_offset, y_max - corner_offset)
    c_br = (x_max - corner_offset, y_max - corner_offset)

    for center in [c_tl, c_tr, c_bl, c_br]:
        draw_disk(center, neuron_radius, gray=label)
        label += 1

    img.save(filename, format="TIFF")
    print(f"Saved {filename}")


if __name__ == "__main__":
    make_neuron_square_with_corners()