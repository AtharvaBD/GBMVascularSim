from PIL import Image
import numpy as np

def make_three_endothelial_ring(
    filename="three_endothelial_ring.tif",
    nx=200, ny=200,
    R_ring=30,
    rad_cell=8,
    values=(80, 160, 240)
):
    """3 endothelial cells on a ring of radius R_ring."""
    cx, cy = nx // 2, ny // 2
    img = np.zeros((ny, nx), dtype=np.uint8)
    Y, X = np.indices((ny, nx))

    angles_deg = [0, 120, 240]
    for val, ang in zip(values, angles_deg):
        theta = np.deg2rad(ang)
        cx_i = int(cx + R_ring * np.cos(theta))
        cy_i = int(cy + R_ring * np.sin(theta))

        r = np.sqrt((X - cx_i)**2 + (Y - cy_i)**2)
        mask = r <= rad_cell
        img[mask] = val

    Image.fromarray(img, mode="L").save(filename)


from PIL import Image
import numpy as np

def make_pericytes_outside_ring(
    filename="pericytes_outside_ring.tif",
    nx=200, ny=200,
    R_ring=30,
    offset=60,
    rad_cell=5,
    values=(50, 100, 150, 200, 250),
    start_angle_deg=45
):
    """
    Five pericyte cells just outside the endothelial ring.

    R_ring: radius of endothelial ring
    offset: how far outside the ring (pixels)
    start_angle_deg: angle of the first pericyte; others are spaced evenly
    """
    cx, cy = nx // 2, ny // 2
    img = np.zeros((ny, nx), dtype=np.uint8)
    Y, X = np.indices((ny, nx))

    R_peri = R_ring + offset
    n_cells = 5
    angle_step = 360.0 / n_cells

    for i in range(n_cells):
        theta = np.deg2rad(start_angle_deg + i * angle_step)
        cx_p = int(cx + R_peri * np.cos(theta))
        cy_p = int(cy + R_peri * np.sin(theta))

        r = np.sqrt((X - cx_p)**2 + (Y - cy_p)**2)
        mask = r <= rad_cell
        img[mask] = values[i]

    Image.fromarray(img, mode="L").save(filename)



if __name__ == "__main__":
    # Toggle these calls as needed
    #make_three_endothelial_ring()          # generates three_endothelial_ring.tif
    make_pericyte_outside_ring()           # generates pericyte_outside_ring.tif
