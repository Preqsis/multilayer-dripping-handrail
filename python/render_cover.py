import os
import h5py
import PIL
import cv2
import numpy as np
from argparse import ArgumentParser
from io import BytesIO
from render.Render import render_disc

def rotate(p, pitch, roll, yaw):
    pitch, roll, yaw = np.radians(pitch), np.radians(roll), np.radians(yaw)
    
    cosa, sina = np.cos(yaw), np.sin(yaw)
    cosb, sinb = np.cos(pitch), np.sin(pitch)
    cosc, sinc = np.cos(roll), np.sin(roll)
    
    Axx = cosa*cosb
    Axy = cosa*sinb*sinc - sina*cosc
    Axz = cosa*sinb*cosc + sina*sinc

    Ayx = sina*cosb
    Ayy = sina*sinb*sinc + cosa*cosc
    Ayz = sina*sinb*cosc - cosa*sinc
    
    Azx = -sinb
    Azy = cosb*sinc
    Azz = cosb*cosc
    
    X = Axx*p[0] + Axy*p[1] + Axz*p[2]
    Y = Ayx*p[0] + Ayy*p[1] + Ayz*p[2]
    Z = Azx*p[0] + Azy*p[1] + Azz*p[2]
    
    return X, Y, Z

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--sim_file", type=str)
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--frame", type=str)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument("--cmap", type=str, default="viridis")
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    with h5py.File(args.sim_file, "r") as f:
        data = f[f"d{args.frame}"][()]

    disc_surface = render_disc(
        data,
        data.shape[0], data.shape[1],
        w=args.width, h=args.height,
        return_surface=True,
        bg_rgb=(1., 1., 1.),
        cmap=args.cmap
    )

    img = np.ndarray(
        shape=(args.width, args.height, 4), 
        dtype=np.uint8, 
        buffer=disc_surface.get_data()
    )

    # input points
    # okraje obrazku
    a = (0, 0, 0)
    b = (0, img.shape[1], 0)
    c = (img.shape[0], 0, 0)
    d = (img.shape[0], img.shape[1], 0)
    input_pts = np.float32([a[:2], b[:2], c[:2], d[:2]])

    # output points
    #  pitch, roll, yaw = 145., -45., 0.
    pitch, roll, yaw = 145., -45., 0.
    A = rotate(a, pitch, roll, yaw)
    B = rotate(b, pitch, roll, yaw)
    C = rotate(c, pitch, roll, yaw)
    D = rotate(d, pitch, roll, yaw)
    output_pts = np.float32([A[:2], B[:2], C[:2], D[:2]])

    # translace
    # aby to neutikalo mimo obrazek
    tns = np.array((np.min(output_pts[:,0]), np.min(output_pts[:,1])))
    output_pts[:] -= tns

    output_pts[:]

    # velikost vysledneho obrazku
    # aby se do nej rotovany vesel
    h = int(np.max(output_pts[:,1]) - np.min(output_pts[:,1]))
    w = int(np.max(output_pts[:,0]) - np.min(output_pts[:,0]))

    # warp perspektivy
    M = cv2.getPerspectiveTransform(input_pts,output_pts)
    out = cv2.warpPerspective(img, M, (w, h), cv2.INTER_LINEAR, cv2.BORDER_CONSTANT,  1)

    # zapis vysledku
    cv2.imwrite(f"{args.outdir}/cover.png", out)

if __name__ == "__main__":
    main()
