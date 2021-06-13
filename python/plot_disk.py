#!/usr/bin/env python

import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy import ndimage
from scipy.interpolate import Rbf

import argparse
import os
import sys

from multiprocessing import Process

def AddLeadingSpaces(number, n=2):
    tmp = str(number)
    while len(tmp) < n:
        tmp = " "+tmp
    return tmp

def create_circular_mask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask

def plot(data, outfile, t=None):
    mlimit = 8.0
    cnorm = mpl.colors.Normalize(0, mlimit)
    colormap = plt.get_cmap("YlOrRd")
    #colormap = plt.get_cmap("Greys")
    idim, jdim = int(data[-1][0]+1), int(data[-1][1]+1)

    # 'podklad'
    rad = np.linspace(0, 1, idim+1)
    azm = np.linspace(0, 2*np.pi, num = jdim+1)
    r, th = np.meshgrid(rad, azm)
    
    #fig = plt.figure()
    #ax = Axes3D(fig)
    #fig.set_size_inches(15, 15)
    
    z = np.empty((idim, jdim), dtype="f8")
    for i in range(idim):
        z[i] = data[i*idim:i*idim+jdim,5]
    z = np.flip(np.transpose(z), 1)
        
    fig = plt.figure(figsize=(10, 10), dpi=200)
    ax = fig.add_subplot(111, polar=True)
    
    ax.pcolormesh(th, r, z, norm=cnorm, cmap=colormap, shading="auto")
    ax.set_rorigin(-1.0)
    #plt.colorbar()
    
    ax.set_xticks([])
    ax.set_yticks([])

    """
    years = int(np.floor(t / (365 * 86400)))
    days = int(np.floor((t - years * 365 * 86400) / 86400))
    hours = int(np.floor((t - years * 365 * 86400 - days * 86400) / 3600))
    """

    ax.set_title("%d" % (t), loc="right")
    ax.plot(azm, r, color=colormap, ls='none')

    plt.savefig(outfile)
    plt.close('all')

def plot2(data, outfile, t=None, mask_radius=None, cmap=None, norm=None, grid_offset=None, grid_step=None, value_threshold=None, output_size=None):
    grid_offset = 0. if grid_offset is None else grid_offset
    grid_step = 1.0 if grid_step is None else grid_step
    cmap = "jet" if cmap is None else cmap
    value_threshold = 0.2 if value_threshold is None else value_threshold
    output_size = (1000, 1000) if output_size is None else output_size

    # prislusne datove sloupce
    radius = 25 - 1 - data[:,0]
    azimut = data[:,9]
    values = data[:,5]
    values[values <= value_threshold] = 0. # vynulovat pod urovni thresholdu

    # do kaztezkych
    x = radius * np.cos(azimut)
    y = radius * np.sin(azimut)

    # rozsah pro grid
    x_min, x_max = np.amin(x), np.amax(x)
    y_min, y_max = np.amin(y), np.amax(y)

    # grid
    grid_x, grid_y = np.mgrid[x_min-grid_offset:x_max+grid_offset:grid_step, y_min-grid_offset:y_max+grid_offset:grid_step]

    # iterpolator
    rbfi = Rbf(x, y, values, smooth=1)
    img = rbfi(grid_x, grid_y)

    #di = np.flip(di, axis=0)
    #di = ndimage.rotate(di, 270)
    #h, w = di.shape
    #fm = np.logical_not(create_circular_mask(h, w))

    if mask_radius is not None:
        h, w = img.shape
        m = np.logical_not(create_circular_mask(h, w, radius=mask_radius))
        img[m] = 0.0

    # zoom na pozadovanou velikost
    img = ndimage.zoom(img, output_size[0] / img.shape[0], mode="nearest")

    img = cmap(norm(img))

    #plt.figure(figsize=(15, 15))
    plt.imsave(outfile, img)
    plt.close("all")

def worker(data, frames, outdir, dt, mask_radius=None, cmap=None, norm=None, grid_offset=None, grid_step=None, value_threshold=None, output_size=None):

    for frame in frames:
        t = int(frame.split("_")[1]) * dt if dt is not None else None
        print("plot: ", frame)
        outfile = "%s/%s.png" % (outdir, frame)
        plot2(data[frame], outfile, t=t, mask_radius=mask_radius, cmap=cmap, norm=norm, grid_offset=grid_offset, grid_step=grid_step, value_threshold=value_threshold, output_size=output_size)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=str, help="HDF5 data file")
    p.add_argument("--output", default="/tmp", type=str, help="Plots output dir.")
    p.add_argument("--clear", action="store_true", help="If set => clear plots output directory")
    p.add_argument("--nproc", type=int, default=4, help="Number of working processes")
    p.add_argument("--nth_frame", type=int, default=1, help="Plot every nth frame (default=100)")
    p.add_argument("--first_frame", type=int, help="First frame to plot")
    p.add_argument("--last_frame", type=int, help="Last frame to plot")
    p.add_argument("--skip_plot", action="store_true", help="If set => skips PNGs ploting")
    p.add_argument("--create_video", action="store_true", help="If set => creates video out of PNGs")
    p.add_argument("--video_framerate", type=int, default=20, help="Sets outputed video framerate (defult=20)")
    args = p.parse_args()

    # vytvorit pokud neexistuje
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # pokud zadano tak smazat obsah plot adresare
    if args.clear and not args.skip_plot:
        os.system("rm -r %s/*" % args.output)

    f =  h5py.File(args.data, "r")
    keys, data = list(f.keys()), {}
    first_frame = 0 if args.first_frame is None else args.first_frame
    last_frame = len(keys) if args.last_frame is None else args.last_frame+1

    idim = 50
    mlimit = 16.0
    mask_radius = idim * 2 + 2
    cmap = plt.get_cmap("magma")
    norm = mpl.colors.Normalize(0, mlimit)
    grid_offset = 5
    grid_step = 0.5
    value_threshold = 0.2
    output_size = (1000, 1000)
    dt = 4847

    if not args.skip_plot:
        batch_size = 40
        f0 = first_frame
        while f0 < last_frame:
            data = {}
            frames = []
            for frame in range(f0, f0+batch_size*args.nth_frame, args.nth_frame):
                if "data_%d" % frame in keys:
                    data["data_%06d" % frame] = f["data_%d" % frame][()]
                    frames.append("data_%06d" % frame)
            
            a, b = int(len(frames) / args.nproc), int(len(frames) % args.nproc)

            jobs = []
            for i in range(args.nproc):
                job_frames = frames[i*a:i*a+a]

                job_data = {}
                for key in job_frames:
                    job_data[key] = data[key] 

                jobs.append([job_frames, job_data])
        
            if b > 0:
                job_frames = frames[-b:]
                jobs[-1][0] += job_frames
                for key in job_frames:
                    jobs[-1][1][key] = data[key] 

            pcs = []
            for j in jobs:
                pcs.append(Process(target=worker, args=(j[1], j[0], args.output, dt, mask_radius, cmap, norm, grid_offset, grid_step, value_threshold, output_size)))

            for p in pcs:
                p.start()

            for p in pcs:
                p.join()

            f0 += batch_size * args.nth_frame

    if args.create_video:
        os.system("ffmpeg -framerate %d -pattern_type glob -i '%s/*.png' -c:v ffv1 %s/disk_animation.avi" % (args.video_framerate, args.output, args.output))

