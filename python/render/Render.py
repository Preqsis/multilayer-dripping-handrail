#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cairo
import IPython
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import PIL

def render_disc(data, idim, jdim, w=1500, h=1500, cmap=None, r_in=0.1, r_out=0.45, mlimit=16., bg_rgb=(0, 0, 0), return_surface=False):
    colormap = plt.cm.get_cmap("magma" if cmap is None else cmap)
    dr = (r_out - r_in) / idim

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
    c = cairo.Context(surface)
    c.scale(h, h)

    # set bacground color (black)
    #c.set_source_rgb(0, 0, 0)
    c.set_source_rgb(*bg_rgb)
    c.paint()
    c.save()
    c.restore()

    c.set_line_width(dr)
    
    dphi = 2. * np.pi / jdim

    r = r_out
    for i in range(idim):
        for j in range(jdim):
            k = i * jdim + j

            m, azm = data[k][5], data[k][9] % (2. * np.pi)

            rgba = colormap(m / mlimit)

            c.arc(0.5, 0.5, r, azm, azm+dphi)
            c.set_source_rgb(rgba[0], rgba[1], rgba[2])  # Solid color
            c.stroke()

        r -= dr

    c.set_line_width(dr * 0.2)
    c.arc(0.5, 0.5, r_out + 0.6 * dr, 0., 2. * np.pi)
    c.set_source_rgb(255, 255, 255)
    c.stroke()

    c.arc(0.5, 0.5, r_in - 0.6 * dr, 0., 2. * np.pi)
    c.set_source_rgb(255, 255, 255)
    c.stroke()

    if return_surface:
        return surface

    buff = BytesIO()
    surface.write_to_png(buff)
    imdata = buff.getvalue()
    buff.close()

    surface.finish()

    return imdata

def render_frame(data_disc, idim, jdim, w=1500, h=1500, cmap=None, r_in=0.1, r_out=0.45, mlimit=16., dpi=150):
    tmp = np.empty((100, 2))
    tmp[:,0] = np.arange(0, tmp.shape[0], 1)
    tmp[:,1] = np.random.rand(tmp.shape[0], 1)[:,0]

    buff = BytesIO()

    fig, axes = plt.subplots(nrows=2, ncols=1)

    fig.set_size_inches(885 / dpi, h / dpi)
    fig.set_dpi(dpi)

    fig.patch.set_facecolor('#000000')

    axes[0].plot(tmp[:,0], tmp[:,1])
    axes[1].plot(tmp[:,0], tmp[:,1])

    for ax in axes:
        ax.patch.set_facecolor('#000000')
        ax.tick_params(color='#ffffff', labelcolor='#ffffff')
        for spine in ax.spines.values():
            spine.set_edgecolor('#ffffff')

    fig.canvas.draw()

    cw, ch = fig.canvas.get_width_height()
    img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    plt.close()

    disc_surface = render_disc(data_disc, idim, jdim, r_in=r_in, r_out=r_out, w=w, h=h, return_surface=True)

    buff = BytesIO()
    disc_surface.write_to_png(buff)
    a = np.array(PIL.Image.open(buff))
    buff.close()
    disc_surface.finish()

    col = np.arange(0, w, 1)
    m = col >= w - cw

    b = np.copy(img)
    a[:,m] = b

    return PIL.Image.fromarray(a)
