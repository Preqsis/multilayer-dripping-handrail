#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cairo
import IPython
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import PIL

def render_disc(data, idim, jdim, w=1920, h=1080, cmap=None, r_in=0.1, r_out=0.45, mlimit=16., bg_rgb=(0, 0, 0), return_surface=False):
    colormap = plt.cm.get_cmap("magma" if cmap is None else cmap)
    dr = (r_out - r_in) / idim

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
    c = cairo.Context(surface)
    c.scale(h, h)

    # set bacground color (black)
    #c.set_source_rgb(255, 0, 0)
    c.set_source_rgb(*bg_rgb)
    c.paint()
    c.save()
    c.restore()

    c.set_line_width(dr)
    
    dphi = 2. * np.pi / jdim

    r = r_out
    for i in range(idim):
        for j in range(jdim):
            #k = i * jdim + j

            m, azm = data[i][j][5], data[i][j][9] % (2. * np.pi)

            rgba = colormap(m / mlimit)

            c.arc(0.495, 0.5, r, azm, azm+dphi)
            c.set_source_rgb(rgba[0], rgba[1], rgba[2])  # Solid color
            c.stroke()

        r -= dr

    c.set_line_width(dr * 0.2)
    c.arc(0.495, 0.5, r_out + 0.6 * dr, 0., 2. * np.pi)
    c.set_source_rgb(255, 255, 255)
    c.stroke()

    c.arc(0.495, 0.5, r_in + 0.4 * dr, 0., 2. * np.pi)
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

def render_frame(data_disc, data_lc, idim, jdim, w=1920, h=1080, cmap=None, r_in=0.1, r_out=0.45, mlimit=16., dpi=150, lc_plot_range=(0., 1.)):
    fig, axes = plt.subplots(nrows=2, ncols=1)

    fig.set_size_inches(885 / dpi, h / dpi)
    fig.set_dpi(dpi)

    fig.patch.set_facecolor('#000000')

    axes[0].plot(data_lc[:,0], data_lc[:,1])
    axes[0].set_xlabel("STEP", color="#ffffff")
    axes[0].set_ylabel("LUMINOSITY", rotation=90, color="#ffffff")
    axes[0].set_ylim(lc_plot_range)

    for ax in axes:
        ax.patch.set_facecolor('#000000')
        ax.tick_params(color='#ffffff', labelcolor='#ffffff')
        for spine in ax.spines.values():
            spine.set_edgecolor('#ffffff')

    fig.canvas.draw()

    iw, ih = fig.get_size_inches()
    img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    #img = img.reshape(fig.canvas.get_width_height()[::-1]*dpi + (3,))
    img = img.reshape((int(ih*dpi), int(iw*dpi), 3))

    plt.close()

    disc_surface = render_disc(data_disc, idim, jdim, r_in=r_in, r_out=r_out, w=w, h=h, return_surface=True)

    buff = BytesIO()
    disc_surface.write_to_png(buff)
    a = np.array(PIL.Image.open(buff))
    buff.close()
    disc_surface.finish()

    col = np.arange(0, w, 1)
    m = col >= w - int(iw * dpi)

    b = np.copy(img)
    a[:,m] = b

    return PIL.Image.fromarray(a)
