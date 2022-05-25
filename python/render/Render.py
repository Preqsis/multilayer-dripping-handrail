#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cairo
import IPython
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
from io import BytesIO, StringIO
import PIL

def zscale_linear(val, val_range):
    return val / val_range[1]

def zscale_exponential(val, limit):
    return 1 / (1 - np.exp(val / limit))

def zscale_log(val, val_range):
    return np.log10(val - val_range[0]) / np.log10(val_range[1] - val_range[0])

def render_disc(data, idim, jdim, w=1920, h=1080, cmap=None, r_in=0.1, r_out=0.45, 
        val_index=5, val_range=(0., 16.), bg_rgb=(0, 0, 0), 
        return_surface=False, skip_empty=False, log_norm=False,
        scf=1.0
        ):
    #colormap = mcol.LinearSegmentedColormap.from_list("custom", ["black", "blue", "yellow", "red"])
    #colormap = plt.cm.get_cmap("magma" if cmap is None else cmap)
    colormap = plt.cm.get_cmap("inferno" if cmap is None else cmap)

    val_range = (val_range[0]*scf, val_range[1]*scf)

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

    c.set_line_width(1.02 * dr)
    
    dphi = 2. * np.pi / jdim

    r = r_out
    for i in range(idim):
        for j in range(jdim):
            #k = i * jdim + j

            # value_index == 5 --> mass
            val, azm = data[i][j][val_index] * scf, data[i][j][8] % (2. * np.pi)

            if val == 0. and skip_empty:
                continue

            rgba = colormap(zscale_linear(val, val_range) if not log_norm else zscale_log(val, val_range))

            c.arc(0.495, 0.5, r, azm, azm+dphi*1.02)
            c.set_source_rgb(rgba[0], rgba[1], rgba[2])  # Solid color
            c.stroke()

        r -= dr

    c.set_line_width(dr * 0.15)
    c.arc(0.495, 0.5, r_out + 0.6 * dr, 0., 2. * np.pi)
    c.set_source_rgb(0, 0, 0)
    c.stroke()

    c.arc(0.495, 0.5, r_in + 0.4 * dr, 0., 2. * np.pi)
    c.set_source_rgb(0, 0, 0)
    c.stroke()

    if return_surface:
        return surface

    buff = BytesIO()
    surface.write_to_png(buff)
    imdata = buff.getvalue()
    buff.close()

    surface.finish()

    return imdata

def render_frame(data_disc, data_obs, idim, jdim, dt, w=1920, h=1080, cmap=None, r_in=0.1, r_out=0.45, 
        val_index=5, val_range=(0., 16.), bg_rgb=(0, 0, 0), 
        return_surface=False, skip_empty=False, log_norm=False,
        scf=1.0, lcdepth=200, lc_ylim=None, dpi=150
        ):
    fig, axes = plt.subplots(nrows=2, ncols=1)

    print(idim, jdim)

    fig.set_size_inches(885 / dpi, h / dpi)
    fig.set_dpi(dpi)

    fig.patch.set_facecolor('#000000')

    axes[0].plot(data_obs[:,0] * dt, data_obs[:,2])
    axes[0].set_xlabel("$t\ (s)$", color="#ffffff")
    axes[0].set_ylabel("$L_{total}\ (erg \cdot s^{-1})$", rotation=90, color="#ffffff")
    axes[0].set_ylim(lc_ylim)
    axes[0].set_xlim((data_obs[:,0].max() * dt - lcdepth * dt, data_obs[:,0].max() * dt))

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
