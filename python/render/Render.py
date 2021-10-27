#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cairo
import IPython
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO, StringIO

def render_disc(data, idim, jdim, w=1500, h=1500, cmap=None, r_in=0.1, r_out=0.45, mlimit=16.):
    colormap = plt.cm.get_cmap("magma" if cmap is None else cmap)
    dr = (r_out - r_in) / idim

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
    c = cairo.Context(surface)
    c.scale(w, h)

    # set bacground color (black)
    #c.set_source_rgb(0, 0, 0)
    c.set_source_rgb(255, 255, 255)
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

    buff = BytesIO()
    surface.write_to_png(buff)
    imdata = buff.getvalue()
    buff.close()

    surface.finish()

    return imdata
