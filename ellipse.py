import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches

plot_size = (640, 480)
plots = [
    (240, 120),
    (240, 200),
    (240,  40),
    (240,   1),
    (  1, 120),
    (120, 240),
    (  1,   2),
    (320, 240),
    (100, 100)]

def solve(semi_major, semi_minor, p):
    px = abs(p[0])
    py = abs(p[1])

    t = math.pi / 4

    a = semi_major
    b = semi_minor

    for x in range(0, 3):
        x = a * math.cos(t)
        y = b * math.sin(t)

        ex = (a*a - b*b) * math.cos(t)**3 / a
        ey = (b*b - a*a) * math.sin(t)**3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)

        delta_c = r * math.asin((rx*qy - ry*qx)/(r*q))
        delta_t = delta_c / math.sqrt(a*a + b*b - x*x - y*y)

        t += delta_t
        t = min(math.pi/2, max(0, t))

    return (math.copysign(x, p[0]), math.copysign(y, p[1]))


xs = np.linspace(-plot_size[0]/2, (plot_size[0]/2)-1, plot_size[0])
ys = np.linspace(-plot_size[1]/2, (plot_size[1]/2)-1, plot_size[1])

def plot(semi_major, semi_minor):
    dist = np.zeros(plot_size)

    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            p = solve(semi_major, semi_minor, (x, y))
            d = math.hypot(y-p[1], x-p[0])
            dist[ix, iy] = d

    return dist

fig, axes = plt.subplots(nrows=3, ncols=3)
for i, ax in enumerate(axes.flat):
    a, b = plots[i]
    dist = plot(a, b)
    title = 'a = {0}, b = {1}'.format(a, b)
    ellipse = matplotlib.patches.Ellipse(xy=(plot_size[0]/2, plot_size[1]/2), width=2*a, height=2*b, angle=0,
                    edgecolor='r', fc='None', linestyle='dashed')

    im = ax.imshow(dist.T)
    im.set_cmap('jet')
    ax.add_patch(ellipse)
    im.axes.set_title(title)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()
