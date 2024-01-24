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


def solve_no_trig_vectorized(semi_major, semi_minor, p):
    px = np.abs(p[:, 0]).reshape([-1, 1])
    py = np.abs(p[:, 1]).reshape([-1, 1])

    npts = p.shape[0]

    tx = np.full((npts, 1), 0.707106)
    ty = np.full((npts, 1), 0.707106)

    a = semi_major
    b = semi_minor

    for iter in range(0, 3):
        x = a * tx
        y = b * ty

        ex = (a * a - b * b) * tx ** 3 / a
        ey = (b * b - a * a) * ty ** 3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = np.hypot(ry, rx)
        q = np.hypot(qy, qx)

        tx = np.clip((qx * r / q + ex) / a, 0, 1)
        ty = np.clip((qy * r / q + ey) / b, 0, 1)
        t = np.hypot(tx, ty)
        tx /= t
        ty /= t

    x = a * tx
    y = b * ty

    return np.hstack((np.copysign(x, p[:, 0].reshape([-1, 1])), np.copysign(y, p[:, 1].reshape([-1, 1]))))


def solve(semi_major, semi_minor, p):
    px = np.abs(p[:, 0]).reshape([-1, 1])
    py = np.abs(p[:, 1]).reshape([-1, 1])

    t = np.full((p.shape[0], 1), math.pi / 4)

    a = semi_major
    b = semi_minor

    for iter in range(0, 3):
        x = a * np.cos(t)
        y = b * np.sin(t)

        ex = (a * a - b * b) * np.cos(t) ** 3 / a
        ey = (b * b - a * a) * np.sin(t) ** 3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = np.hypot(ry, rx)
        q = np.hypot(qy, qx)

        delta_c = r * np.arcsin((rx * qy - ry * qx) / (r * q))
        delta_t = delta_c / np.sqrt(a * a + b * b - x * x - y * y)

        t += delta_t
        t = np.clip(t, 0, math.pi / 2)

    return np.hstack((np.copysign(x, p[:, 0].reshape([-1, 1])), np.copysign(y, p[:, 1].reshape([-1, 1]))))

def plot(semi_major, semi_minor):
    xs = np.linspace(-plot_size[0]/2, (plot_size[0]/2)-1, plot_size[0])
    ys = np.linspace(-plot_size[1]/2, (plot_size[1]/2)-1, plot_size[1])

    xxs, yys = np.meshgrid(xs, ys)
    pts = np.hstack((xxs.reshape([-1,1]), yys.reshape([-1,1])))
    p = solve_no_trig_vectorized(semi_major, semi_minor, pts)
    #p = solve(semi_major, semi_minor, pts)

    return np.sqrt(np.sum((pts - p) ** 2, axis=1)).reshape([plot_size[1],plot_size[0]]).T

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
