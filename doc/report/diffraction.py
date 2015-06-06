from pylab import *
from astropy.modeling import models
from matplotlib.colors import LogNorm
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = 'CMU Sans Serif'
rcParams['font.serif'] = 'CMU Serif'
rcParams['font.monospace'] = 'CMU Typewriter'

figure(figsize=(11,4))
a = models.AiryDisk2D(1, 0, 0, 50)

def norm_z(z):
    return z**0.38

subplot(131)

y, x = mgrid[-150:150, -150:150]
z = a(x, y)

im = imshow(norm_z(z), interpolation='bicubic', cmap='gray')
xticks([])
yticks([])

subplot(132)

a = models.AiryDisk2D(1, -25, 0, 50)
b = models.AiryDisk2D(1, 25, 0, 50)

z = (a + b)(x, y)
im = imshow(norm_z(z), interpolation='bicubic', cmap='gray')
xticks([])
yticks([])

subplot(133)

a = models.AiryDisk2D(1, -15, 0, 50)
b = models.AiryDisk2D(1, 15, 0, 50)

z = (a + b)(x, y)
im = imshow(norm_z(z), interpolation='bicubic', cmap='gray')
xticks([])
yticks([])
tight_layout()
#savefig('figures/airy.pdf')
show()
