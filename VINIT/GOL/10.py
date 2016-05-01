# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 01:43:58 2016

@author: Vinit Dedhiya
"""

size = np.array(Z.shape)
dpi = 72.0
figsize= size[1]/float(dpi),size[0]/float(dpi)
fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
plt.imshow(Z,interpolation='nearest', cmap=plt.cm.gray_r)
plt.xticks([]), plt.yticks([])
plt.show()