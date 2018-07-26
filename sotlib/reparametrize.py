#!/usr/bin/python
# Modified by Zheng Ruan on 4 Oct 2017
#
# Grant Rotskoff, 12 July 2012
# Bjorn Wesen, June-Oct 2014
#
# Implement reparametrization method based on Luca Maragliano et al
# The Journal of Chemical Physics 125, 024106 (2006); doi: 10.1063/1.2212942
#

# This file is part of Copernicus
# http://www.copernicus-computing.org/
# 
# Copyright (C) 2011-2015, Sander Pronk, Iman Pouya, Grant Rotskoff, Bjorn Wesen, Erik Lindahl and others.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as published 
# by the Free Software Foundation
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# helper functions

def add(x, y): return x + y

def sub(x, y): return y - x

def scale(k, v): return [k*vi for vi in v]

def dist(v1, v2, norm):
    v = map(sub, v1, v2)
    return sum([(x/n)**2 for x, n in zip(v, norm)])**(0.5)

def L(n, path, norm):
    if n==0: 
        return 1
    else:
        pathlength = 0
        for i in range(n - 1):
            pathlength += dist(path[i], path[i + 1], norm)
        return pathlength

def s(m, path, norm):
    R = len(path) - 1
    return (m - 1) * L(R, path, norm) / (R - 1)

def dir(v1, v2, norm):
    normed = []
    d = dist(v1, v2, norm)
    for x in map(sub, v1, v2):
        normed += [x / d]
    return normed

# Reparametrize the points
def rep_pts(newpts):
    norm = []
    for i in range(newpts.shape[1]):
        coords = newpts[:,i]
        norm.append(max(coords)-min(coords))
    adjusted = [ newpts[0], newpts[ len(newpts) - 1 ] ]
    for i in range(2, len(newpts)): 
            k = 2
            while (L(k - 1, newpts, norm) >= s(i, newpts, norm) or s(i, newpts, norm) > L(k, newpts, norm)):
               k += 1
            v = dir(newpts[k - 2], newpts[k - 1], norm)
            reppt = (map(add, newpts[k - 2], scale((s(i, newpts, norm) - L(k - 1, newpts, norm)), v)))
            adjusted.insert(i - 1, reppt)
    return adjusted

if __name__ == '__main__':
    pass
