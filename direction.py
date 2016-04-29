from enum import IntEnum

class Drn(IntEnum):
  e = 0
  n = 1
  w = 2
  s = 3
  ne = 4
  nw = 5
  sw = 6
  se = 7

  num = 8

def drn_opposite(drn):
  return ((drn + 2) % 4) + 4*(drn > 3)

def drn_real(drn):
  return (drn + 1)

def drn_move(x, y, drn):
  x += (-1)**(drn in [Drn.w, Drn.nw, Drn.sw]) if drn not in [Drn.n, Drn.s] else 0
  y += (-1)**(drn in [Drn.s, Drn.se, Drn.sw]) if drn not in [Drn.e, Drn.w] else 0
  return (x, y)