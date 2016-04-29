import numpy
from direction import *

class Node(object):

  def __init__(self, pos):
    self.pos = pos
    self.neighbours = [None for _ in range(Drn.num)]
    self.rank = LBMGraph.neutral_rank

  def __getitem__(self, key):
    return self.neighbours[key]

  def __repr__(self):
    return str(self.pos) + " -> " + str([n.pos if n is not None else None for n in self.neighbours])

  def add_neighbour(self, drn, node):
    self.neighbours[drn] = node
    node.neighbours[drn_opposite(drn)] = self

  def del_neighbour(self, drn):
    nb = self.neighbours[drn]
    if nb != None:
      self.neighbours[drn].neighbours[drn_opposite(drn)] = None
      self.neighbours[drn] = None


  def cut_dist(edge1, edge2):
    return numpy.linalg.norm(numpy.array(edge1) - numpy.array(edge2))

class LBMGraph(object):

  neutral_rank = 65

  lrank = neutral_rank
  rrank = neutral_rank

  def __init__(self, rank = neutral_rank, x = 0, y = 0, w = 0, h = 0):
    self.nodes = dict()
    self.w = w
    self.h = h
    self.rank = rank

    self.x = x
    self.y = y
    self.w = w
    self.h = h

    self.right_halos = list()
    self.left_halos = list()

  def __getitem__(self, key):
    try:
      node = self.nodes[key]
    except KeyError:
      node = self.nodes[key] = Node(key)
    return node

  def has(self, key):
    return key in self.nodes

  def __repr__(self):
    return "|\n".join([" ".join([chr(40 + self.nodes[(x,y)].char) if (x,y) in self.nodes else " " \
                          for x in range(self.x, self.x + self.w)]) \
                          for y in range(self.y + self.h - 1, self.y - 1, - 1)]) \


  def add(self, node):
    self.nodes[node.pos] = node
    node.rank = self.rank
    node.char = self.rank
    return node

  def remove(self, node):
    try:
      for i in range(Drn.num):
        node.del_neighbour(i)
      self.nodes.pop(node.pos)
    except KeyError:
      pass

  def add_all(self, list):
    for node in list:
      self.add(node)

  def add_vert_halo(self, x, drn, ystart, yend):
    try:
      col = self.halos[y]


    except KeyError:
      self.halos[y] = set(node.pos)
    node.char = -2

  def timestep_code(self):
    return "return;"

class LBMData(object):

  def __init__(self, pf, num_ranks):
    self.box_size = 100.0
    self.w = int(pf.readline())
    self.h = int(pf.readline())
    self.num_iters = int(pf.readline())
    self.reynolds_dim = int(pf.readline())
    self.fluid_density = float(pf.readline())
    self.accel = float(pf.readline()) 
    self.omega = float(pf.readline())
    _, self.dim, self.idx = pf.readline().split()
    self.idx = int((self.w if self.dim == "column" else self.h) * int(self.idx)/self.box_size)
    self.n_obs = int(pf.readline()[0])
    self.objects = [map(int, pf.readline().split()) for _ in range(self.n_obs)]

    self.graph = LBMGraph(138, 0, 0, self.w, self.h)
    self.regions = []
    self.num_ranks = num_ranks

    self.w0 = self.fluid_density * 4 / 9
    self.w1 = self.fluid_density / 9
    self.w2 = self.fluid_density / 36

  def speed(self, j, size, id):
    return j * size + id

  def in_object(self, x, y, object_):
    if x >= object_[0] and x < object_[2] and y >= object_[1] and y < object_[3]:
      return True
    return False

  def obs_type(self, x, y):
    for s in range(x - 1, x + 2):
      for t in range(y - 1, y + 2):
        s = s % self.w
        t = t % self.h
        if self.obstacles[t * self.w + s] == 0:
          self.tot_edges += 1
          return 2
    return 1

  def build_cells(self):
    size = self.w * self.h
    self.obstacles = list([0]) * self.w * self.h
    cells = list([0]) * self.w * self.h * 9

    self.tot_cells = 0
    self.tot_edges = 0
    for y in range(self.h):
      for x in range(self.w):
        idi = y * self.w + x
        cells[idi] = self.w0
        for j in range(1, 5): cells[self.speed(j, size, idi)] = self.w1
        for j in range(5, 9): cells[self.speed(j, size, idi)] = self.w2
        self.obstacles[idi] = 0
    for x in range(self.w):
      for y in range(self.h):
        for k in range(self.n_obs):
          x_pos = x * self.box_size / self.w
          y_pos = y * self.box_size / self.h
          if self.in_object(x_pos, y_pos, self.objects[k]):
            self.obstacles[y * self.w + x] = 1
        if not self.obstacles[y * self.w + x]:
          self.tot_cells += 1
    for y in range(self.h - 1, - 1, - 1):
      for x in range(self.w):
        if self.obstacles[y * self.w + x] == 1:
          self.obstacles[y * self.w + x] = self.obs_type(x, y)
    self.old_obstacles = list(self.obstacles)
    self.cells_per_region = (self.tot_cells + self.tot_edges) / self.num_ranks
  
  def build_graph(self):
    for y in range(self.h):
      for x in range(self.w):
        if self.old_obstacles[y * self.w + x] in [0, 2]:
          node = self.graph[(x,y)]
          for drn in range(Drn.num):
            opos = drn_move(x, y, drn)
            opos = (opos[0] % self.w, opos[1] % self.h)
            if self.old_obstacles[opos[1] * self.w + opos[0]] in [0, 2]:
              other = self.graph[opos]
              node.add_neighbour(drn, other)

  def __repr__(self):
    output = ""
    for y in range(self.h - 1, - 1, - 1):
      for x in range(self.w):
        output += str(self.obstacles[y * self.w + x]) + " "
      output += "\n"
    return output

class LBMVerticleDistribution(LBMData):
  
  def make_column(self, x):
    column = []
    ys = -1
    h = 0
    for y in range(self.h):
      if self.graph.has((x,y)):
        if ys == -1:
          ys = y
        column.append(self.graph[(x, y)])
        h += 1
    return column, ys, h 

  def make_right_halos(self, region, column):
    for drn in [Drn.ne, Drn.e, Drn.se]:
      send_min = self.h
      recv_min = self.h
      cnt = 0
      for node in column:
        if node[drn] and node[drn].rank not in [LBMGraph.neutral_rank, region.rank]:
          region.rrank = node[drn].rank
          self.regions[region.rrank].lrank = region.rank
          send_min = min(send_min, node.pos[1])
          recv_min = min(recv_min, node[drn].pos[1])
          cnt += 1
      send = dict({'drn': drn_real(drn), 'x': region.w, 'y': send_min - (region.y - 1), 'cnt': cnt, 'rank' : region.rrank})
      recv = dict({'drn': drn_real(drn_opposite(drn)), 'x': region.w + 1, 'y': recv_min - (region.y - 1), 'cnt': cnt, 'rank' : region.rank})
      if cnt != 0:
        region.right_halos.append((send, recv))

  def make_left_halos(self, region, column):
    for drn in [Drn.sw, Drn.w, Drn.nw]:
      send_min = self.h
      recv_min = self.h
      cnt = 0
      for node in column:
        if node[drn] and node[drn].rank not in [LBMGraph.neutral_rank, region.rank]:
          region.lrank = node[drn].rank
          self.regions[region.lrank].rrank = region.rank
          send_min = min(send_min, node.pos[1])
          recv_min = min(recv_min, node[drn].pos[1])
          cnt += 1
      send = dict({'drn': drn_real(drn), 'x': 1, 'y': send_min - (region.y - 1), 'cnt': cnt, 'rank' : region.lrank})
      recv = dict({'drn': drn_real(drn_opposite(drn)), 'x': 0, 'y': recv_min - (region.y - 1), 'cnt': cnt, 'rank' : region.rank})
      if cnt != 0:
        region.left_halos.append((send, recv))

  def build_regions(self):
    rank = claimed = init_height = tot_claimed = x = 0
    region = LBMGraph(rank)
    self.regions.append(region)
    enforce_rect = True
    while x < self.w:
      column, ys, h = self.make_column(x)
      next_rank = False
      bad_col = False
      if init_height == 0:
        init_height = len(column)
        region.x = x
        region.y = ys
        region.h = h
      else:
        if len(column) != init_height and enforce_rect:
          next_rank = True
          bad_col = True
      if not bad_col:
        region.add_all(column)
        if region.w == 0:
          self.make_left_halos(region, column)
        region.w += 1
        claimed += len(column)
        tot_claimed += len(column)
        if rank != self.num_ranks - 1:
          self.cells_per_region = 1 + (self.tot_cells + self.tot_edges - tot_claimed) / (self.num_ranks - (rank + 1))
        x += 1
      if claimed > self.cells_per_region:
        next_rank = True
      if next_rank or x == self.w:
        if rank != 0:
          self.make_right_halos(self.regions[rank-1], self.make_column(x-region.w-1)[0])
        init_height = 0
        claimed = 0
        rank += 1
        if x != self.w:
          region = LBMGraph(rank)
          self.regions.append(region)
        else:
          self.make_right_halos(region, column)
          self.make_left_halos(self.regions[0], self.make_column(0)[0])
    for r in self.regions:
      # print r.rank
      # print r.lrank
      # print r.rrank
      # print "lefts\n" + "\n".join(map(str, [h for h in r.left_halos]))
      # print "rihts\n" + "\n".join(map(str, [h for h in r.right_halos]))
      # print "\n"
      r.x = r.x - 1
      r.y = r.y - 1
      r.w = r.w + 2
      r.h = r.h + 2

    return
