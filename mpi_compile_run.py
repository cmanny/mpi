from __future__ import print_function
import os
import sys
import subprocess
import re
import time
#import itertools
import importlib
from lbm_structures import *

class HackyDict(dict):
  def __getitem__(self, key):
    try:
      return super(HackyDict, self).__getitem__(key)
    except KeyError:
      try:
        self[key] = importlib.import_module(key)
      except ImportError:
        pass
    return super(HackyDict, self).__getitem__(key)

class Pympi(object):
  
  def __init__(self):
    self.pympi_pat = "/\*!pympi\s*([\S\s]*?)\s*[^!/*]\*/"

  def compile(self, in_file, context):
    with open(in_file) as input_file:
      input_text = input_file.read()
      matches = re.findall(self.pympi_pat, input_text, re.MULTILINE)
      values = []
      for match in matches:
        try:
          values.append(str(eval(match, globals(), HackyDict(context))))
        except:
          print(match)
          raise
      output_text = input_text
      for val in values:
          output_text = re.sub(self.pympi_pat, val, output_text, count=1)
      out_file = "processes/src/" + str(context["r"].rank) + ".c"
      with open(out_file, 'w') as output_file:
          output_file.write(output_text)
      return out_file

class LBMRunner(object):

  def __init__(self, params_file, main_file, machine_file, num_processors, using_openmp, do_run):
    self.num_processors = int(num_processors)
    self.using_openmp = using_openmp
    self.main_file = main_file
    self.params_file = params_file
    self.ext = ".px"
    self.g_opts = dict()
    self.p_opts = []
    self.define_list = ["MPI_RANK=" + str(i) for i in range(0, self.num_processors)]
    self.machine_file = machine_file
    self.do_run = do_run == "true"

    try:
      os.makedirs("processes/src/")
    except OSError:
      pass
    try:
      os.makedirs("processes/bin/")
    except OSError:
      pass


  def run(self):
    """
    Parses, pymps, compiles and executes the processes using the arguments from the constructor
    """
    built = False
    try:
      with open("processes/src/build.txt", "r") as in_file:
        mod_time = int(in_file.readline().strip())
        params_file = in_file.readline().strip()
        params_mod_time = int(in_file.readline().strip())
        np = int(in_file.readline().strip())
      if mod_time >= int(os.path.getmtime(self.main_file)) and \
          params_file == self.params_file and \
          params_mod_time >= int(os.path.getmtime(self.params_file)) and \
          np == self.num_processors:
        print("Build valid, skipped.")
        built = True
    except IOError:
      pass
    if not built:
      self.parse_params()
      self.pympi_processes()
      self.compile_processes()
    if self.do_run:
      return self.spawn_mpi()

  def parse_params(self):
    """
    Parse params file and load the cells into the process source.
    """
    try:
      with open(self.params_file) as pf:
        self.parser = LBMVerticleDistribution(pf, self.num_processors)
        self.parser.build_cells()
    except IOError as e:
      print("Error parsing params file: " + e)

  def pympi_processes(self):
    self.parser.build_graph()
    self.parser.build_regions()
    #print self.parser.graph

  def compile_processes(self):
    try:
      os.remove("processes/src/build.txt")
    except OSError:
      pass
    self.pympi = Pympi()
    error_id = -1
    completed = 0
    comps = []
    print("Building source...\n[" + ("".join(" " for x in xrange(100))) + "]", end="\r")
    sys.stdout.flush()
    try:
      for ni in range(0, self.num_processors):
        try:
          context = dict({"p" : self.parser, "r" : self.parser.regions[ni]})
          compile_str = "mpicc -Ofast -xavx -qopt-report=5 -qopt-report-phase=vec -Wall -Wno-comment " + self.pympi.compile(self.main_file, context) + " -D " + self.define_list[ni] + " -o processes/bin/" + str(ni) + self.ext + " -lm"
          comps.append(subprocess.Popen(compile_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
          print("[" + ("".join(("="  if (ni >= x * self.num_processors / 100) else " ") for x in xrange(100))) + "] (%d/%d)" % (ni + 1, self.num_processors), end="\r")
          sys.stdout.flush()
        except:
          print("Could not build procedd %d" % ni)
          error_id = ni
          raise
      print("\nBuilt source, building programs")
      print("\nBuilding binaries...\n[" + ("".join(" " for x in xrange(100))) + "]", end="\r")
      try:
        iter = 0
        while completed < ni and iter < 50:
          to_remove = []
          for c in comps:
            if c.poll() is not None:
              if c.returncode != 0:
                raise RuntimeError
              completed += 1
              to_remove.append(c)
              print("[" + ("".join(("="  if (completed >= x * self.num_processors / 100) else " ") for x in xrange(100))) + "] (%d/%d)" % (completed, self.num_processors), end="\r")
              sys.stdout.flush()
          comps = [c for c in comps if c not in to_remove]
          time.sleep(0.1)
          iter += 1
        for c in comps:
          c.communicate()
          if c.poll() is not None:
            if c.returncode != 0:
              raise RuntimeError
            completed += 1
            to_remove.append(c)
            print("[" + ("".join(("="  if (completed >= x * self.num_processors / 100) else " ") for x in xrange(100))) + "] (%d/%d)" % (completed, self.num_processors), end="\r")
            sys.stdout.flush()

        print("\nBuild complete")
        #for ni in range(0, self.num_processors):
        #  print(comps[ni].stdout.read())
        #  print(comps[ni].stderr.read())
      except RuntimeError:
        print("Could not compile program")
        raise
    finally:
      for i, c in enumerate(comps):
        try:
          c.kill()
        except OSError:
          pass
        if i == error_id:
          print(comps[i].stderr.read())
          print(comps[i].stdout.read())
    # Generate compile report
    with open("processes/src/build.txt", "w") as out_file:
      out_file.write("\n".join(map(str, [
        int(os.path.getmtime(self.main_file)),
        self.params_file,
        int(os.path.getmtime(self.params_file)),
        self.num_processors,
        ])))

  def spawn_mpi(self):
    with open(self.machine_file) as mf:
      hosts = [n.rstrip("\n") for n in mf.readlines()]
    mpi_str = "mpirun " + " : ".join(["-host " + hosts[i] + " -np 1 processes/bin/" + str(i) + self.ext + " -x OMPI_COMM_WORLD_RANK=" + str(i) for i in range(0, self.num_processors)])
    subprocess.call(mpi_str, shell=True)

if __name__ == "__main__":
  if(len(sys.argv) < 5):
    print("Not enough arguments supplied, syntax is:\n\nmpi_compile_run.py [params_file] [mpi_template.c] [machine_file] [num_processors] [using_openmp(true/false)] [mpirun(true/false)]")
    exit()
  lbm = LBMRunner(*sys.argv[1:])
  lbm.run()
