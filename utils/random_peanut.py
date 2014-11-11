# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

# Utility code for building random initial configuration on a peanut

from datetime import *
from random import uniform 
from math import *
import argparse
import numpy as np

from particle import *

def projection_matrix(axis):
  n = axis/np.sqrt(np.dot(axis,axis))
  return (np.identity(3) - np.outer(n,n))


def rotation_matrix(axis,theta):
  n = axis/np.sqrt(np.dot(axis,axis))
  a = np.cos(theta/2)
  b,c,d = -n*np.sin(theta/2)
  return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

class Peanut:
  
  def __init__(self, a, b, N, v):
    self.a = a
    self.b = b
    self.N = N
    self.particles = [Particle(i) for i in range(N)]
    self.__generate_pos()
    self.__generate_vel(v)
    self.__generate_director()
    
  def __generate_pos(self):
    for i in range(self.N):
      u = uniform(0,pi)
      v = uniform(0,2*pi)
      f = sqrt(self.a**2*cos(2*u)+sqrt(self.b**4-self.a**4*(sin(2*u))**2))
      x = cos(u)*f
      y = cos(v)*sin(u)*f
      z = sin(v)*sin(u)*f
      self.particles[i].r = [x,y,z]

  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      x, y, z = p.r
      fact_1 = 4.0*self.a*self.a;
      fact_2 = 4.0*(x*x + y*y + z*z);
      Nx = (-fact_1 + fact_2)*x;
      Ny = ( fact_1 + fact_2)*y; 
      Nz = ( fact_1 + fact_2)*z;
      axis = np.array([Nx,Ny,Nz])
      v = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      v = np.dot(projection_matrix(axis),v)
      vlen = sqrt(sum(v**2))
      v *= vav/vlen
      theta = uniform(0,2*pi)
      v = np.dot(rotation_matrix(axis,theta),v)
      p.v = v
 
  def __generate_director(self):
    for p in self.particles:
      x, y, z = p.r
      fact_1 = 4.0*self.a*self.a;
      fact_2 = 4.0*(x*x + y*y + z*z);
      Nx = (-fact_1 + fact_2)*x;
      Ny = ( fact_1 + fact_2)*y; 
      Nz = ( fact_1 + fact_2)*z;
      axis = np.array([Nx,Ny,Nz])
      n = np.array([uniform(-0.5,0.5), uniform(-0.5,0.5), uniform(-0.5,0.5)])
      n = np.dot(projection_matrix(axis),n)
      nlen = sqrt(sum(n**2))
      n *= 1.0/nlen
      theta = uniform(0,2*pi)
      n = np.dot(rotation_matrix(axis,theta),n)
      p.n = n
 
  def write(self,outfile):
    gentime = datetime.now()
    out = open(outfile,'w')
    out.write('# Total of %d particles\n' % self.N)
    out.write('# Generated on : %s\n' % str(gentime))
    out.write('# id  type radius  x   y   z   vx   vy   vz   nx   ny   nz   omega\n')
    for p in self.particles:
      x, y, z = p.r
      vx, vy, vz = p.v
      nx, ny, nz = p.n
      out.write('%d  %d  %f %f  %f  %f  %f  %f  %f  %f  %f  %f   %f\n' % (p.idx,p.tp,p.R,x,y,z,vx,vy,vz,nx,ny,nz,p.omega))
    out.close()
    


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--a", type=float, default=8.0, help="parameter a")
parser.add_argument("-b", "--b", type=float, default=10.0, help="parameter b")
parser.add_argument("-f", "--phi",  type=float, default=0.5, help="packing fraction")
parser.add_argument("-o", "--output", type=str, default='out.dat', help="output file")
parser.add_argument("-v", "--vavr", type=float, default=1.0, help="average velocity")
args = parser.parse_args()

print
print "\tActive Particles on Curved Spaces (APCS)"
print "\tBuilding of a random configuration on a peanut"
print 
print "\tRastko Sknepnek"
print "\tUniversity of Dundee"
print "\t(c) 2013, 2014"
print "\t----------------------------------------------"
print 
print "\ta : ", args.a
print "\tb : ", args.b
print "\tPacking fraction : ", args.phi
print "\tAverage velocity : ", args.vavr
###### WRONG !!!!#########
N=int(round(4.0*args.a**2*args.phi))   
##########################
print "\tTotal number of particles : ", N
print "\tOutput file : ", args.output
print 

start = datetime.now()


p = Peanut(args.a, args.b, N, args.vavr)
p.write(args.output)

end = datetime.now()

total = end - start

print 
print "  *** Completed in ", total.total_seconds(), " seconds *** "
print
