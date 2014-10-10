/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_morse_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Oct-2014
 * \brief Implementation of PairMorsePotential class
 */ 

#include "pair_morse_potential.hpp"

void PairMorsePotential::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double D = m_D;
  double a = m_a;
  double re = m_re;
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
 
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("morse",0.0);
    }
  }

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    double ai = pi.get_radius();
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      if (periodic)
      {
        if (dx > box->xhi) dx -= box->Lx;
        else if (dx < box->xlo) dx += box->Lx;
        if (dy > box->yhi) dy -= box->Ly;
        else if (dy < box->ylo) dy += box->Ly;
        if (dz > box->zhi) dz -= box->Lz;
        else if (dz < box->zlo) dz += box->Lz;
      }
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= rcut_sq)
      {
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          D = m_pair_params[pi_t][pj_t].D;
          a = m_pair_params[pi_t][pj_t].a;
          re = m_pair_params[pi_t][pj_t].re;
        }
        if (m_use_particle_radii)
          re = ai + pj.get_radius();
        double r = sqrt(r_sq);
        double exp_fact = exp(-a*(r-re));
        double pot_fact = exp_fact - 1.0;
        // Handle potential 
        double potential_energy = D*(pot_fact*pot_fact-1.0);
        if (m_shifted)
        {
          double shift_fact = exp(-a*(rcut-re)) - 1.0;
          potential_energy -= D*(shift_fact*shift_fact-1.0);
        }
        m_potential_energy += potential_energy;
        // Handle force
        double force_factor = 2.0*D*a*exp_fact*pot_fact/r;
        pi.fx += force_factor*dx;
        pi.fy += force_factor*dy;
        pi.fz += force_factor*dz;
        // Use 3d Newton's law
        pj.fx -= force_factor*dx;
        pj.fy -= force_factor*dy;
        pj.fz -= force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("morse",potential_energy);
          pj.add_pot_energy("morse",potential_energy);
        }
      }
    }
  }
}