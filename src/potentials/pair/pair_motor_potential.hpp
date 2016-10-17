/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file pair_motor_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Sep-2016
 * \brief Declaration of PairMotorPotential class
 */ 

#ifndef __PAIR_MOTOR_POTENTIAL_HPP__
#define __PAIR_MOTOR_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the motor pair potential
struct MotorParameters
{
  double alpha;
  double a; 
};

/*! PairMotorPotential is a rather crude attempt to model effects of molecular motors 
 *  acting between two filaments (e.g., actin cables) in a particle based model. The idea here is 
 *  to intorduce an additional force that acts on a pair of particles if their directors (\f$ \vec n \f$), roughly
 *  speaking,  point in opposite directions. 
 *  More precisely, for a pair of particles \f$ i \f$ and \f$ j \f$ within distance \f$ a \f$, if \f$ \zeta = \vec n_i \cdot \vec n_j < 0 \f$,
 *  each particle receives a force of magnitude \f$ \alpha |\zeta| \f$. Force on each particle points in the direction of the corresponding director. 
 *  We note that this is quite artifical as it introduces motion of the centre of mass of the pair.  
 */
class PairMotorPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairMotorPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("alpha");
    m_known_params.push_back("a");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for motor pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in motor pair potential.");
    }
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No activiry (alpha) specified for soft pair potential. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global activity (alpha) for motor pair potential is set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.pair.motor.alpha",lexical_cast<string>(m_alpha));
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential range (a) specified for motor pair potential. Setting it to 2.");
      m_a = 2.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential range (a) for motor pair potential is set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
    }
    m_msg->write_config("potential.pair.motor.a",lexical_cast<string>(m_a));
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Motor pair potential is set to use particle radii to control its range. Parameter \"a\" will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.motor.use_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Motor pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.motor.phase_in","true");
    }    
    
    m_pair_params = new MotorParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new MotorParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].alpha = m_alpha;
        m_pair_params[i][j].a = m_a;
      }
    }
    
  }
  
  virtual ~PairMotorPotential()
  {
    for (int i = 0; i < m_ntypes; i++)
      delete [] m_pair_params[i];
    delete [] m_pair_params;
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in motor potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in motor potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Motor pair potential. Setting strength to "+pair_param["alpha"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Motor pair potential. Using default strength ("+lexical_cast<string>(m_alpha)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = m_alpha;
    }
    m_msg->write_config("potential.pair.soft.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".alpha",lexical_cast<string>(param["alpha"]));
    if (pair_param.find("a") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Motor pair potential. Setting range to "+pair_param["a"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = lexical_cast<double>(pair_param["a"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Motor pair potential. Using default range ("+lexical_cast<string>(m_a)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = m_a;
    }
    m_msg->write_config("potential.pair.soft.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".a",lexical_cast<string>(param["a"]));
        
    m_pair_params[type_1-1][type_2-1].alpha = param["alpha"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].alpha = param["alpha"];
    m_pair_params[type_1-1][type_2-1].a = param["a"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].a = param["a"];
    
    m_has_pair_params = true;
  }
  
  //! Returns true since soft potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_alpha;                   //!< activity
  double m_a;                       //!< potential range
  MotorParameters** m_pair_params;  //!< type specific pair parameters 
     
};

typedef shared_ptr<PairMotorPotential> PairMotorPotentialPtr;

#endif
