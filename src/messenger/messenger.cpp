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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file messenger.coo
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of Messenger class
 */ 

#include "messenger.hpp"

//! Construct Messenger object
//! \param name file name to which to send messages or "terminal" to send it to terminal
Messenger::Messenger(string& name) : m_file_name(name), m_out( (to_lower_copy(name) == "terminal")  ? cout : *(new ofstream(name.c_str(),ios_base::out)) )
{ 
  m_msg_type[Messenger::INFO] = "INFO: ";
  m_msg_type[Messenger::WARNING] = "WARNING: ";
  m_msg_type[Messenger::ERROR] = "ERROR: ";
}

//! Destructor
Messenger::~Messenger() 
{
  if (to_lower_copy(m_file_name) != "terminal")
    delete (&m_out);
}

//! output message
//! Message format: TYPE: DATE TIME: MESSAGE
//! \param type Type of the message (INFO, WARNING, or ERROR)
//! \param m The message to print
void Messenger::msg(const Messenger::MSG_TYPE& type,const string& m) const
{
  ptime now = second_clock::local_time();
  string val = m_msg_type.find(type)->second; // Workaround the fact that std::map::operator[] cannot pass constant arguments
  m_out << val << now << ": " << m << endl;
}