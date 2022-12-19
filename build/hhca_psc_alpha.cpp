/*
 *  hhca_psc_alpha.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Generated from NESTML at time: 2021-09-17 08:49:06.755693
**/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "hhca_psc_alpha.h"

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<hhca_psc_alpha> hhca_psc_alpha::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<hhca_psc_alpha>::create()
  {
    // add state variables to recordables map
    insert_("V_m", &hhca_psc_alpha::get_V_m);
    insert_("Ca_i", &hhca_psc_alpha::get_Ca_i);
    insert_("alpha_n_init", &hhca_psc_alpha::get_alpha_n_init);
    insert_("beta_n_init", &hhca_psc_alpha::get_beta_n_init);
    insert_("alpha_m_init", &hhca_psc_alpha::get_alpha_m_init);
    insert_("beta_m_init", &hhca_psc_alpha::get_beta_m_init);
    insert_("alpha_h_init", &hhca_psc_alpha::get_alpha_h_init);
    insert_("beta_h_init", &hhca_psc_alpha::get_beta_h_init);
    insert_("Act_m", &hhca_psc_alpha::get_Act_m);
    insert_("Inact_h", &hhca_psc_alpha::get_Inact_h);
    insert_("Act_n", &hhca_psc_alpha::get_Act_n);
    insert_("h_Ca", &hhca_psc_alpha::get_h_Ca);
    insert_("m_Ca", &hhca_psc_alpha::get_m_Ca);
    insert_("m_AHP", &hhca_psc_alpha::get_m_AHP);
    insert_("I_syn_ex__X__spikeExc", &hhca_psc_alpha::get_I_syn_ex__X__spikeExc);
    insert_("I_syn_ex__X__spikeExc__d", &hhca_psc_alpha::get_I_syn_ex__X__spikeExc__d);
    insert_("I_syn_in__X__spikeInh", &hhca_psc_alpha::get_I_syn_in__X__spikeInh);
    insert_("I_syn_in__X__spikeInh__d", &hhca_psc_alpha::get_I_syn_in__X__spikeInh__d);
  }
}

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of hhca_psc_alpha's constructor.
// ---------------------------------------------------------------------------

hhca_psc_alpha::Parameters_::Parameters_()
{
}

hhca_psc_alpha::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

hhca_psc_alpha::Buffers_::Buffers_(hhca_psc_alpha &n):
  logger_(n)
  , __s( 0 ), __c( 0 ), __e( 0 )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

hhca_psc_alpha::Buffers_::Buffers_(const Buffers_ &, hhca_psc_alpha &n):
  logger_(n)
  , __s( 0 ), __c( 0 ), __e( 0 )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

hhca_psc_alpha::hhca_psc_alpha():ArchivingNode(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();

  calibrate();

  // use a default "good enough" value for the absolute error. It can be adjusted via `node.set()`
  P_.__gsl_error_tol = 1e-3;

  // initial values for parameters
  P_.t_ref = 2.0; // as ms
  P_.g_Na = 12000.0; // as nS
  P_.g_K = 3600.0; // as nS
  P_.g_L = 30; // as nS
  P_.C_m = 100.0; // as pF
  P_.E_Na = 50; // as mV
  P_.E_K = (-77.0); // as mV
  P_.E_L = (-54.402); // as mV
  P_.tau_syn_ex = 0.2; // as ms
  P_.tau_syn_in = 2.0; // as ms
  P_.g_Ca = 0.001; // as nS
  P_.Ca_env = 1.0 / 1.0; // as mmol / l
  P_.Ca_i_eq = 5e-05 / 1.0; // as mmol / l
  P_.tau_Ca = 500.0; // as ms
  P_.V_hmCa = (-27.5); // as mV
  P_.k_mCa = 5.7; // as mV
  P_.tau_mCa = 0.5; // as ms
  P_.V_hhCa = (-52.4); // as mV
  P_.k_hCa = 5.2; // as mV
  P_.tau_hCa = 18.0; // as ms
  P_.g_AHP = 0.2; // as nS
  P_.K_AHP = 0.002 / 1.0; // as mmol / l
  P_.b_AHP = 2.5; // as real
  P_.tau_AHP = 30.0; // as ms
  P_.E_0 = 26.64; // as mV
  P_.Vcell = 4 / 3.0 * 3.14159 * pow((5 * 1.0), 3); // as um3
  P_.I_e = 0; // as pA

  // initial values for state variables
  S_.ode_state[State_::r] = 0; // as integer
  S_.ode_state[State_::V_m] = (-65.0); // as mV
  S_.ode_state[State_::Ca_i] = P_.Ca_i_eq; // as mmol / l
  S_.ode_state[State_::alpha_n_init] = (0.01 * (S_.ode_state[State_::V_m] / 1.0 + 55.0)) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 55.0)) / 10.0)); // as real
  S_.ode_state[State_::beta_n_init] = 0.125 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 80.0); // as real
  S_.ode_state[State_::alpha_m_init] = (0.1 * (S_.ode_state[State_::V_m] / 1.0 + 40.0)) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 40.0)) / 10.0)); // as real
  S_.ode_state[State_::beta_m_init] = 4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0); // as real
  S_.ode_state[State_::alpha_h_init] = 0.07 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 20.0); // as real
  S_.ode_state[State_::beta_h_init] = 1.0 / (1.0 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 35.0)) / 10.0)); // as real
  S_.ode_state[State_::Act_m] = S_.ode_state[State_::alpha_m_init] / (S_.ode_state[State_::alpha_m_init] + S_.ode_state[State_::beta_m_init]); // as real
  S_.ode_state[State_::Inact_h] = S_.ode_state[State_::alpha_h_init] / (S_.ode_state[State_::alpha_h_init] + S_.ode_state[State_::beta_h_init]); // as real
  S_.ode_state[State_::Act_n] = S_.ode_state[State_::alpha_n_init] / (S_.ode_state[State_::alpha_n_init] + S_.ode_state[State_::beta_n_init]); // as real
  S_.ode_state[State_::h_Ca] = y_shape(S_.ode_state[State_::V_m], P_.V_hhCa, (-P_.k_hCa)); // as real
  S_.ode_state[State_::m_Ca] = y_shape(S_.ode_state[State_::V_m], P_.V_hmCa, P_.k_mCa); // as real
  S_.ode_state[State_::m_AHP] = 0.0; // as real
  S_.ode_state[State_::I_syn_ex__X__spikeExc] = 0; // as real
  S_.ode_state[State_::I_syn_ex__X__spikeExc__d] = 0; // as real
  S_.ode_state[State_::I_syn_in__X__spikeInh] = 0; // as real
  S_.ode_state[State_::I_syn_in__X__spikeInh__d] = 0; // as real
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

hhca_psc_alpha::hhca_psc_alpha(const hhca_psc_alpha& __n):
  ArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this) {
  // copy parameter struct P_
  P_.t_ref = __n.P_.t_ref;
  P_.g_Na = __n.P_.g_Na;
  P_.g_K = __n.P_.g_K;
  P_.g_L = __n.P_.g_L;
  P_.C_m = __n.P_.C_m;
  P_.E_Na = __n.P_.E_Na;
  P_.E_K = __n.P_.E_K;
  P_.E_L = __n.P_.E_L;
  P_.tau_syn_ex = __n.P_.tau_syn_ex;
  P_.tau_syn_in = __n.P_.tau_syn_in;
  P_.g_Ca = __n.P_.g_Ca;
  P_.Ca_env = __n.P_.Ca_env;
  P_.Ca_i_eq = __n.P_.Ca_i_eq;
  P_.tau_Ca = __n.P_.tau_Ca;
  P_.V_hmCa = __n.P_.V_hmCa;
  P_.k_mCa = __n.P_.k_mCa;
  P_.tau_mCa = __n.P_.tau_mCa;
  P_.V_hhCa = __n.P_.V_hhCa;
  P_.k_hCa = __n.P_.k_hCa;
  P_.tau_hCa = __n.P_.tau_hCa;
  P_.g_AHP = __n.P_.g_AHP;
  P_.K_AHP = __n.P_.K_AHP;
  P_.b_AHP = __n.P_.b_AHP;
  P_.tau_AHP = __n.P_.tau_AHP;
  P_.E_0 = __n.P_.E_0;
  P_.Vcell = __n.P_.Vcell;
  P_.I_e = __n.P_.I_e;

  // copy state struct S_
  S_.ode_state[State_::r] = __n.S_.ode_state[State_::r];
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::Ca_i] = __n.S_.ode_state[State_::Ca_i];
  S_.ode_state[State_::alpha_n_init] = __n.S_.ode_state[State_::alpha_n_init];
  S_.ode_state[State_::beta_n_init] = __n.S_.ode_state[State_::beta_n_init];
  S_.ode_state[State_::alpha_m_init] = __n.S_.ode_state[State_::alpha_m_init];
  S_.ode_state[State_::beta_m_init] = __n.S_.ode_state[State_::beta_m_init];
  S_.ode_state[State_::alpha_h_init] = __n.S_.ode_state[State_::alpha_h_init];
  S_.ode_state[State_::beta_h_init] = __n.S_.ode_state[State_::beta_h_init];
  S_.ode_state[State_::Act_m] = __n.S_.ode_state[State_::Act_m];
  S_.ode_state[State_::Inact_h] = __n.S_.ode_state[State_::Inact_h];
  S_.ode_state[State_::Act_n] = __n.S_.ode_state[State_::Act_n];
  S_.ode_state[State_::h_Ca] = __n.S_.ode_state[State_::h_Ca];
  S_.ode_state[State_::m_Ca] = __n.S_.ode_state[State_::m_Ca];
  S_.ode_state[State_::m_AHP] = __n.S_.ode_state[State_::m_AHP];
  S_.ode_state[State_::I_syn_ex__X__spikeExc] = __n.S_.ode_state[State_::I_syn_ex__X__spikeExc];
  S_.ode_state[State_::I_syn_ex__X__spikeExc__d] = __n.S_.ode_state[State_::I_syn_ex__X__spikeExc__d];
  S_.ode_state[State_::I_syn_in__X__spikeInh] = __n.S_.ode_state[State_::I_syn_in__X__spikeInh];
  S_.ode_state[State_::I_syn_in__X__spikeInh__d] = __n.S_.ode_state[State_::I_syn_in__X__spikeInh__d];
  S_.ode_state[State_::r] = __n.S_.ode_state[State_::r];
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::Ca_i] = __n.S_.ode_state[State_::Ca_i];
  S_.ode_state[State_::alpha_n_init] = __n.S_.ode_state[State_::alpha_n_init];
  S_.ode_state[State_::beta_n_init] = __n.S_.ode_state[State_::beta_n_init];
  S_.ode_state[State_::alpha_m_init] = __n.S_.ode_state[State_::alpha_m_init];
  S_.ode_state[State_::beta_m_init] = __n.S_.ode_state[State_::beta_m_init];
  S_.ode_state[State_::alpha_h_init] = __n.S_.ode_state[State_::alpha_h_init];
  S_.ode_state[State_::beta_h_init] = __n.S_.ode_state[State_::beta_h_init];
  S_.ode_state[State_::Act_m] = __n.S_.ode_state[State_::Act_m];
  S_.ode_state[State_::Inact_h] = __n.S_.ode_state[State_::Inact_h];
  S_.ode_state[State_::Act_n] = __n.S_.ode_state[State_::Act_n];
  S_.ode_state[State_::h_Ca] = __n.S_.ode_state[State_::h_Ca];
  S_.ode_state[State_::m_Ca] = __n.S_.ode_state[State_::m_Ca];
  S_.ode_state[State_::m_AHP] = __n.S_.ode_state[State_::m_AHP];
  S_.ode_state[State_::I_syn_ex__X__spikeExc] = __n.S_.ode_state[State_::I_syn_ex__X__spikeExc];
  S_.ode_state[State_::I_syn_ex__X__spikeExc__d] = __n.S_.ode_state[State_::I_syn_ex__X__spikeExc__d];
  S_.ode_state[State_::I_syn_in__X__spikeInh] = __n.S_.ode_state[State_::I_syn_in__X__spikeInh];
  S_.ode_state[State_::I_syn_in__X__spikeInh__d] = __n.S_.ode_state[State_::I_syn_in__X__spikeInh__d];

  // copy internals V_
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__h = __n.V_.__h;
  V_.charge = __n.V_.charge;
  V_.avogadro = __n.V_.avogadro;
  V_.c2c = __n.V_.c2c;
  V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc = __n.V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc;
  V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d = __n.V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d;
  V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc = __n.V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc;
  V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d = __n.V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d;
  V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh = __n.V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh;
  V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d = __n.V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d;
  V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh = __n.V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh;
  V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d = __n.V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

hhca_psc_alpha::~hhca_psc_alpha()
{
  // GSL structs may not have been allocated, so we need to protect destruction

  if (B_.__s)
  {
    gsl_odeiv_step_free( B_.__s );
  }

  if (B_.__c)
  {
    gsl_odeiv_control_free( B_.__c );
  }

  if (B_.__e)
  {
    gsl_odeiv_evolve_free( B_.__e );
  }
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------

void hhca_psc_alpha::init_buffers_()
{
  get_spikeInh().clear(); //includes resize
  get_spikeExc().clear(); //includes resize
  get_I_stim().clear(); //includes resize
  B_.logger_.reset(); // includes resize
  ArchivingNode::clear_history();

  if ( B_.__s == 0 )
  {
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 19 );
  }
  else
  {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( B_.__c == 0 )
  {
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.__e == 0 )
  {
    B_.__e = gsl_odeiv_evolve_alloc( 19 );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = hhca_psc_alpha_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 19;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void hhca_psc_alpha::calibrate()
{
  B_.logger_.init();

  // internals V_
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps();
  V_.__h =nest::Time::get_resolution().get_ms();
  V_.charge =1.6e-19;
  V_.avogadro =6.02e+23 / 1.0;
  V_.c2c =1000.0000000000005 * (1.0 / (V_.charge * V_.avogadro * P_.Vcell));
  V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc =1.0 * (V_.__h + P_.tau_syn_ex) * std::exp((-V_.__h) / P_.tau_syn_ex) / P_.tau_syn_ex;
  V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d =1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_ex);
  V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc =(-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_ex) / pow(P_.tau_syn_ex, 2);
  V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d =1.0 * ((-V_.__h) + P_.tau_syn_ex) * std::exp((-V_.__h) / P_.tau_syn_ex) / P_.tau_syn_ex;
  V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh =1.0 * (V_.__h + P_.tau_syn_in) * std::exp((-V_.__h) / P_.tau_syn_in) / P_.tau_syn_in;
  V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d =1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_in);
  V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh =(-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_in) / pow(P_.tau_syn_in, 2);
  V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d =1.0 * ((-V_.__h) + P_.tau_syn_in) * std::exp((-V_.__h) / P_.tau_syn_in) / P_.tau_syn_in;

  // state S_

  // buffers B_
}

// ---------------------------------------------------------------------------
//   Functions defined in the NESTML model
// ---------------------------------------------------------------------------
//
double hhca_psc_alpha::y_shape(double V_m, double Vh, double k_y) const

{  
    double arg_exp = std::min((-(V_m - Vh)) / k_y, 50.0);
    double res = std::min(1.0, std::max(0.0, 1.0 / (1.0 + std::exp(arg_exp))));
    return res;
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int hhca_psc_alpha_dynamics(double, const double ode_state[], double f[], void* pnode)
{
  typedef hhca_psc_alpha::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const hhca_psc_alpha& node = *( reinterpret_cast< hhca_psc_alpha* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = (pow(ode_state[State_::Act_m], 3) * node.get_E_Na() * ode_state[State_::Inact_h] * node.get_g_Na() - pow(ode_state[State_::Act_m], 3) * ode_state[State_::Inact_h] * ode_state[State_::V_m] * node.get_g_Na() + pow(ode_state[State_::Act_n], 4) * node.get_E_K() * node.get_g_K() - pow(ode_state[State_::Act_n], 4) * ode_state[State_::V_m] * node.get_g_K() + node.get_E_K() * node.get_g_AHP() * pow(ode_state[State_::m_AHP], 2) + node.get_E_L() * node.get_g_L() + node.get_I_e() + node.B_.I_stim_grid_sum_ + ode_state[State_::I_syn_ex__X__spikeExc] + ode_state[State_::I_syn_in__X__spikeInh] - ode_state[State_::V_m] * node.get_g_AHP() * pow(ode_state[State_::m_AHP], 2) - ode_state[State_::V_m] * node.get_g_L()) / node.get_C_m();
  f[State_::Ca_i] = -(ode_state[State_::Ca_i]) / node.get_tau_Ca() + node.get_Ca_i_eq() / node.get_tau_Ca() + 2 * node.get_E_0() * node.get_c2c() * node.get_g_Ca() * ode_state[State_::h_Ca] * ode_state[State_::m_Ca] * std::log(node.get_Ca_env() / ode_state[State_::Ca_i]) - ode_state[State_::V_m] * node.get_c2c() * node.get_g_Ca() * ode_state[State_::h_Ca] * ode_state[State_::m_Ca];
  f[State_::Act_n] = 0.01 * ode_state[State_::Act_n] * ode_state[State_::V_m] * std::exp(0.1125 * ode_state[State_::V_m]) / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m])) + 0.055468413760135 * ode_state[State_::Act_n] * std::exp(0.1 * ode_state[State_::V_m]) / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m])) + 0.55 * ode_state[State_::Act_n] * std::exp(0.1125 * ode_state[State_::V_m]) / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m])) - 0.000226686729091827 * ode_state[State_::Act_n] / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m])) - 0.01 * ode_state[State_::V_m] * std::exp(0.1125 * ode_state[State_::V_m]) / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m])) - 0.55 * std::exp(0.1125 * ode_state[State_::V_m]) / (0.00408677143846407 * std::exp(0.0125 * ode_state[State_::V_m]) - 1.0 * std::exp(0.1125 * ode_state[State_::V_m]));
  f[State_::Act_m] = 0.1 * ode_state[State_::Act_m] * ode_state[State_::V_m] * std::exp(0.155555555555556 * ode_state[State_::V_m]) / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m])) + 0.108087223804836 * ode_state[State_::Act_m] * std::exp(0.1 * ode_state[State_::V_m]) / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m])) + 4.0 * ode_state[State_::Act_m] * std::exp(0.155555555555556 * ode_state[State_::V_m]) / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m])) - 0.00197968655969517 * ode_state[State_::Act_m] / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m])) - 0.1 * ode_state[State_::V_m] * std::exp(0.155555555555556 * ode_state[State_::V_m]) / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m])) - 4.0 * std::exp(0.155555555555556 * ode_state[State_::V_m]) / (0.0183156388887342 * std::exp(0.0555555555555556 * ode_state[State_::V_m]) - 1.0 * std::exp(0.155555555555556 * ode_state[State_::V_m]));
  f[State_::Inact_h] = -(0.00271419454822054) * ode_state[State_::Inact_h] * std::exp(0.1 * ode_state[State_::V_m]) / (0.0301973834223185 * std::exp(0.05 * ode_state[State_::V_m]) + 1.0 * std::exp(0.15 * ode_state[State_::V_m])) - 1.0 * ode_state[State_::Inact_h] * std::exp(0.15 * ode_state[State_::V_m]) / (0.0301973834223185 * std::exp(0.05 * ode_state[State_::V_m]) + 1.0 * std::exp(0.15 * ode_state[State_::V_m])) - 8.19615734553822e-05 * ode_state[State_::Inact_h] / (0.0301973834223185 * std::exp(0.05 * ode_state[State_::V_m]) + 1.0 * std::exp(0.15 * ode_state[State_::V_m])) + 0.00271419454822054 * std::exp(0.1 * ode_state[State_::V_m]) / (0.0301973834223185 * std::exp(0.05 * ode_state[State_::V_m]) + 1.0 * std::exp(0.15 * ode_state[State_::V_m])) + 8.19615734553822e-05 / (0.0301973834223185 * std::exp(0.05 * ode_state[State_::V_m]) + 1.0 * std::exp(0.15 * ode_state[State_::V_m]));
  f[State_::h_Ca] = (-(ode_state[State_::h_Ca]) + node.y_shape(ode_state[State_::V_m], node.get_V_hhCa(), -(node.get_k_hCa()))) / node.get_tau_hCa();
  f[State_::m_Ca] = (-(ode_state[State_::m_Ca]) + node.y_shape(ode_state[State_::V_m], node.get_V_hmCa(), node.get_k_mCa())) / node.get_tau_mCa();
  f[State_::m_AHP] = -(pow(ode_state[State_::Ca_i], 2)) * ode_state[State_::m_AHP] / (pow(ode_state[State_::Ca_i], 2) * node.get_tau_AHP() + pow(node.get_K_AHP(), 2) * node.get_b_AHP() * node.get_tau_AHP()) + pow(ode_state[State_::Ca_i], 2) / (pow(ode_state[State_::Ca_i], 2) * node.get_tau_AHP() + pow(node.get_K_AHP(), 2) * node.get_b_AHP() * node.get_tau_AHP()) - pow(node.get_K_AHP(), 2) * node.get_b_AHP() * ode_state[State_::m_AHP] / (pow(ode_state[State_::Ca_i], 2) * node.get_tau_AHP() + pow(node.get_K_AHP(), 2) * node.get_b_AHP() * node.get_tau_AHP());
  f[State_::I_syn_ex__X__spikeExc] = 1.0 * ode_state[State_::I_syn_ex__X__spikeExc__d];
  f[State_::I_syn_ex__X__spikeExc__d] = -(ode_state[State_::I_syn_ex__X__spikeExc]) / pow(node.get_tau_syn_ex(), 2) - 2 * ode_state[State_::I_syn_ex__X__spikeExc__d] / node.get_tau_syn_ex();
  f[State_::I_syn_in__X__spikeInh] = 1.0 * ode_state[State_::I_syn_in__X__spikeInh__d];
  f[State_::I_syn_in__X__spikeInh__d] = -(ode_state[State_::I_syn_in__X__spikeInh]) / pow(node.get_tau_syn_in(), 2) - 2 * ode_state[State_::I_syn_in__X__spikeInh__d] / node.get_tau_syn_in();
  f[State_::r] = 0.;
  f[State_::alpha_n_init] = 0.;
  f[State_::beta_n_init] = 0.;
  f[State_::alpha_m_init] = 0.;
  f[State_::beta_m_init] = 0.;
  f[State_::alpha_h_init] = 0.;
  f[State_::beta_h_init] = 0.;

  return GSL_SUCCESS;
}

void hhca_psc_alpha::update(nest::Time const & origin,const long from, const long to)
{
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag )
  {
    B_.spikeInh_grid_sum_ = get_spikeInh().get_value(lag);
    B_.spikeExc_grid_sum_ = get_spikeExc().get_value(lag);
    B_.I_stim_grid_sum_ = get_I_stim().get_value(lag);

    // NESTML generated code for the update block:
    double U_old = S_.ode_state[State_::V_m];
      double I_syn_ex__X__spikeExc__tmp = S_.ode_state[State_::I_syn_ex__X__spikeExc] * V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc + S_.ode_state[State_::I_syn_ex__X__spikeExc__d] * V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d;
      double I_syn_ex__X__spikeExc__d__tmp = S_.ode_state[State_::I_syn_ex__X__spikeExc] * V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc + S_.ode_state[State_::I_syn_ex__X__spikeExc__d] * V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d;
      double I_syn_in__X__spikeInh__tmp = S_.ode_state[State_::I_syn_in__X__spikeInh] * V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh + S_.ode_state[State_::I_syn_in__X__spikeInh__d] * V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d;
      double I_syn_in__X__spikeInh__d__tmp = S_.ode_state[State_::I_syn_in__X__spikeInh] * V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh + S_.ode_state[State_::I_syn_in__X__spikeInh__d] * V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d;
    __t = 0;
    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( __t < B_.__step )
    {
      const int status = gsl_odeiv_evolve_apply(B_.__e,
                                                B_.__c,
                                                B_.__s,
                                                &B_.__sys,              // system of ODE
                                                &__t,                   // from t
                                                B_.__step,              // to t <= step
                                                &B_.__integration_step, // integration step size
                                                S_.ode_state);          // neuronal state

      if ( status != GSL_SUCCESS )
      {
        throw nest::GSLSolverFailure( get_name(), status );
      }
    }
    /* replace analytically solvable variables with precisely integrated values  */
    S_.ode_state[State_::I_syn_ex__X__spikeExc] = I_syn_ex__X__spikeExc__tmp;
    S_.ode_state[State_::I_syn_ex__X__spikeExc__d] = I_syn_ex__X__spikeExc__d__tmp;
    S_.ode_state[State_::I_syn_in__X__spikeInh] = I_syn_in__X__spikeInh__tmp;
    S_.ode_state[State_::I_syn_in__X__spikeInh__d] = I_syn_in__X__spikeInh__d__tmp;
    S_.ode_state[State_::I_syn_ex__X__spikeExc__d] += (B_.spikeExc_grid_sum_) * (numerics::e / P_.tau_syn_ex) / (1.0);
    S_.ode_state[State_::I_syn_in__X__spikeInh__d] += (B_.spikeInh_grid_sum_) * (numerics::e / P_.tau_syn_in) / (1.0);
    if (S_.ode_state[State_::r]>0)
    {
      S_.ode_state[State_::r] -= 1;
    }
    else if (S_.ode_state[State_::V_m]>0&&U_old>S_.ode_state[State_::V_m])
    {
      S_.ode_state[State_::r] = V_.RefractoryCounts;
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    }

    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void hhca_psc_alpha::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

void hhca_psc_alpha::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  if ( weight < 0.0 )
  {
    // inhibitory
    get_spikeInh().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                      
                       weight * multiplicity );
  }
  if ( weight >= 0.0 )
  {
    // excitatory
    get_spikeExc().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
  }
}

void hhca_psc_alpha::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}
