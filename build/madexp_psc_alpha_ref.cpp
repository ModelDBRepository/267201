/*
 *  madexp_psc_alpha_ref.cpp
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
 *  Generated from NESTML at time: 2021-09-17 08:49:05.952704
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

#include "madexp_psc_alpha_ref.h"

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<madexp_psc_alpha_ref> madexp_psc_alpha_ref::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<madexp_psc_alpha_ref>::create()
  {
    // add state variables to recordables map
    insert_("epsilon", &madexp_psc_alpha_ref::get_epsilon);
    insert_("V_m", &madexp_psc_alpha_ref::get_V_m);
    insert_("w", &madexp_psc_alpha_ref::get_w);
    insert_("I_syn_ex__X__spikesExc", &madexp_psc_alpha_ref::get_I_syn_ex__X__spikesExc);
    insert_("I_syn_ex__X__spikesExc__d", &madexp_psc_alpha_ref::get_I_syn_ex__X__spikesExc__d);
    insert_("I_syn_in__X__spikesInh", &madexp_psc_alpha_ref::get_I_syn_in__X__spikesInh);
    insert_("I_syn_in__X__spikesInh__d", &madexp_psc_alpha_ref::get_I_syn_in__X__spikesInh__d);
  }
}

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of madexp_psc_alpha_ref's constructor.
// ---------------------------------------------------------------------------

madexp_psc_alpha_ref::Parameters_::Parameters_()
{
}

madexp_psc_alpha_ref::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

madexp_psc_alpha_ref::Buffers_::Buffers_(madexp_psc_alpha_ref &n):
  logger_(n)
  , __s( 0 ), __c( 0 ), __e( 0 )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

madexp_psc_alpha_ref::Buffers_::Buffers_(const Buffers_ &, madexp_psc_alpha_ref &n):
  logger_(n)
  , __s( 0 ), __c( 0 ), __e( 0 )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

madexp_psc_alpha_ref::madexp_psc_alpha_ref():ArchivingNode(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();

  calibrate();

  // use a default "good enough" value for the absolute error. It can be adjusted via `node.set()`
  P_.__gsl_error_tol = 1e-3;

  // initial values for parameters
  P_.C_m = 100.0; // as pF
  P_.g_L = 9.0; // as nS
  P_.Delta_T = 1.0; // as mV
  P_.V_peak = 0; // as mV
  P_.E_0 = (-65.0); // as mV
  P_.E_u = (-58.0); // as mV
  P_.E_d = (-50.0); // as mV
  P_.E_f = (-60.0); // as mV
  P_.epsilon_0 = 0.5; // as real
  P_.epsilon_c = 0.2; // as real
  P_.alpha = 1.0; // as real
  P_.delta = 0.1; // as real
  P_.tau_e = 1000.0; // as ms
  P_.I_e = 0.0; // as pA
  P_.V_th = (-55.0); // as mV
  P_.a = 1.0; // as nS
  P_.b = 2.0; // as pA
  P_.gamma = 1000.0; // as pA
  P_.I_KATP = 10.0; // as pA
  P_.tau_w = 300.0; // as ms
  P_.V_reset = (-62.0); // as mV
  P_.tau_syn_ex = 0.2; // as ms
  P_.tau_syn_in = 2.0; // as ms
  P_.t_ref = 2.0; // as ms

  // initial values for state variables
  S_.ode_state[State_::r] = 0; // as integer
  S_.ode_state[State_::epsilon] = P_.alpha * P_.epsilon_0; // as real
  S_.ode_state[State_::V_m] = P_.E_0 + ((-P_.E_0) + P_.E_u) * ((-S_.ode_state[State_::epsilon]) / P_.epsilon_0 + 1); // as mV
  S_.ode_state[State_::w] = 0; // as pA
  S_.ode_state[State_::I_syn_ex__X__spikesExc] = 0; // as real
  S_.ode_state[State_::I_syn_ex__X__spikesExc__d] = 0; // as real
  S_.ode_state[State_::I_syn_in__X__spikesInh] = 0; // as real
  S_.ode_state[State_::I_syn_in__X__spikesInh__d] = 0; // as real
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

madexp_psc_alpha_ref::madexp_psc_alpha_ref(const madexp_psc_alpha_ref& __n):
  ArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this) {
  // copy parameter struct P_
  P_.C_m = __n.P_.C_m;
  P_.g_L = __n.P_.g_L;
  P_.Delta_T = __n.P_.Delta_T;
  P_.V_peak = __n.P_.V_peak;
  P_.E_0 = __n.P_.E_0;
  P_.E_u = __n.P_.E_u;
  P_.E_d = __n.P_.E_d;
  P_.E_f = __n.P_.E_f;
  P_.epsilon_0 = __n.P_.epsilon_0;
  P_.epsilon_c = __n.P_.epsilon_c;
  P_.alpha = __n.P_.alpha;
  P_.delta = __n.P_.delta;
  P_.tau_e = __n.P_.tau_e;
  P_.I_e = __n.P_.I_e;
  P_.V_th = __n.P_.V_th;
  P_.a = __n.P_.a;
  P_.b = __n.P_.b;
  P_.gamma = __n.P_.gamma;
  P_.I_KATP = __n.P_.I_KATP;
  P_.tau_w = __n.P_.tau_w;
  P_.V_reset = __n.P_.V_reset;
  P_.tau_syn_ex = __n.P_.tau_syn_ex;
  P_.tau_syn_in = __n.P_.tau_syn_in;
  P_.t_ref = __n.P_.t_ref;

  // copy state struct S_
  S_.ode_state[State_::r] = __n.S_.ode_state[State_::r];
  S_.ode_state[State_::epsilon] = __n.S_.ode_state[State_::epsilon];
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::w] = __n.S_.ode_state[State_::w];
  S_.ode_state[State_::I_syn_ex__X__spikesExc] = __n.S_.ode_state[State_::I_syn_ex__X__spikesExc];
  S_.ode_state[State_::I_syn_ex__X__spikesExc__d] = __n.S_.ode_state[State_::I_syn_ex__X__spikesExc__d];
  S_.ode_state[State_::I_syn_in__X__spikesInh] = __n.S_.ode_state[State_::I_syn_in__X__spikesInh];
  S_.ode_state[State_::I_syn_in__X__spikesInh__d] = __n.S_.ode_state[State_::I_syn_in__X__spikesInh__d];
  S_.ode_state[State_::r] = __n.S_.ode_state[State_::r];
  S_.ode_state[State_::epsilon] = __n.S_.ode_state[State_::epsilon];
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::w] = __n.S_.ode_state[State_::w];
  S_.ode_state[State_::I_syn_ex__X__spikesExc] = __n.S_.ode_state[State_::I_syn_ex__X__spikesExc];
  S_.ode_state[State_::I_syn_ex__X__spikesExc__d] = __n.S_.ode_state[State_::I_syn_ex__X__spikesExc__d];
  S_.ode_state[State_::I_syn_in__X__spikesInh] = __n.S_.ode_state[State_::I_syn_in__X__spikesInh];
  S_.ode_state[State_::I_syn_in__X__spikesInh__d] = __n.S_.ode_state[State_::I_syn_in__X__spikesInh__d];

  // copy internals V_
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__h = __n.V_.__h;
  V_.Vspike = __n.V_.Vspike;
  V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc = __n.V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc;
  V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d = __n.V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d;
  V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc = __n.V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc;
  V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d = __n.V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d;
  V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh = __n.V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh;
  V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d = __n.V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d;
  V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh = __n.V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh;
  V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d = __n.V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

madexp_psc_alpha_ref::~madexp_psc_alpha_ref()
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

void madexp_psc_alpha_ref::init_buffers_()
{
  get_spikesInh().clear(); //includes resize
  get_spikesExc().clear(); //includes resize
  get_currents().clear(); //includes resize
  B_.logger_.reset(); // includes resize
  ArchivingNode::clear_history();

  if ( B_.__s == 0 )
  {
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 8 );
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
    B_.__e = gsl_odeiv_evolve_alloc( 8 );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = madexp_psc_alpha_ref_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 8;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void madexp_psc_alpha_ref::calibrate()
{
  B_.logger_.init();

  // internals V_
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps();
  V_.__h =nest::Time::get_resolution().get_ms();
  V_.Vspike =(P_.Delta_T==0) ? (std::min(P_.V_th, P_.V_peak)) : (P_.V_peak);
  V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc =1.0 * (V_.__h + P_.tau_syn_ex) * std::exp((-V_.__h) / P_.tau_syn_ex) / P_.tau_syn_ex;
  V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d =1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_ex);
  V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc =(-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_ex) / pow(P_.tau_syn_ex, 2);
  V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d =1.0 * ((-V_.__h) + P_.tau_syn_ex) * std::exp((-V_.__h) / P_.tau_syn_ex) / P_.tau_syn_ex;
  V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh =1.0 * (V_.__h + P_.tau_syn_in) * std::exp((-V_.__h) / P_.tau_syn_in) / P_.tau_syn_in;
  V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d =1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_in);
  V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh =(-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_in) / pow(P_.tau_syn_in, 2);
  V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d =1.0 * ((-V_.__h) + P_.tau_syn_in) * std::exp((-V_.__h) / P_.tau_syn_in) / P_.tau_syn_in;

  // state S_

  // buffers B_
}

// ---------------------------------------------------------------------------
//   Functions defined in the NESTML model
// ---------------------------------------------------------------------------
//
double madexp_psc_alpha_ref::I_spike(double epsilon, double V_m) const

{  
    double Ispk = 0.0;
    double arg_exp = 0.0;
    if (P_.Delta_T>0.0)
    {
      arg_exp = std::min(20.0, std::max((-20.0), (V_m - P_.V_th) / P_.Delta_T));
      Ispk = (epsilon - P_.epsilon_c) * P_.g_L * P_.Delta_T * std::exp(arg_exp) / P_.epsilon_0;
    }
    return Ispk;
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int madexp_psc_alpha_ref_dynamics(double, const double ode_state[], double f[], void* pnode)
{
  typedef madexp_psc_alpha_ref::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const madexp_psc_alpha_ref& node = *( reinterpret_cast< madexp_psc_alpha_ref* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = node.get_g_L() * (node.get_E_0() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_C_m() * node.get_epsilon_0()) + node.get_E_u() / node.get_C_m() - node.get_E_u() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_C_m() * node.get_epsilon_0()) - std::min(ode_state[State_::V_m], node.get_V_peak()) / node.get_C_m()) + (node.get_I_e() + ode_state[State_::I_syn_ex__X__spikesExc] + ode_state[State_::I_syn_in__X__spikesInh] + node.B_.currents_grid_sum_ - ode_state[State_::w] + node.I_spike(std::max(0.0, ode_state[State_::epsilon]), std::min(ode_state[State_::V_m], node.get_V_peak()))) / node.get_C_m();
  f[State_::w] = node.get_E_0() * (-(node.get_a()) * node.get_epsilon_c() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_epsilon_0() * (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]))) - node.get_a() * pow(std::max(0.0, ode_state[State_::epsilon]), 2) / (node.get_epsilon_0() * (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon])))) + node.get_a() * (-(node.get_E_u()) * node.get_epsilon_c() / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon])) - node.get_E_u() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon])) + node.get_E_u() * node.get_epsilon_c() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_epsilon_0() * (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]))) + node.get_E_u() * pow(std::max(0.0, ode_state[State_::epsilon]), 2) / (node.get_epsilon_0() * (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]))) + node.get_epsilon_c() * std::min(ode_state[State_::V_m], node.get_V_peak()) / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon])) + std::max(0.0, ode_state[State_::epsilon]) * std::min(ode_state[State_::V_m], node.get_V_peak()) / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]))) + node.get_epsilon_c() * (node.get_I_KATP() / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon])) - ode_state[State_::w] / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]))) - ode_state[State_::w] * std::max(0.0, ode_state[State_::epsilon]) / (node.get_epsilon_c() * node.get_tau_w() + node.get_tau_w() * std::max(0.0, ode_state[State_::epsilon]));
  f[State_::epsilon] = pow(node.get_epsilon_0(), 3) * (node.get_E_d() * node.get_gamma() / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) - node.get_E_d() * ode_state[State_::w] / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) + node.get_E_f() * ode_state[State_::w] / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) - node.get_gamma() * std::min(ode_state[State_::V_m], node.get_V_peak()) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e())) + (-(3) * node.get_E_d() * pow(node.get_epsilon_0(), 2) * node.get_gamma() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) + 3 * node.get_E_f() * pow(node.get_epsilon_0(), 2) * node.get_gamma() * std::max(0.0, ode_state[State_::epsilon]) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e())) / node.get_alpha() + (3 * node.get_E_d() * node.get_epsilon_0() * node.get_gamma() * pow(std::max(0.0, ode_state[State_::epsilon]), 2) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) - 3 * node.get_E_f() * node.get_epsilon_0() * node.get_gamma() * pow(std::max(0.0, ode_state[State_::epsilon]), 2) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e())) / pow(node.get_alpha(), 2) + (-(node.get_E_d()) * node.get_gamma() * pow(std::max(0.0, ode_state[State_::epsilon]), 3) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e()) + node.get_E_f() * node.get_gamma() * pow(std::max(0.0, ode_state[State_::epsilon]), 3) / (node.get_E_d() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e() - node.get_E_f() * pow(node.get_epsilon_0(), 3) * node.get_gamma() * node.get_tau_e())) / pow(node.get_alpha(), 3);
  f[State_::I_syn_ex__X__spikesExc] = 1.0 * ode_state[State_::I_syn_ex__X__spikesExc__d];
  f[State_::I_syn_ex__X__spikesExc__d] = -(ode_state[State_::I_syn_ex__X__spikesExc]) / pow(node.get_tau_syn_ex(), 2) - 2 * ode_state[State_::I_syn_ex__X__spikesExc__d] / node.get_tau_syn_ex();
  f[State_::I_syn_in__X__spikesInh] = 1.0 * ode_state[State_::I_syn_in__X__spikesInh__d];
  f[State_::I_syn_in__X__spikesInh__d] = -(ode_state[State_::I_syn_in__X__spikesInh]) / pow(node.get_tau_syn_in(), 2) - 2 * ode_state[State_::I_syn_in__X__spikesInh__d] / node.get_tau_syn_in();
  f[State_::r] = 0.;

  return GSL_SUCCESS;
}

void madexp_psc_alpha_ref::update(nest::Time const & origin,const long from, const long to)
{
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag )
  {
    B_.spikesInh_grid_sum_ = get_spikesInh().get_value(lag);
    B_.spikesExc_grid_sum_ = get_spikesExc().get_value(lag);
    B_.currents_grid_sum_ = get_currents().get_value(lag);

    // NESTML generated code for the update block:
      double I_syn_ex__X__spikesExc__tmp = S_.ode_state[State_::I_syn_ex__X__spikesExc] * V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc + S_.ode_state[State_::I_syn_ex__X__spikesExc__d] * V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d;
      double I_syn_ex__X__spikesExc__d__tmp = S_.ode_state[State_::I_syn_ex__X__spikesExc] * V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc + S_.ode_state[State_::I_syn_ex__X__spikesExc__d] * V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d;
      double I_syn_in__X__spikesInh__tmp = S_.ode_state[State_::I_syn_in__X__spikesInh] * V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh + S_.ode_state[State_::I_syn_in__X__spikesInh__d] * V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d;
      double I_syn_in__X__spikesInh__d__tmp = S_.ode_state[State_::I_syn_in__X__spikesInh] * V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh + S_.ode_state[State_::I_syn_in__X__spikesInh__d] * V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d;
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
    S_.ode_state[State_::I_syn_ex__X__spikesExc] = I_syn_ex__X__spikesExc__tmp;
    S_.ode_state[State_::I_syn_ex__X__spikesExc__d] = I_syn_ex__X__spikesExc__d__tmp;
    S_.ode_state[State_::I_syn_in__X__spikesInh] = I_syn_in__X__spikesInh__tmp;
    S_.ode_state[State_::I_syn_in__X__spikesInh__d] = I_syn_in__X__spikesInh__d__tmp;
    S_.ode_state[State_::I_syn_ex__X__spikesExc__d] += (B_.spikesExc_grid_sum_) * (numerics::e / P_.tau_syn_ex) / (1.0);
    S_.ode_state[State_::I_syn_in__X__spikesInh__d] += (B_.spikesInh_grid_sum_) * (numerics::e / P_.tau_syn_in) / (1.0);
    if (S_.ode_state[State_::epsilon]<0.0)
    {
      S_.ode_state[State_::epsilon] = 0.0;
    }
    if (S_.ode_state[State_::r]>0)
    {
      S_.ode_state[State_::r] -= 1;
      S_.ode_state[State_::V_m] = P_.V_reset;
    }
    else if (S_.ode_state[State_::V_m]>=V_.Vspike&&S_.ode_state[State_::epsilon]>P_.epsilon_c)
    {
      S_.ode_state[State_::V_m] = P_.V_reset;
      S_.ode_state[State_::epsilon] -= P_.delta;
      S_.ode_state[State_::w] += P_.b;
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
      S_.ode_state[State_::r] = V_.RefractoryCounts;
    }

    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void madexp_psc_alpha_ref::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

void madexp_psc_alpha_ref::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  if ( weight < 0.0 )
  {
    // inhibitory
    get_spikesInh().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                      
                       weight * multiplicity );
  }
  if ( weight >= 0.0 )
  {
    // excitatory
    get_spikesExc().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
  }
}

void madexp_psc_alpha_ref::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_currents().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}
