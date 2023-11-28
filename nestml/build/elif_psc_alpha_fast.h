/**
 *  elif_psc_alpha_fast.h
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
 *  Generated from NESTML at time: 2021-09-17 08:49:04.849510
**/
#ifndef ELIF_PSC_ALPHA_FAST
#define ELIF_PSC_ALPHA_FAST

#include "config.h"

#ifndef HAVE_GSL
#error "The GSL library is required for neurons that require a numerical solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"


/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
**/
extern "C" inline int elif_psc_alpha_fast_dynamics( double, const double y[], double f[], void* pnode );


/* BeginDocumentation
  Name: elif_psc_alpha_fast.

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
C_m [pF]  Membrane capacitance
g_L [nS]  leak conductance
E_0 [mV]  resting potential
E_u [mV]  upper potential
E_d [mV]  energy depletion potential
E_f [mV]  energy inflexion potential
epsilon_0 [real]  standard resting energy level
epsilon_c [real]  standard resting energy level
alpha [real]  energetic health
delta [real]  energy consumption per spike
tau_e [ms]  time constant for energy production
I_e [pA]  Constant input current
V_th [mV]  Spike detection threshold (reset condition)
V_reset [mV]  Reset potential
tau_syn_ex [ms]  Synaptic Time Constant Excitatory Synapse
tau_syn_in [ms]  Synaptic Time Constant for Inhibitory Synapse
t_ref [ms]  Refractory period


  Dynamic state variables:
r [integer]  number of steps for refractory phase
V_m [mV]  Membrane potential
epsilon [real]  Energy


  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class elif_psc_alpha_fast : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  elif_psc_alpha_fast();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c calibrate().
  **/
  elif_psc_alpha_fast(const elif_psc_alpha_fast &);

  /**
   * Destructor.
  **/
  ~elif_psc_alpha_fast();

  // -------------------------------------------------------------------------
  //   Import sets of overloaded virtual functions.
  //   See: Technical Issues / Virtual Functions: Overriding, Overloading,
  //        and Hiding
  // -------------------------------------------------------------------------

  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
   * Used to validate that we can send nest::SpikeEvent to desired target:port.
  **/
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  // -------------------------------------------------------------------------
  //   Functions handling incoming events.
  //   We tell nest that we can handle incoming events of various types by
  //   defining handle() for the given event.
  // -------------------------------------------------------------------------


  void handle(nest::SpikeEvent &);        //! accept spikes
  void handle(nest::CurrentEvent &);      //! accept input current
  void handle(nest::DataLoggingRequest &);//! allow recording with multimeter
  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:

  /**
   * Reset internal buffers of neuron.
  **/
  void init_buffers_();

  /**
   * Initialize auxiliary quantities, leave parameters and state untouched.
  **/
  void calibrate();

  /**
   * Take neuron through given time interval
  **/
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<elif_psc_alpha_fast>;
  friend class nest::UniversalDataLogger<elif_psc_alpha_fast>;

  /**
   * Free parameters of the neuron.
   *
   *
   *
   * These are the parameters that can be set by the user through @c `node.set()`.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update() and are not reset by
   * @c ResetNetwork.
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct Parameters_
  {    
    //!  Membrane capacitance
    double C_m;
    //!  leak conductance
    double g_L;
    //!  resting potential
    double E_0;
    //!  upper potential
    double E_u;
    //!  energy depletion potential
    double E_d;
    //!  energy inflexion potential
    double E_f;
    //!  standard resting energy level
    double epsilon_0;
    //!  standard resting energy level
    double epsilon_c;
    //!  energetic health
    double alpha;
    //!  energy consumption per spike
    double delta;
    //!  time constant for energy production
    double tau_e;
    //!  Constant input current
    double I_e;
    //!  Spike detection threshold (reset condition)
    double V_th;
    //!  Reset potential
    double V_reset;
    //!  Synaptic Time Constant Excitatory Synapse
    double tau_syn_ex;
    //!  Synaptic Time Constant for Inhibitory Synapse
    double tau_syn_in;
    //!  Refractory period
    double t_ref;

    double __gsl_error_tol;

    /**
     * Initialize parameters to their default values.
    **/
    Parameters_();
  };

  /**
   * Dynamic state of the neuron.
   *
   *
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c `node.set()`. The state variables are initialized from the model
   * prototype when the node is created. State variables are reset by @c ResetNetwork.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct State_
  {
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      epsilon,
      I_syn_ex__X__spikesExc,
      I_syn_ex__X__spikesExc__d,
      I_syn_in__X__spikesInh,
      I_syn_in__X__spikesInh__d,
      r,
      STATE_VEC_SIZE
    };

    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];

    State_();
  };

  /**
   * Internal variables of the neuron.
   *
   *
   *
   * These variables must be initialized by @c calibrate, which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c calibrate(). If Variables_ has members that
   *       cannot destroy themselves, Variables_ will need a destructor.
  **/
  struct Variables_
  {
    //!  refractory time in steps
    long RefractoryCounts;
    double __h;
    double invae;
    double inveps;
    double invEdEf;
    double invte;
    double invCm;
    double __P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc;
    double __P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d;
    double __P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc;
    double __P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d;
    double __P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh;
    double __P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d;
    double __P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh;
    double __P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d;
  };

  /**
   * Buffers of the neuron.
   * Usually buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c calibrate() on the first call to @c Simulate after the start of NEST,
   * ResetKernel or ResetNetwork.
   * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
   *       cannot destroy themselves, Buffers_ will need a destructor.
  **/
  struct Buffers_
  {
    Buffers_(elif_psc_alpha_fast &);
    Buffers_(const Buffers_ &, elif_psc_alpha_fast &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<elif_psc_alpha_fast> logger_;

    inline nest::RingBuffer& get_spikesInh() {return spikesInh;}
    //!< Buffer for input (type: pA)
    nest::RingBuffer spikesInh;
    double spikesInh_grid_sum_;
    inline nest::RingBuffer& get_spikesExc() {return spikesExc;}
    //!< Buffer for input (type: pA)
    nest::RingBuffer spikesExc;
    double spikesExc_grid_sum_;
    //!< Buffer for input (type: pA)
    nest::RingBuffer currents;
    inline nest::RingBuffer& get_currents() {return currents;}
    double currents_grid_sum_;

    // -----------------------------------------------------------------------
    //   GSL ODE solver data structures
    // -----------------------------------------------------------------------
    
    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // __integration_step should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline long get_r() const
  {
    return S_.ode_state[State_::r];
  }

  inline void set_r(const long __v)
  {
    S_.ode_state[State_::r] = __v;
  }

  inline double get_V_m() const
  {
    return S_.ode_state[State_::V_m];
  }

  inline void set_V_m(const double __v)
  {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_epsilon() const
  {
    return S_.ode_state[State_::epsilon];
  }

  inline void set_epsilon(const double __v)
  {
    S_.ode_state[State_::epsilon] = __v;
  }

  inline double get_I_syn_ex__X__spikesExc() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikesExc];
  }

  inline void set_I_syn_ex__X__spikesExc(const double __v)
  {
    S_.ode_state[State_::I_syn_ex__X__spikesExc] = __v;
  }

  inline double get_I_syn_ex__X__spikesExc__d() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikesExc__d];
  }

  inline void set_I_syn_ex__X__spikesExc__d(const double __v)
  {
    S_.ode_state[State_::I_syn_ex__X__spikesExc__d] = __v;
  }

  inline double get_I_syn_in__X__spikesInh() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikesInh];
  }

  inline void set_I_syn_in__X__spikesInh(const double __v)
  {
    S_.ode_state[State_::I_syn_in__X__spikesInh] = __v;
  }

  inline double get_I_syn_in__X__spikesInh__d() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikesInh__d];
  }

  inline void set_I_syn_in__X__spikesInh__d(const double __v)
  {
    S_.ode_state[State_::I_syn_in__X__spikesInh__d] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_C_m() const
  {
    return P_.C_m;
  }

  inline void set_C_m(const double __v)
  {
    P_.C_m = __v;
  }

  inline double get_g_L() const
  {
    return P_.g_L;
  }

  inline void set_g_L(const double __v)
  {
    P_.g_L = __v;
  }

  inline double get_E_0() const
  {
    return P_.E_0;
  }

  inline void set_E_0(const double __v)
  {
    P_.E_0 = __v;
  }

  inline double get_E_u() const
  {
    return P_.E_u;
  }

  inline void set_E_u(const double __v)
  {
    P_.E_u = __v;
  }

  inline double get_E_d() const
  {
    return P_.E_d;
  }

  inline void set_E_d(const double __v)
  {
    P_.E_d = __v;
  }

  inline double get_E_f() const
  {
    return P_.E_f;
  }

  inline void set_E_f(const double __v)
  {
    P_.E_f = __v;
  }

  inline double get_epsilon_0() const
  {
    return P_.epsilon_0;
  }

  inline void set_epsilon_0(const double __v)
  {
    P_.epsilon_0 = __v;
  }

  inline double get_epsilon_c() const
  {
    return P_.epsilon_c;
  }

  inline void set_epsilon_c(const double __v)
  {
    P_.epsilon_c = __v;
  }

  inline double get_alpha() const
  {
    return P_.alpha;
  }

  inline void set_alpha(const double __v)
  {
    P_.alpha = __v;
  }

  inline double get_delta() const
  {
    return P_.delta;
  }

  inline void set_delta(const double __v)
  {
    P_.delta = __v;
  }

  inline double get_tau_e() const
  {
    return P_.tau_e;
  }

  inline void set_tau_e(const double __v)
  {
    P_.tau_e = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_V_reset() const
  {
    return P_.V_reset;
  }

  inline void set_V_reset(const double __v)
  {
    P_.V_reset = __v;
  }

  inline double get_tau_syn_ex() const
  {
    return P_.tau_syn_ex;
  }

  inline void set_tau_syn_ex(const double __v)
  {
    P_.tau_syn_ex = __v;
  }

  inline double get_tau_syn_in() const
  {
    return P_.tau_syn_in;
  }

  inline void set_tau_syn_in(const double __v)
  {
    P_.tau_syn_in = __v;
  }

  inline double get_t_ref() const
  {
    return P_.t_ref;
  }

  inline void set_t_ref(const double __v)
  {
    P_.t_ref = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

  inline long get_RefractoryCounts() const
  {
    return V_.RefractoryCounts;
  }

  inline void set_RefractoryCounts(const long __v)
  {
    V_.RefractoryCounts = __v;
  }

  inline double get___h() const
  {
    return V_.__h;
  }

  inline void set___h(const double __v)
  {
    V_.__h = __v;
  }

  inline double get_invae() const
  {
    return V_.invae;
  }

  inline void set_invae(const double __v)
  {
    V_.invae = __v;
  }

  inline double get_inveps() const
  {
    return V_.inveps;
  }

  inline void set_inveps(const double __v)
  {
    V_.inveps = __v;
  }

  inline double get_invEdEf() const
  {
    return V_.invEdEf;
  }

  inline void set_invEdEf(const double __v)
  {
    V_.invEdEf = __v;
  }

  inline double get_invte() const
  {
    return V_.invte;
  }

  inline void set_invte(const double __v)
  {
    V_.invte = __v;
  }

  inline double get_invCm() const
  {
    return V_.invCm;
  }

  inline void set_invCm(const double __v)
  {
    V_.invCm = __v;
  }

  inline double get___P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc() const
  {
    return V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc;
  }

  inline void set___P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc(const double __v)
  {
    V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc = __v;
  }

  inline double get___P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d() const
  {
    return V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d;
  }

  inline void set___P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d(const double __v)
  {
    V_.__P__I_syn_ex__X__spikesExc__I_syn_ex__X__spikesExc__d = __v;
  }

  inline double get___P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc() const
  {
    return V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc;
  }

  inline void set___P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc(const double __v)
  {
    V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc = __v;
  }

  inline double get___P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d() const
  {
    return V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d;
  }

  inline void set___P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d(const double __v)
  {
    V_.__P__I_syn_ex__X__spikesExc__d__I_syn_ex__X__spikesExc__d = __v;
  }

  inline double get___P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh() const
  {
    return V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh;
  }

  inline void set___P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh(const double __v)
  {
    V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh = __v;
  }

  inline double get___P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d() const
  {
    return V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d;
  }

  inline void set___P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d(const double __v)
  {
    V_.__P__I_syn_in__X__spikesInh__I_syn_in__X__spikesInh__d = __v;
  }

  inline double get___P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh() const
  {
    return V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh;
  }

  inline void set___P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh(const double __v)
  {
    V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh = __v;
  }

  inline double get___P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d() const
  {
    return V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d;
  }

  inline void set___P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d(const double __v)
  {
    V_.__P__I_syn_in__X__spikesInh__d__I_syn_in__X__spikesInh__d = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  
  inline double get_eps_bound() const
  {
    return std::max(S_.ode_state[State_::epsilon], 0.0);
  }

  inline double get_E_L() const
  {
    return P_.E_0 + (P_.E_u - P_.E_0) * (1 - (std::max(S_.ode_state[State_::epsilon], 0.0)) * V_.inveps);
  }

  inline double get_I_in() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikesInh];
  }

  inline double get_I_ex() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikesExc];
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------
  
  inline nest::RingBuffer& get_spikesInh() {return B_.get_spikesInh();};
  inline nest::RingBuffer& get_spikesExc() {return B_.get_spikesExc();};
  inline nest::RingBuffer& get_currents() {return B_.get_currents();};

  // -------------------------------------------------------------------------
  //   Member variables of neuron model.
  //   Each model neuron should have precisely the following four data members,
  //   which are one instance each of the parameters, state, buffers and variables
  //   structures. Experience indicates that the state and variables member should
  //   be next to each other to achieve good efficiency (caching).
  //   Note: Devices require one additional data member, an instance of the 
  //   ``Device`` child class they belong to.
  // -------------------------------------------------------------------------

  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<elif_psc_alpha_fast> recordablesMap_;
  friend int elif_psc_alpha_fast_dynamics( double, const double y[], double f[], void* pnode );

}; /* neuron elif_psc_alpha_fast */

inline nest::port elif_psc_alpha_fast::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port elif_psc_alpha_fast::handles_test_event(nest::SpikeEvent&, nest::port receptor_type)
{
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
    {
      throw nest::UnknownReceptorType(receptor_type, get_name());
    }
    return 0;
}

inline nest::port elif_psc_alpha_fast::handles_test_event(nest::CurrentEvent&, nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }
  return 0;
}

inline nest::port elif_psc_alpha_fast::handles_test_event(nest::DataLoggingRequest& dlr, nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void elif_psc_alpha_fast::get_status(DictionaryDatum &__d) const
{
  // parameters
  def<double>(__d, "C_m", get_C_m());
  def<double>(__d, "g_L", get_g_L());
  def<double>(__d, "E_0", get_E_0());
  def<double>(__d, "E_u", get_E_u());
  def<double>(__d, "E_d", get_E_d());
  def<double>(__d, "E_f", get_E_f());
  def<double>(__d, "epsilon_0", get_epsilon_0());
  def<double>(__d, "epsilon_c", get_epsilon_c());
  def<double>(__d, "alpha", get_alpha());
  def<double>(__d, "delta", get_delta());
  def<double>(__d, "tau_e", get_tau_e());
  def<double>(__d, "I_e", get_I_e());
  def<double>(__d, "V_th", get_V_th());
  def<double>(__d, "V_reset", get_V_reset());
  def<double>(__d, "tau_syn_ex", get_tau_syn_ex());
  def<double>(__d, "tau_syn_in", get_tau_syn_in());
  def<double>(__d, "t_ref", get_t_ref());

  // initial values for state variables in ODE or kernel
  def<long>(__d, "r", get_r());
  def<double>(__d, "V_m", get_V_m());
  def<double>(__d, "epsilon", get_epsilon());
  def<double>(__d, "I_syn_ex__X__spikesExc", get_I_syn_ex__X__spikesExc());
  def<double>(__d, "I_syn_ex__X__spikesExc__d", get_I_syn_ex__X__spikesExc__d());
  def<double>(__d, "I_syn_in__X__spikesInh", get_I_syn_in__X__spikesInh());
  def<double>(__d, "I_syn_in__X__spikesInh__d", get_I_syn_in__X__spikesInh__d());

  ArchivingNode::get_status( __d );

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void elif_psc_alpha_fast::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);
  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);
  double tmp_E_0 = get_E_0();
  updateValue<double>(__d, "E_0", tmp_E_0);
  double tmp_E_u = get_E_u();
  updateValue<double>(__d, "E_u", tmp_E_u);
  double tmp_E_d = get_E_d();
  updateValue<double>(__d, "E_d", tmp_E_d);
  double tmp_E_f = get_E_f();
  updateValue<double>(__d, "E_f", tmp_E_f);
  double tmp_epsilon_0 = get_epsilon_0();
  updateValue<double>(__d, "epsilon_0", tmp_epsilon_0);
  double tmp_epsilon_c = get_epsilon_c();
  updateValue<double>(__d, "epsilon_c", tmp_epsilon_c);
  double tmp_alpha = get_alpha();
  updateValue<double>(__d, "alpha", tmp_alpha);
  double tmp_delta = get_delta();
  updateValue<double>(__d, "delta", tmp_delta);
  double tmp_tau_e = get_tau_e();
  updateValue<double>(__d, "tau_e", tmp_tau_e);
  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);
  double tmp_V_th = get_V_th();
  updateValue<double>(__d, "V_th", tmp_V_th);
  double tmp_V_reset = get_V_reset();
  updateValue<double>(__d, "V_reset", tmp_V_reset);
  double tmp_tau_syn_ex = get_tau_syn_ex();
  updateValue<double>(__d, "tau_syn_ex", tmp_tau_syn_ex);
  double tmp_tau_syn_in = get_tau_syn_in();
  updateValue<double>(__d, "tau_syn_in", tmp_tau_syn_in);
  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);

  // initial values for state variables in ODE or kernel
  long tmp_r = get_r();
  updateValue<long>(__d, "r", tmp_r);
  double tmp_V_m = get_V_m();
  updateValue<double>(__d, "V_m", tmp_V_m);
  double tmp_epsilon = get_epsilon();
  updateValue<double>(__d, "epsilon", tmp_epsilon);
  double tmp_I_syn_ex__X__spikesExc = get_I_syn_ex__X__spikesExc();
  updateValue<double>(__d, "I_syn_ex__X__spikesExc", tmp_I_syn_ex__X__spikesExc);
  double tmp_I_syn_ex__X__spikesExc__d = get_I_syn_ex__X__spikesExc__d();
  updateValue<double>(__d, "I_syn_ex__X__spikesExc__d", tmp_I_syn_ex__X__spikesExc__d);
  double tmp_I_syn_in__X__spikesInh = get_I_syn_in__X__spikesInh();
  updateValue<double>(__d, "I_syn_in__X__spikesInh", tmp_I_syn_in__X__spikesInh);
  double tmp_I_syn_in__X__spikesInh__d = get_I_syn_in__X__spikesInh__d();
  updateValue<double>(__d, "I_syn_in__X__spikesInh__d", tmp_I_syn_in__X__spikesInh__d);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_C_m(tmp_C_m);
  set_g_L(tmp_g_L);
  set_E_0(tmp_E_0);
  set_E_u(tmp_E_u);
  set_E_d(tmp_E_d);
  set_E_f(tmp_E_f);
  set_epsilon_0(tmp_epsilon_0);
  set_epsilon_c(tmp_epsilon_c);
  set_alpha(tmp_alpha);
  set_delta(tmp_delta);
  set_tau_e(tmp_tau_e);
  set_I_e(tmp_I_e);
  set_V_th(tmp_V_th);
  set_V_reset(tmp_V_reset);
  set_tau_syn_ex(tmp_tau_syn_ex);
  set_tau_syn_in(tmp_tau_syn_in);
  set_t_ref(tmp_t_ref);
  set_r(tmp_r);
  set_V_m(tmp_V_m);
  set_epsilon(tmp_epsilon);
  set_I_syn_ex__X__spikesExc(tmp_I_syn_ex__X__spikesExc);
  set_I_syn_ex__X__spikesExc__d(tmp_I_syn_ex__X__spikesExc__d);
  set_I_syn_in__X__spikesInh(tmp_I_syn_in__X__spikesInh);
  set_I_syn_in__X__spikesInh__d(tmp_I_syn_in__X__spikesInh__d);


  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
};

#endif /* #ifndef ELIF_PSC_ALPHA_FAST */
