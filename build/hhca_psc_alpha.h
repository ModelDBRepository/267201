/**
 *  hhca_psc_alpha.h
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
 *  Generated from NESTML at time: 2021-09-17 08:49:06.519180
**/
#ifndef HHCA_PSC_ALPHA
#define HHCA_PSC_ALPHA

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
extern "C" inline int hhca_psc_alpha_dynamics( double, const double y[], double f[], void* pnode );


/* BeginDocumentation
  Name: hhca_psc_alpha.

  Description:

    """
  Name: hhca_psc_alpha - Hodgkin Huxley neuron model.

  Description:

   hhca_psc_alpha is an implementation of a spiking neuron using the Hodkin-Huxley
   formalism.

   (1) Post-syaptic currents
   Incoming spike events induce a post-synaptic change of current modelled
   by an alpha function. The alpha function is normalised such that an event of
   weight 1.0 results in a peak current of 1 pA.


   (2) Spike Detection
   Spike detection is done by a combined threshold-and-local-maximum search: if
   there is a local maximum above a certain threshold of the membrane potential,
   it is considered a spike.

  Problems/Todo:

   better spike detection
   initial wavelet/spike at simulation onset

  References:

   Spiking Neuron Models:
   Single Neurons, Populations, Plasticity
   Wulfram Gerstner, Werner Kistler,  Cambridge University Press

   Theoretical Neuroscience:
   Computational and Mathematical Modeling of Neural Systems
   Peter Dayan, L. F. Abbott, MIT Press (parameters taken from here)

   Hodgkin, A. L. and Huxley, A. F.,
   A Quantitative Description of Membrane Current
   and Its Application to Conduction and Excitation in Nerve,
   Journal of Physiology, 117, 500-544 (1952)

  Sends: SpikeEvent

  Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

  Authors: Tanguy Fardet
  SeeAlso: hh_psc_alpha
  """


  Parameters:
  The following parameters can be set in the status dictionary.
t_ref [ms]  Refractory period
g_Na [nS]  Sodium peak conductance
g_K [nS]  Potassium peak conductance
g_L [nS]  Leak conductance
C_m [pF]  Membrane Capacitance
E_Na [mV]  Sodium reversal potential
E_K [mV]  Potassium reversal potential
E_L [mV]  Leak reversal Potential (aka resting potential)
tau_syn_ex [ms]  Rise time of the excitatory synaptic alpha function i
tau_syn_in [ms]  Rise time of the inhibitory synaptic alpha function
I_e [pA]  constant external input current


  Dynamic state variables:
r [integer]  number of steps in the current refractory phase
V_m [mV]  Membrane potential
Act_m [real]  Activation variable m for Na
Inact_h [real]  Inactivation variable h for Na
Act_n [real]  Activation variable n for K


  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class hhca_psc_alpha : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  hhca_psc_alpha();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c calibrate().
  **/
  hhca_psc_alpha(const hhca_psc_alpha &);

  /**
   * Destructor.
  **/
  ~hhca_psc_alpha();

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
  friend class nest::RecordablesMap<hhca_psc_alpha>;
  friend class nest::UniversalDataLogger<hhca_psc_alpha>;

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
    //!  Refractory period
    double t_ref;
    //!  Sodium peak conductance
    double g_Na;
    //!  Potassium peak conductance
    double g_K;
    //!  Leak conductance
    double g_L;
    //!  Membrane Capacitance
    double C_m;
    //!  Sodium reversal potential
    double E_Na;
    //!  Potassium reversal potential
    double E_K;
    //!  Leak reversal Potential (aka resting potential)
    double E_L;
    //!  Rise time of the excitatory synaptic alpha function i
    double tau_syn_ex;
    //!  Rise time of the inhibitory synaptic alpha function
    double tau_syn_in;
    double g_Ca;
    double Ca_env;
    double Ca_i_eq;
    double tau_Ca;
    double V_hmCa;
    double k_mCa;
    double tau_mCa;
    double V_hhCa;
    double k_hCa;
    double tau_hCa;
    double g_AHP;
    double K_AHP;
    double b_AHP;
    double tau_AHP;
    double E_0;
    double Vcell;
    //!  constant external input current
    double I_e;

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
      Ca_i,
      Act_n,
      Act_m,
      Inact_h,
      h_Ca,
      m_Ca,
      m_AHP,
      I_syn_ex__X__spikeExc,
      I_syn_ex__X__spikeExc__d,
      I_syn_in__X__spikeInh,
      I_syn_in__X__spikeInh__d,
      r,
      alpha_n_init,
      beta_n_init,
      alpha_m_init,
      beta_m_init,
      alpha_h_init,
      beta_h_init,
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
    double charge;
    double avogadro;
    //!  convert charge flux to concentration roughly 1e-2 mmol/L/fC
    double c2c;
    double __P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc;
    double __P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d;
    double __P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc;
    double __P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d;
    double __P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh;
    double __P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d;
    double __P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh;
    double __P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d;
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
    Buffers_(hhca_psc_alpha &);
    Buffers_(const Buffers_ &, hhca_psc_alpha &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<hhca_psc_alpha> logger_;

    inline nest::RingBuffer& get_spikeInh() {return spikeInh;}
    //!< Buffer for input (type: pA)
    nest::RingBuffer spikeInh;
    double spikeInh_grid_sum_;
    inline nest::RingBuffer& get_spikeExc() {return spikeExc;}
    //!< Buffer for input (type: pA)
    nest::RingBuffer spikeExc;
    double spikeExc_grid_sum_;
    //!< Buffer for input (type: pA)
    nest::RingBuffer I_stim;
    inline nest::RingBuffer& get_I_stim() {return I_stim;}
    double I_stim_grid_sum_;

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

  inline double get_Ca_i() const
  {
    return S_.ode_state[State_::Ca_i];
  }

  inline void set_Ca_i(const double __v)
  {
    S_.ode_state[State_::Ca_i] = __v;
  }

  inline double get_alpha_n_init() const
  {
    return S_.ode_state[State_::alpha_n_init];
  }

  inline void set_alpha_n_init(const double __v)
  {
    S_.ode_state[State_::alpha_n_init] = __v;
  }

  inline double get_beta_n_init() const
  {
    return S_.ode_state[State_::beta_n_init];
  }

  inline void set_beta_n_init(const double __v)
  {
    S_.ode_state[State_::beta_n_init] = __v;
  }

  inline double get_alpha_m_init() const
  {
    return S_.ode_state[State_::alpha_m_init];
  }

  inline void set_alpha_m_init(const double __v)
  {
    S_.ode_state[State_::alpha_m_init] = __v;
  }

  inline double get_beta_m_init() const
  {
    return S_.ode_state[State_::beta_m_init];
  }

  inline void set_beta_m_init(const double __v)
  {
    S_.ode_state[State_::beta_m_init] = __v;
  }

  inline double get_alpha_h_init() const
  {
    return S_.ode_state[State_::alpha_h_init];
  }

  inline void set_alpha_h_init(const double __v)
  {
    S_.ode_state[State_::alpha_h_init] = __v;
  }

  inline double get_beta_h_init() const
  {
    return S_.ode_state[State_::beta_h_init];
  }

  inline void set_beta_h_init(const double __v)
  {
    S_.ode_state[State_::beta_h_init] = __v;
  }

  inline double get_Act_m() const
  {
    return S_.ode_state[State_::Act_m];
  }

  inline void set_Act_m(const double __v)
  {
    S_.ode_state[State_::Act_m] = __v;
  }

  inline double get_Inact_h() const
  {
    return S_.ode_state[State_::Inact_h];
  }

  inline void set_Inact_h(const double __v)
  {
    S_.ode_state[State_::Inact_h] = __v;
  }

  inline double get_Act_n() const
  {
    return S_.ode_state[State_::Act_n];
  }

  inline void set_Act_n(const double __v)
  {
    S_.ode_state[State_::Act_n] = __v;
  }

  inline double get_h_Ca() const
  {
    return S_.ode_state[State_::h_Ca];
  }

  inline void set_h_Ca(const double __v)
  {
    S_.ode_state[State_::h_Ca] = __v;
  }

  inline double get_m_Ca() const
  {
    return S_.ode_state[State_::m_Ca];
  }

  inline void set_m_Ca(const double __v)
  {
    S_.ode_state[State_::m_Ca] = __v;
  }

  inline double get_m_AHP() const
  {
    return S_.ode_state[State_::m_AHP];
  }

  inline void set_m_AHP(const double __v)
  {
    S_.ode_state[State_::m_AHP] = __v;
  }

  inline double get_I_syn_ex__X__spikeExc() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikeExc];
  }

  inline void set_I_syn_ex__X__spikeExc(const double __v)
  {
    S_.ode_state[State_::I_syn_ex__X__spikeExc] = __v;
  }

  inline double get_I_syn_ex__X__spikeExc__d() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikeExc__d];
  }

  inline void set_I_syn_ex__X__spikeExc__d(const double __v)
  {
    S_.ode_state[State_::I_syn_ex__X__spikeExc__d] = __v;
  }

  inline double get_I_syn_in__X__spikeInh() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikeInh];
  }

  inline void set_I_syn_in__X__spikeInh(const double __v)
  {
    S_.ode_state[State_::I_syn_in__X__spikeInh] = __v;
  }

  inline double get_I_syn_in__X__spikeInh__d() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikeInh__d];
  }

  inline void set_I_syn_in__X__spikeInh__d(const double __v)
  {
    S_.ode_state[State_::I_syn_in__X__spikeInh__d] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_t_ref() const
  {
    return P_.t_ref;
  }

  inline void set_t_ref(const double __v)
  {
    P_.t_ref = __v;
  }

  inline double get_g_Na() const
  {
    return P_.g_Na;
  }

  inline void set_g_Na(const double __v)
  {
    P_.g_Na = __v;
  }

  inline double get_g_K() const
  {
    return P_.g_K;
  }

  inline void set_g_K(const double __v)
  {
    P_.g_K = __v;
  }

  inline double get_g_L() const
  {
    return P_.g_L;
  }

  inline void set_g_L(const double __v)
  {
    P_.g_L = __v;
  }

  inline double get_C_m() const
  {
    return P_.C_m;
  }

  inline void set_C_m(const double __v)
  {
    P_.C_m = __v;
  }

  inline double get_E_Na() const
  {
    return P_.E_Na;
  }

  inline void set_E_Na(const double __v)
  {
    P_.E_Na = __v;
  }

  inline double get_E_K() const
  {
    return P_.E_K;
  }

  inline void set_E_K(const double __v)
  {
    P_.E_K = __v;
  }

  inline double get_E_L() const
  {
    return P_.E_L;
  }

  inline void set_E_L(const double __v)
  {
    P_.E_L = __v;
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

  inline double get_g_Ca() const
  {
    return P_.g_Ca;
  }

  inline void set_g_Ca(const double __v)
  {
    P_.g_Ca = __v;
  }

  inline double get_Ca_env() const
  {
    return P_.Ca_env;
  }

  inline void set_Ca_env(const double __v)
  {
    P_.Ca_env = __v;
  }

  inline double get_Ca_i_eq() const
  {
    return P_.Ca_i_eq;
  }

  inline void set_Ca_i_eq(const double __v)
  {
    P_.Ca_i_eq = __v;
  }

  inline double get_tau_Ca() const
  {
    return P_.tau_Ca;
  }

  inline void set_tau_Ca(const double __v)
  {
    P_.tau_Ca = __v;
  }

  inline double get_V_hmCa() const
  {
    return P_.V_hmCa;
  }

  inline void set_V_hmCa(const double __v)
  {
    P_.V_hmCa = __v;
  }

  inline double get_k_mCa() const
  {
    return P_.k_mCa;
  }

  inline void set_k_mCa(const double __v)
  {
    P_.k_mCa = __v;
  }

  inline double get_tau_mCa() const
  {
    return P_.tau_mCa;
  }

  inline void set_tau_mCa(const double __v)
  {
    P_.tau_mCa = __v;
  }

  inline double get_V_hhCa() const
  {
    return P_.V_hhCa;
  }

  inline void set_V_hhCa(const double __v)
  {
    P_.V_hhCa = __v;
  }

  inline double get_k_hCa() const
  {
    return P_.k_hCa;
  }

  inline void set_k_hCa(const double __v)
  {
    P_.k_hCa = __v;
  }

  inline double get_tau_hCa() const
  {
    return P_.tau_hCa;
  }

  inline void set_tau_hCa(const double __v)
  {
    P_.tau_hCa = __v;
  }

  inline double get_g_AHP() const
  {
    return P_.g_AHP;
  }

  inline void set_g_AHP(const double __v)
  {
    P_.g_AHP = __v;
  }

  inline double get_K_AHP() const
  {
    return P_.K_AHP;
  }

  inline void set_K_AHP(const double __v)
  {
    P_.K_AHP = __v;
  }

  inline double get_b_AHP() const
  {
    return P_.b_AHP;
  }

  inline void set_b_AHP(const double __v)
  {
    P_.b_AHP = __v;
  }

  inline double get_tau_AHP() const
  {
    return P_.tau_AHP;
  }

  inline void set_tau_AHP(const double __v)
  {
    P_.tau_AHP = __v;
  }

  inline double get_E_0() const
  {
    return P_.E_0;
  }

  inline void set_E_0(const double __v)
  {
    P_.E_0 = __v;
  }

  inline double get_Vcell() const
  {
    return P_.Vcell;
  }

  inline void set_Vcell(const double __v)
  {
    P_.Vcell = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
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

  inline double get_charge() const
  {
    return V_.charge;
  }

  inline void set_charge(const double __v)
  {
    V_.charge = __v;
  }

  inline double get_avogadro() const
  {
    return V_.avogadro;
  }

  inline void set_avogadro(const double __v)
  {
    V_.avogadro = __v;
  }

  inline double get_c2c() const
  {
    return V_.c2c;
  }

  inline void set_c2c(const double __v)
  {
    V_.c2c = __v;
  }

  inline double get___P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc() const
  {
    return V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc;
  }

  inline void set___P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc(const double __v)
  {
    V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc = __v;
  }

  inline double get___P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d() const
  {
    return V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d;
  }

  inline void set___P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d(const double __v)
  {
    V_.__P__I_syn_ex__X__spikeExc__I_syn_ex__X__spikeExc__d = __v;
  }

  inline double get___P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc() const
  {
    return V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc;
  }

  inline void set___P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc(const double __v)
  {
    V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc = __v;
  }

  inline double get___P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d() const
  {
    return V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d;
  }

  inline void set___P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d(const double __v)
  {
    V_.__P__I_syn_ex__X__spikeExc__d__I_syn_ex__X__spikeExc__d = __v;
  }

  inline double get___P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh() const
  {
    return V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh;
  }

  inline void set___P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh(const double __v)
  {
    V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh = __v;
  }

  inline double get___P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d() const
  {
    return V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d;
  }

  inline void set___P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d(const double __v)
  {
    V_.__P__I_syn_in__X__spikeInh__I_syn_in__X__spikeInh__d = __v;
  }

  inline double get___P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh() const
  {
    return V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh;
  }

  inline void set___P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh(const double __v)
  {
    V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh = __v;
  }

  inline double get___P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d() const
  {
    return V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d;
  }

  inline void set___P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d(const double __v)
  {
    V_.__P__I_syn_in__X__spikeInh__d__I_syn_in__X__spikeInh__d = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  
  inline double get_m_infAHP() const
  {
    return pow((S_.ode_state[State_::Ca_i] / P_.K_AHP), 2) / (pow((S_.ode_state[State_::Ca_i] / P_.K_AHP), 2) + P_.b_AHP);
  }

  inline double get_m_infCa() const
  {
    return y_shape(S_.ode_state[State_::V_m], P_.V_hmCa, P_.k_mCa);
  }

  inline double get_h_infCa() const
  {
    return y_shape(S_.ode_state[State_::V_m], P_.V_hhCa, (-P_.k_hCa));
  }

  inline double get_E_Ca() const
  {
    return 2 * P_.E_0 * std::log(P_.Ca_env / S_.ode_state[State_::Ca_i]);
  }

  inline double get_I_syn_exc() const
  {
    return S_.ode_state[State_::I_syn_ex__X__spikeExc];
  }

  inline double get_I_syn_inh() const
  {
    return S_.ode_state[State_::I_syn_in__X__spikeInh];
  }

  inline double get_I_Na() const
  {
    return P_.g_Na * S_.ode_state[State_::Act_m] * S_.ode_state[State_::Act_m] * S_.ode_state[State_::Act_m] * S_.ode_state[State_::Inact_h] * (S_.ode_state[State_::V_m] - P_.E_Na);
  }

  inline double get_I_K() const
  {
    return (P_.g_K * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] + P_.g_AHP * pow(S_.ode_state[State_::m_AHP], 2)) * (S_.ode_state[State_::V_m] - P_.E_K);
  }

  inline double get_I_Ca() const
  {
    return (-P_.g_Ca) * S_.ode_state[State_::m_Ca] * S_.ode_state[State_::h_Ca] * (S_.ode_state[State_::V_m] - (2 * P_.E_0 * std::log(P_.Ca_env / S_.ode_state[State_::Ca_i])));
  }

  inline double get_I_L() const
  {
    return P_.g_L * (S_.ode_state[State_::V_m] - P_.E_L);
  }

  inline double get_alpha_n() const
  {
    return (0.01 * (S_.ode_state[State_::V_m] / 1.0 + 55.0)) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 55.0)) / 10.0));
  }

  inline double get_beta_n() const
  {
    return 0.125 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 80.0);
  }

  inline double get_alpha_m() const
  {
    return (0.1 * (S_.ode_state[State_::V_m] / 1.0 + 40.0)) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 40.0)) / 10.0));
  }

  inline double get_beta_m() const
  {
    return 4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0);
  }

  inline double get_alpha_h() const
  {
    return 0.07 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 20.0);
  }

  inline double get_beta_h() const
  {
    return 1.0 / (1.0 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 35.0)) / 10.0));
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------
  
  inline nest::RingBuffer& get_spikeInh() {return B_.get_spikeInh();};
  inline nest::RingBuffer& get_spikeExc() {return B_.get_spikeExc();};
  inline nest::RingBuffer& get_I_stim() {return B_.get_I_stim();};
  // -------------------------------------------------------------------------
  //   Function declarations
  // -------------------------------------------------------------------------


  //
  double y_shape(double, double, double) const
  ;

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
  static nest::RecordablesMap<hhca_psc_alpha> recordablesMap_;
  friend int hhca_psc_alpha_dynamics( double, const double y[], double f[], void* pnode );

}; /* neuron hhca_psc_alpha */

inline nest::port hhca_psc_alpha::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port hhca_psc_alpha::handles_test_event(nest::SpikeEvent&, nest::port receptor_type)
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

inline nest::port hhca_psc_alpha::handles_test_event(nest::CurrentEvent&, nest::port receptor_type)
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

inline nest::port hhca_psc_alpha::handles_test_event(nest::DataLoggingRequest& dlr, nest::port receptor_type)
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

inline void hhca_psc_alpha::get_status(DictionaryDatum &__d) const
{
  // parameters
  def<double>(__d, "t_ref", get_t_ref());
  def<double>(__d, "g_Na", get_g_Na());
  def<double>(__d, "g_K", get_g_K());
  def<double>(__d, "g_L", get_g_L());
  def<double>(__d, "C_m", get_C_m());
  def<double>(__d, "E_Na", get_E_Na());
  def<double>(__d, "E_K", get_E_K());
  def<double>(__d, "E_L", get_E_L());
  def<double>(__d, "tau_syn_ex", get_tau_syn_ex());
  def<double>(__d, "tau_syn_in", get_tau_syn_in());
  def<double>(__d, "g_Ca", get_g_Ca());
  def<double>(__d, "Ca_env", get_Ca_env());
  def<double>(__d, "Ca_i_eq", get_Ca_i_eq());
  def<double>(__d, "tau_Ca", get_tau_Ca());
  def<double>(__d, "V_hmCa", get_V_hmCa());
  def<double>(__d, "k_mCa", get_k_mCa());
  def<double>(__d, "tau_mCa", get_tau_mCa());
  def<double>(__d, "V_hhCa", get_V_hhCa());
  def<double>(__d, "k_hCa", get_k_hCa());
  def<double>(__d, "tau_hCa", get_tau_hCa());
  def<double>(__d, "g_AHP", get_g_AHP());
  def<double>(__d, "K_AHP", get_K_AHP());
  def<double>(__d, "b_AHP", get_b_AHP());
  def<double>(__d, "tau_AHP", get_tau_AHP());
  def<double>(__d, "E_0", get_E_0());
  def<double>(__d, "Vcell", get_Vcell());
  def<double>(__d, "I_e", get_I_e());

  // initial values for state variables in ODE or kernel
  def<long>(__d, "r", get_r());
  def<double>(__d, "V_m", get_V_m());
  def<double>(__d, "Ca_i", get_Ca_i());
  def<double>(__d, "alpha_n_init", get_alpha_n_init());
  def<double>(__d, "beta_n_init", get_beta_n_init());
  def<double>(__d, "alpha_m_init", get_alpha_m_init());
  def<double>(__d, "beta_m_init", get_beta_m_init());
  def<double>(__d, "alpha_h_init", get_alpha_h_init());
  def<double>(__d, "beta_h_init", get_beta_h_init());
  def<double>(__d, "Act_m", get_Act_m());
  def<double>(__d, "Inact_h", get_Inact_h());
  def<double>(__d, "Act_n", get_Act_n());
  def<double>(__d, "h_Ca", get_h_Ca());
  def<double>(__d, "m_Ca", get_m_Ca());
  def<double>(__d, "m_AHP", get_m_AHP());
  def<double>(__d, "I_syn_ex__X__spikeExc", get_I_syn_ex__X__spikeExc());
  def<double>(__d, "I_syn_ex__X__spikeExc__d", get_I_syn_ex__X__spikeExc__d());
  def<double>(__d, "I_syn_in__X__spikeInh", get_I_syn_in__X__spikeInh());
  def<double>(__d, "I_syn_in__X__spikeInh__d", get_I_syn_in__X__spikeInh__d());

  ArchivingNode::get_status( __d );

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void hhca_psc_alpha::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);
  double tmp_g_Na = get_g_Na();
  updateValue<double>(__d, "g_Na", tmp_g_Na);
  double tmp_g_K = get_g_K();
  updateValue<double>(__d, "g_K", tmp_g_K);
  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);
  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);
  double tmp_E_Na = get_E_Na();
  updateValue<double>(__d, "E_Na", tmp_E_Na);
  double tmp_E_K = get_E_K();
  updateValue<double>(__d, "E_K", tmp_E_K);
  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);
  double tmp_tau_syn_ex = get_tau_syn_ex();
  updateValue<double>(__d, "tau_syn_ex", tmp_tau_syn_ex);
  double tmp_tau_syn_in = get_tau_syn_in();
  updateValue<double>(__d, "tau_syn_in", tmp_tau_syn_in);
  double tmp_g_Ca = get_g_Ca();
  updateValue<double>(__d, "g_Ca", tmp_g_Ca);
  double tmp_Ca_env = get_Ca_env();
  updateValue<double>(__d, "Ca_env", tmp_Ca_env);
  double tmp_Ca_i_eq = get_Ca_i_eq();
  updateValue<double>(__d, "Ca_i_eq", tmp_Ca_i_eq);
  double tmp_tau_Ca = get_tau_Ca();
  updateValue<double>(__d, "tau_Ca", tmp_tau_Ca);
  double tmp_V_hmCa = get_V_hmCa();
  updateValue<double>(__d, "V_hmCa", tmp_V_hmCa);
  double tmp_k_mCa = get_k_mCa();
  updateValue<double>(__d, "k_mCa", tmp_k_mCa);
  double tmp_tau_mCa = get_tau_mCa();
  updateValue<double>(__d, "tau_mCa", tmp_tau_mCa);
  double tmp_V_hhCa = get_V_hhCa();
  updateValue<double>(__d, "V_hhCa", tmp_V_hhCa);
  double tmp_k_hCa = get_k_hCa();
  updateValue<double>(__d, "k_hCa", tmp_k_hCa);
  double tmp_tau_hCa = get_tau_hCa();
  updateValue<double>(__d, "tau_hCa", tmp_tau_hCa);
  double tmp_g_AHP = get_g_AHP();
  updateValue<double>(__d, "g_AHP", tmp_g_AHP);
  double tmp_K_AHP = get_K_AHP();
  updateValue<double>(__d, "K_AHP", tmp_K_AHP);
  double tmp_b_AHP = get_b_AHP();
  updateValue<double>(__d, "b_AHP", tmp_b_AHP);
  double tmp_tau_AHP = get_tau_AHP();
  updateValue<double>(__d, "tau_AHP", tmp_tau_AHP);
  double tmp_E_0 = get_E_0();
  updateValue<double>(__d, "E_0", tmp_E_0);
  double tmp_Vcell = get_Vcell();
  updateValue<double>(__d, "Vcell", tmp_Vcell);
  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);

  // initial values for state variables in ODE or kernel
  long tmp_r = get_r();
  updateValue<long>(__d, "r", tmp_r);
  double tmp_V_m = get_V_m();
  updateValue<double>(__d, "V_m", tmp_V_m);
  double tmp_Ca_i = get_Ca_i();
  updateValue<double>(__d, "Ca_i", tmp_Ca_i);
  double tmp_alpha_n_init = get_alpha_n_init();
  updateValue<double>(__d, "alpha_n_init", tmp_alpha_n_init);
  double tmp_beta_n_init = get_beta_n_init();
  updateValue<double>(__d, "beta_n_init", tmp_beta_n_init);
  double tmp_alpha_m_init = get_alpha_m_init();
  updateValue<double>(__d, "alpha_m_init", tmp_alpha_m_init);
  double tmp_beta_m_init = get_beta_m_init();
  updateValue<double>(__d, "beta_m_init", tmp_beta_m_init);
  double tmp_alpha_h_init = get_alpha_h_init();
  updateValue<double>(__d, "alpha_h_init", tmp_alpha_h_init);
  double tmp_beta_h_init = get_beta_h_init();
  updateValue<double>(__d, "beta_h_init", tmp_beta_h_init);
  double tmp_Act_m = get_Act_m();
  updateValue<double>(__d, "Act_m", tmp_Act_m);
  double tmp_Inact_h = get_Inact_h();
  updateValue<double>(__d, "Inact_h", tmp_Inact_h);
  double tmp_Act_n = get_Act_n();
  updateValue<double>(__d, "Act_n", tmp_Act_n);
  double tmp_h_Ca = get_h_Ca();
  updateValue<double>(__d, "h_Ca", tmp_h_Ca);
  double tmp_m_Ca = get_m_Ca();
  updateValue<double>(__d, "m_Ca", tmp_m_Ca);
  double tmp_m_AHP = get_m_AHP();
  updateValue<double>(__d, "m_AHP", tmp_m_AHP);
  double tmp_I_syn_ex__X__spikeExc = get_I_syn_ex__X__spikeExc();
  updateValue<double>(__d, "I_syn_ex__X__spikeExc", tmp_I_syn_ex__X__spikeExc);
  double tmp_I_syn_ex__X__spikeExc__d = get_I_syn_ex__X__spikeExc__d();
  updateValue<double>(__d, "I_syn_ex__X__spikeExc__d", tmp_I_syn_ex__X__spikeExc__d);
  double tmp_I_syn_in__X__spikeInh = get_I_syn_in__X__spikeInh();
  updateValue<double>(__d, "I_syn_in__X__spikeInh", tmp_I_syn_in__X__spikeInh);
  double tmp_I_syn_in__X__spikeInh__d = get_I_syn_in__X__spikeInh__d();
  updateValue<double>(__d, "I_syn_in__X__spikeInh__d", tmp_I_syn_in__X__spikeInh__d);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_t_ref(tmp_t_ref);
  set_g_Na(tmp_g_Na);
  set_g_K(tmp_g_K);
  set_g_L(tmp_g_L);
  set_C_m(tmp_C_m);
  set_E_Na(tmp_E_Na);
  set_E_K(tmp_E_K);
  set_E_L(tmp_E_L);
  set_tau_syn_ex(tmp_tau_syn_ex);
  set_tau_syn_in(tmp_tau_syn_in);
  set_g_Ca(tmp_g_Ca);
  set_Ca_env(tmp_Ca_env);
  set_Ca_i_eq(tmp_Ca_i_eq);
  set_tau_Ca(tmp_tau_Ca);
  set_V_hmCa(tmp_V_hmCa);
  set_k_mCa(tmp_k_mCa);
  set_tau_mCa(tmp_tau_mCa);
  set_V_hhCa(tmp_V_hhCa);
  set_k_hCa(tmp_k_hCa);
  set_tau_hCa(tmp_tau_hCa);
  set_g_AHP(tmp_g_AHP);
  set_K_AHP(tmp_K_AHP);
  set_b_AHP(tmp_b_AHP);
  set_tau_AHP(tmp_tau_AHP);
  set_E_0(tmp_E_0);
  set_Vcell(tmp_Vcell);
  set_I_e(tmp_I_e);
  set_r(tmp_r);
  set_V_m(tmp_V_m);
  set_Ca_i(tmp_Ca_i);
  set_alpha_n_init(tmp_alpha_n_init);
  set_beta_n_init(tmp_beta_n_init);
  set_alpha_m_init(tmp_alpha_m_init);
  set_beta_m_init(tmp_beta_m_init);
  set_alpha_h_init(tmp_alpha_h_init);
  set_beta_h_init(tmp_beta_h_init);
  set_Act_m(tmp_Act_m);
  set_Inact_h(tmp_Inact_h);
  set_Act_n(tmp_Act_n);
  set_h_Ca(tmp_h_Ca);
  set_m_Ca(tmp_m_Ca);
  set_m_AHP(tmp_m_AHP);
  set_I_syn_ex__X__spikeExc(tmp_I_syn_ex__X__spikeExc);
  set_I_syn_ex__X__spikeExc__d(tmp_I_syn_ex__X__spikeExc__d);
  set_I_syn_in__X__spikeInh(tmp_I_syn_in__X__spikeInh);
  set_I_syn_in__X__spikeInh__d(tmp_I_syn_in__X__spikeInh__d);


  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
};

#endif /* #ifndef HHCA_PSC_ALPHA */
