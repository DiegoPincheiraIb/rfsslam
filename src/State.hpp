// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

/**
 * \class State
 * \brief An abstract base class for defining a state
 * \tparam State container
 * \author Keith Leung
 */
template<class StateType>
class State
{

public:

  /** Default constructor */
  State(){}

  /** Default destructor */
  ~State(){}

  /** 
   * Function for setting the pose state
   * \param x state to be set
   */
  void set( StateType x ){x_ = x;}
  
  /** 
   * Function for getting the pose state 
   * \param x state [overwritten]
   */
  void get( StateType &x ){x = x_;}

protected:

  StateType x_;

};

/**
 * \class StateWithUncertainty
 * \brief An abstract base class for defining the vehicle pose state
 * \tparam State container for the state
 * \author Keith Leung
 */
template<class StateType, class UncertaintyType>
class StateWithUncertainty : public State<StateType>
{

public:

  /** Default constructor */
  StateWithUncertainty(){};

  /** Default destructor */
  ~StateWithUncertainty(){};

  /** 
   * Function for setting the pose state with uncertainty
   * \param x state to be set
   * \param Sx uncertainty to be set
   */
  void set( StateType x, UncertaintyType Sx){
    this->x_ = x;
    Sx_ = Sx;
  }
  
  /** 
   * Function for getting the pose state with uncertianty
   * \param x state [overwritten]
   * \param Sx uncertainty [overwritten]
   */
  void get( StateType &x, UncertaintyType &Sx){
    x = this->x_;
    Sx = Sx_;
  }

protected:

  UncertaintyType Sx_;

};

#endif
