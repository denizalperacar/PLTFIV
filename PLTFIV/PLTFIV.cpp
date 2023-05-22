#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <vector>


template <typename T, typename F>
T runge_kutta_4_step(const T& tn, const T& h, const T& yn, F& der) {
  T half_h = h / 2;
  T k1 = h * der(tn, yn);
  T k2 = h * der(tn + half_h, yn + k1 * half_h);
  T k3 = h * der(tn + half_h, yn + k2 * half_h);
  T k4 = h * der(tn + h, yn + k3 * h);
  return yn + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

template <typename T>
class volume_multiplier_function {

public:
  volume_multiplier_function() = default;
  void set_variables(T rho, T eta) {
    rho_ = rho;
    eta_ = eta;
	}
  volume_multiplier_function(T rho, T eta) : rho_(rho), eta_(eta) {}

  T operator()(T r, T v) {
		return 0 * (rho_ < (T)(fabs(eta_ - r))) 
      + (M_PI) * (rho_ > (eta_ + r))
      + (acos((eta_ * eta_ + r * r - rho_ * rho_) / (2 * eta_ * r))) 
      * (rho_ > (T)fabs(eta_ - r) && rho_ < (eta_ + r));
  }

private:
  T rho_;
  T eta_;
};

template <typename T>
class intersection_volume {
public:
  intersection_volume(T rho, T eta) {
    f.set_variables(rho, eta);
  }

  T operator()(T r, T v) {
    return r * (T)sqrt(1 - r * r) * f(r, v);
  }

 public:
  volume_multiplier_function<T> f;

};



int main()
{
  intersection_volume<double> f(0.05, 0.025);
  double dr = 0.01;
  double result = 0;
  for (double r = 0; r <= 1; r += dr) {
    result = runge_kutta_4_step(r, dr, 0., f);
    std::cout << f.f(r, 0.) << std::endl;
  }
  

}

