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
    T result = 0;
    if ((rho_ < (T)(fabs(eta_ - r)))) result += 0.;
    else if ((rho_ > (eta_ + r))) result += M_PI;
    else result += (acos((eta_ * eta_ + r * r - rho_ * rho_) / (2 * eta_ * r)));
    return result;
  }

private:
  T rho_;
  T eta_;
};

template <typename T>
class intersection_volume_unit_sphere {
public:
  intersection_volume_unit_sphere(T rho, T eta) {
    f.set_variables(rho, eta);
  }

  T operator()(T r, T v) {
    return r * (T)sqrt(1 - r * r) * f(r, v);
  }

  T integrate(T dr) {
    double result = 0.;
    for (double r = 0; r <= 1 - dr; r += dr) {
      result += runge_kutta_4_step(r, dr, 0., *this);
    }
    return 4.0 * result;
  }

 public:
  volume_multiplier_function<T> f;
};



int main()
{
  intersection_volume_unit_sphere<double> f(0.05, 0.025);
  double dr = 0.0001;

  std::cout << f.integrate(dr) << std::endl;

}

