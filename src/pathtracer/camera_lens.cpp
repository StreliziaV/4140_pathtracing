#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Assignment 7: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
  double c_x = 2.0 * tan(M_PI * 0.5 * hFov / 180.0) * x - tan(M_PI * 0.5 * hFov / 180.0);
  double c_y = 2.0 * tan(M_PI * 0.5 * vFov / 180.0) * y - tan(M_PI * 0.5 * vFov / 180.0);
  Vector3D c_target(c_x, c_y, -1);

  Vector3D pLens(this->lensRadius * sqrt(rndR) * cos(rndTheta), this->lensRadius * sqrt(rndR) * sin(rndTheta), 0);

  double f_t = (- this->focalDistance - 0.0) / -1.0;
  Vector3D c_pf = Vector3D(0, 0, 0) + f_t * c_target;
  Vector3D c_blue = c_pf - pLens;
  c_blue.normalize();

  Vector3D w_pLens = c2w * pLens + pos;
  Vector3D w_blue = (c2w * c_blue);
  w_blue.normalize();

  Ray g_ray(w_pLens, w_blue);
  g_ray.max_t = fClip;
  g_ray.min_t = nClip;
  return g_ray;

  // double c_x = 2.0 * tan(M_PI * 0.5 * hFov / 180.0) * x - tan(M_PI * 0.5 * hFov / 180.0);
  // double c_y = 2.0 * tan(M_PI * 0.5 * vFov / 180.0) * y - tan(M_PI * 0.5 * vFov / 180.0);
  // Vector3D c_target(c_x, c_y, -1);
  // Vector3D my_result = c_target;
  // my_result = (c2w * my_result);
  // my_result.normalize();
  // Ray g_ray(pos, my_result);
  // g_ray.max_t = fClip;
  // g_ray.min_t = nClip;
  // return g_ray;
  
  // return Ray(pos, Vector3D(0, 0, -1));
}


} // namespace CGL
