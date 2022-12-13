#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Implement MirrorBSDF
  this->reflect(wo, wi);
  *pdf = 1;
  return reflectance / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Assignment 7: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  Vector3D macro_normal(0, 0, 1);
  double cos_theta_h = dot(macro_normal, h);
  double tan_theta_h_2 = (1 - pow(cos_theta_h, 2)) / pow(cos_theta_h, 2);

  double D_h = pow(exp(1), - tan_theta_h_2 / pow(this->alpha, 2)) / (M_PI * pow(this->alpha, 2) * pow(cos_theta_h, 4));
  // if (D_h >= 1) std::cout << D_h << std::endl;
  return D_h;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
  Vector3D my_eta = this->eta;
  Vector3D my_k = this->k;
  double cos_theta_i = wi.z;

  Vector3D helper_term = (my_eta * my_eta) + (my_k * my_k);
  double R_s_x = (helper_term.x - 2 * my_eta.x * cos_theta_i + pow(cos_theta_i, 2)) / (helper_term.x + 2 * my_eta.x * cos_theta_i + pow(cos_theta_i, 2));
  double R_p_x = (helper_term.x * pow(cos_theta_i, 2) - 2 * my_eta.x * cos_theta_i + 1) / (helper_term.x * pow(cos_theta_i, 2) + 2 * my_eta.x * cos_theta_i + 1);
  
  double R_s_y = (helper_term.y - 2 * my_eta.y * cos_theta_i + pow(cos_theta_i, 2)) / (helper_term.y + 2 * my_eta.y * cos_theta_i + pow(cos_theta_i, 2));
  double R_p_y = (helper_term.y * pow(cos_theta_i, 2) - 2 * my_eta.y * cos_theta_i + 1) / (helper_term.y * pow(cos_theta_i, 2) + 2 * my_eta.y * cos_theta_i + 1);
  
  double R_s_z = (helper_term.z - 2 * my_eta.z * cos_theta_i + pow(cos_theta_i, 2)) / (helper_term.z + 2 * my_eta.z * cos_theta_i + pow(cos_theta_i, 2));
  double R_p_z = (helper_term.z * pow(cos_theta_i, 2) - 2 * my_eta.z * cos_theta_i + 1) / (helper_term.z * pow(cos_theta_i, 2) + 2 * my_eta.z * cos_theta_i + 1);
  Vector3D Fres_coe((R_s_x + R_p_x) / 2.0, (R_s_y + R_p_y) / 2.0, (R_s_z + R_p_z) / 2.0);
  return Fres_coe;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Implement microfacet model here.
  if (wo.z <= 0 || wi.z <= 0) return Vector3D(0, 0, 0);
  Vector3D my_h = (wo + wi) / 2.0;
  // std::cout << my_h.norm() << ' ' << wo.norm() << ' ' << wi.norm() << std::endl;
  my_h.normalize();
  Vector3D macro_normal(0, 0, 1);
  Vector3D my_f = F(wi) * G(wo, wi) * D(my_h) / (4 * dot(macro_normal, wo) * dot(macro_normal, wi));
  return my_f;
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  Vector2D my_r = this->sampler.get_sample();
  double theta_h = atan(sqrt(- pow(this->alpha, 2) * log(1 - my_r.x)));
  double phi_h = 2 * M_PI * my_r.y;
  Vector3D my_h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
  my_h.normalize();
  *wi = 2 * my_h - wo;
  (*wi).normalize();
  if ((*wi).z <= 0) {
    *pdf = 0.0;
    return Vector3D(0, 0, 0);
  }

  double p_theta_h = 2 * sin(theta_h) * exp(- pow(tan(theta_h), 2)/ pow(this->alpha, 2)) / (pow(this->alpha, 2) * pow(cos(theta_h), 3));
  double p_phi_h = 1 / (2 * M_PI);
  double p_h = (p_theta_h * p_phi_h) / sin(theta_h);
  *pdf = p_h / (4 * dot(*wi, my_h));

  // *wi = cosineHemisphereSampler.get_sample(pdf);
  Vector3D f_result = MicrofacetBSDF::f(wo, *wi);
  return f_result;
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement RefractionBSDF
  *pdf = 1;
  if (!this->refract(wo, wi, this->ior)) return Vector3D();

  double eta = 0.0;
  if (wo.z >= 0) eta = 1.0 / ior;
  else eta = ior / 1.0;
  return this->transmittance / abs_cos_theta(*wi) / pow(eta, 2);
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

  double eta = 0.0;
  int dir = 1;
  bool if_full_refl = false;
  if (wo.z >= 0) {
    eta = 1.0 / this->ior;
    dir = -1;
  }
  else eta = this->ior / 1.0;
  double my_test = 1 - pow(eta, 2) * (1 - pow(wo.z, 2));
  if (my_test <= 0) if_full_refl = true;
  else *wi = Vector3D(- eta * wo.x, - eta * wo.y, dir * sqrt(my_test));
  (*wi).normalize();

  if (if_full_refl) {
    *pdf = 1;
    Vector3D refl(-wo.x, -wo.y, wo.z);
    *wi = refl;
    return reflectance / abs_cos_theta(*wi);
  }

  double n1 = 1.0;
  double n2 = this->ior;
  if (wo.z >= 0) {
    n1 = this->ior;
    n2 = 1.0;
  }
  double r0 = pow(((n1 - n2) / (n1 + n2)), 2);
  double Fres_coe = r0 + (1 - r0) * pow((1 - abs_cos_theta(*wi)), 5);

  if (coin_flip(Fres_coe)) {
    Vector3D refl(-wo.x, -wo.y, wo.z);
    *wi = refl;
    *pdf = Fres_coe;
    return Fres_coe * this->reflectance / abs_cos_theta(*wi);
  }

  *pdf = 1 - Fres_coe;
  return (1 - Fres_coe) * this->transmittance / abs_cos_theta(*wi) / pow(eta, 2);
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Assignment 7: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  Vector3D refl(-wo.x, -wo.y, wo.z);
  *wi = refl;
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Assignment 7: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta = 0.0;
  int dir = 1;
  if (wo.z >= 0) {
    eta = 1.0 / ior;
    dir = -1;
  }
  else eta = ior / 1.0;
  double my_test = 1 - pow(eta, 2) * (1 - pow(wo.z, 2));
  if (my_test < 0) return false;

  *wi = Vector3D(- eta * wo.x, - eta * wo.y, dir * sqrt(my_test));
  return true;

}

} // namespace CGL
