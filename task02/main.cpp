#include <filesystem>
// #include <experimental/filesystem> // uncomment here if the <filesystem> cannot be included above
//
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Eigen/Core"
//
#include "parse_svg.h"

/***
 * signed area of a triangle connecting points (p0, p1, p2) in counter-clockwise order.
 * @param p0 1st point xy-coordinate
 * @param p1 2nd point xy-coordinate
 * @param p2 3rd point xy-coordinate
 * @return signed area (float)
 */
float area(
    const Eigen::Vector2f &p0,
    const Eigen::Vector2f &p1,
    const Eigen::Vector2f &p2) {
  const auto v01 = p1 - p0;
  const auto v02 = p2 - p0;
  // return 0.5f * (v01[0] * v02[1] - v01[1] * v02[0]); // right handed coordinate
  return 0.5f * (v01[1] * v02[0] - v01[0] * v02[1]); // left-handed coordinate (because pixel y-coordinate is going down)
}

/***
 * compute number of intersection of a ray against a line segment
 * @param org ray origin
 * @param dir ray direction (unit normal)
 * @param ps one of the two end points
 * @param pe the other end point
 * @return number of intersection
 */
int number_of_intersection_ray_against_edge(
    const Eigen::Vector2f &org,
    const Eigen::Vector2f &dir,
    const Eigen::Vector2f &ps,
    const Eigen::Vector2f &pe) {
  auto a = area(org, org + dir, ps);
  auto b = area(org, pe, org + dir);
  auto c = area(org, ps, pe);
  auto d = area(dir+ps, ps, pe);
  if (a * b > 0.f && d * c < 0.f) { return 1; }
  return 0;
  // the following code was a bug
  //auto d = area(org + dir, ps, pe);
  //if (a * b > 0.f && d * c > 0.f && fabs(d) > fabs(c)) { return 1; }
}

Eigen::VectorXf remove_zero_in_high_order(Eigen::VectorXf &p) {
  int size = p.size();
  while (size > 1 && p[size - 1] == 0) {
    --size;
  }
  p.conservativeResize(size);
  return p;
}

/* calculate the coefficients of Linear Interpolation and return the
 * coefficients. p[0] is the coefficient for the constant term, p[1] is the
 * coefficient for the linear term, and so on.
 * */
Eigen::VectorXf lerp(Eigen::VectorXf &p0, Eigen::VectorXf &p1) {
  Eigen::VectorXf p0_term = Eigen::VectorXf::Zero(p0.size() + 1);
  p0_term[0] = p0[0];
  for (int i = 1; i < p0.size(); i++) {
    p0_term[i] = p0[i] - p0[i - 1];
  }
  p0_term[p0.size()] = -p0[p0.size() - 1];

  Eigen::VectorXf p1_term = Eigen::VectorXf::Zero(p1.size() + 1);
  for (int i = 0; i < p1.size(); i++) {
    p1_term[i + 1] = p1[i];
  }

  int max_order = std::max(p0_term.size(), p1_term.size());
  Eigen::VectorXf p = Eigen::VectorXf::Zero(max_order);
  for (int i = 0; i < max_order; i++) {
    if (i < p0_term.size()) {
      p[i] += p0_term[i];
    }
    if (i < p1_term.size()) {
      p[i] += p1_term[i];
    }
  }
  return p;
}

float calculate_polynomial(const Eigen::VectorXf &p, float t) {
  float result = 0;
  // Start from the highest degree coefficient
  for (int i = p.size() - 1; i >= 0; --i) {
    result = result * t + p[i];  // Horner's method
  }
  return result;
}

/***
 *
 * @param org ray origin
 * @param dir ray direction (unit vector)
 * @param ps one of the two end points
 * @param pc control point
 * @param pe the other end point
 * @return the number of intersections
 */
int number_of_intersection_ray_against_quadratic_bezier(
    const Eigen::Vector2f &org,
    const Eigen::Vector2f &dir,
    const Eigen::Vector2f &ps,
    const Eigen::Vector2f &pc,
    const Eigen::Vector2f &pe) {
  // comment out below to do the assignment
//  return number_of_intersection_ray_against_edge(org, dir, ps, pe);
  // write some code below to find the intersection between ray and the quadratic

  /* construct the equation with t */
  // get p(t)
  Eigen::VectorXf ps0 = Eigen::VectorXf::Zero(1);
  ps0[0] = ps[0];
  Eigen::VectorXf pc0 = Eigen::VectorXf::Zero(1);
  pc0[0] = pc[0];
  Eigen::VectorXf pe0 = Eigen::VectorXf::Zero(1);
  pe0[0] = pe[0];
  Eigen::VectorXf a0 = lerp(ps0, pc0);
  Eigen::VectorXf b0 = lerp(pc0, pe0);
  Eigen::VectorXf p0 = lerp(a0, b0);

  Eigen::VectorXf ps1 = Eigen::VectorXf::Zero(1);
  ps1[0] = ps[1];
  Eigen::VectorXf pc1 = Eigen::VectorXf::Zero(1);
  pc1[0] = pc[1];
  Eigen::VectorXf pe1 = Eigen::VectorXf::Zero(1);
  pe1[0] = pe[1];
  Eigen::VectorXf a1 = lerp(ps1, pc1);
  Eigen::VectorXf b1 = lerp(pc1, pe1);
  Eigen::VectorXf p1 = lerp(a1, b1);

  assert((p0.size() == p1.size()) && "p0 must be the same dimension with p1");

  Eigen::MatrixXf p(2, p0.size());
  p << p0.transpose(), p1.transpose();

  // get w, where w^T * v =0
  Eigen::Vector2f dir_perp = Eigen::Vector2f(-dir[1], dir[0]);

  //  get w^T * p(t) - w^T * q = 0
  Eigen::VectorXf lhe = dir_perp.transpose() * p;
  lhe[0] += -dir_perp.dot(org);

  /* find the root for the equation */
  /* Since the equation is quadratic, we can get the analytic solution.
   * Otherwise, we need to implement the newton method or other numeric methods to
   * find the solution. */
  float a = lhe[2];
  float b = lhe[1];
  float c = lhe[0];
  float delta = b * b - 4 * a * c;

  if (delta < 0) {
    return 0; // no solution
  }

  float t_0 = (-b + sqrt(delta)) / (2 * a);
  float t_1 = (-b - sqrt(delta)) / (2 * a);

  /* calculate s and check the solution: 0<=t<=1 and s>0 */
  int num_solution = 0;

  if (t_0 > 0 && t_0 < 1) {
    float s_0 = (calculate_polynomial(p0, t_0) - org[0]) / dir[0];
    float s_0_other = (calculate_polynomial(p1, t_0) - org[1]) / dir[1];
    //  std::cout << " s_0 " << s_0 << " s_0_other " << s_0_other << std::endl;
    assert((s_0 - s_0_other < 1e-3) &&
           "something wrong with root finding, please check!");
    if (s_0 > 0) {
      num_solution++;
    }
  }

  if (delta == 0) {
    return num_solution; // only one solution
  }

  if (t_1 > 0 && t_1 < 1) {
    float s_1 = (calculate_polynomial(p0, t_1) - org[0]) / dir[0];
    float s_1_other = (calculate_polynomial(p1, t_1) - org[1]) / dir[1];
    //  std::cout << " s_1 " << s_1 << " s_1_other " << s_1_other << std::endl;
    assert((s_1 - s_1_other < 1e-3) &&
           "something wrong with root finding, please check!");
    if (s_1 > 0) {
      num_solution++;
    }
  }

  return num_solution;
}

int main() {
  const auto input_file_path = std::filesystem::path(PROJECT_SOURCE_DIR) / ".." / "asset" / "r.svg";
  const auto [width, height, shape] = acg::svg_get_image_size_and_shape(input_file_path);
  if (width == 0) { // something went wrong in loading the function
    std::cout << "file open failure" << std::endl;
    abort();
  }
  const std::vector<std::string> outline_path = acg::svg_outline_path_from_shape(shape);
  const std::vector<std::vector<acg::Edge>> loops = acg::svg_loops_from_outline_path(outline_path);
  //
  std::vector<unsigned char> img_data(width * height, 255); // grayscale image initialized white
  for (unsigned int ih = 0; ih < height; ++ih) {
    for (unsigned int iw = 0; iw < width; ++iw) {
      const auto org = Eigen::Vector2f(iw + 0.5, ih + 0.5); // pixel center
      const auto dir = Eigen::Vector2f(60., 20.); // search direction
      int count_cross = 0;
      for (const auto &loop: loops) { // loop over loop (letter R have internal/external loops)
        for (const auto &edge: loop) { // loop over edge in the loop
          if (edge.is_bezier) { // in case the edge is a quadratic Bézier
            count_cross += number_of_intersection_ray_against_quadratic_bezier(
                org, dir,
                edge.ps, edge.pc, edge.pe);
          } else { // in case the edge is a line segment
            count_cross += number_of_intersection_ray_against_edge(
                org, dir,
                edge.ps, edge.pe);
          }
        }
      }
      if (count_cross % 2 == 1) { // Jordan's curve theory
        img_data[ih * width + iw] = 0; // paint black if it is inside
      }
    }
  }
  const auto output_file_path = std::filesystem::path(PROJECT_SOURCE_DIR) / "output.png";
  stbi_write_png(output_file_path.string().c_str(), width, height, 1, img_data.data(), width);
}
