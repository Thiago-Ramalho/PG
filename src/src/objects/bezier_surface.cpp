#include "Prism/objects/bezier_surface.hpp"

#include "Prism/core/matrix.hpp"
#include "Prism/core/utils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Prism {

BezierSurface::BezierSurface(const std::vector<std::vector<Point3>>& control_points, 
                             std::shared_ptr<Material> material, int u_segments, int v_segments)
    : control_points(control_points), material(std::move(material)), 
      u_segments(u_segments), v_segments(v_segments) {
    if (control_points.empty() || control_points[0].empty()) {
        throw std::invalid_argument("Control points grid cannot be empty");
    }
    
    v_count = static_cast<int>(control_points.size()); //Number of control points in v direction
    u_count = static_cast<int>(control_points[0].size()); // Number of control points in u direction
    
    // Verify that all rows have the same number of control points
    for (const auto& row : control_points) {
        if (static_cast<int>(row.size()) != u_count) {
            throw std::invalid_argument("All rows must have the same number of control points");
        }
    }
    
    // Tessellate the surface into triangles
    tessellateSurface();
    
    // Compute bounding box
    computeBoundingBox();
}

bool BezierSurface::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
    Ray transformed_ray = ray.transform(inverseTransform);
    
    // First, test ray against bounding box for early rejection
    if (!rayIntersectsBoundingBox(transformed_ray, t_min, t_max)) {
        return false;
    }
    
    // Test intersection against all triangles in the mesh
    HitRecord closest_hit;
    closest_hit.t = t_max;
    bool hit_found = false;
    
    for (const auto& triangle : triangulated_mesh) {
        HitRecord temp_rec;
        if (triangle.hit(transformed_ray, t_min, closest_hit.t, temp_rec)) {
            closest_hit = temp_rec;
            hit_found = true;
        }
    }
    
    if (!hit_found) {
        return false;
    }
    
    // Transform hit point back to world space
    Point3 world_hit = transform * transformed_ray.at(closest_hit.t);
    double world_t = (world_hit - ray.origin()).magnitude();
    
    if (world_t < t_min || world_t > t_max) {
        return false;
    }
    
    rec.t = world_t;
    rec.p = world_hit;
    
    // Transform normal to world space
    Vector3 world_normal = (inverseTransposeTransform * closest_hit.normal).normalize();
    rec.set_face_normal(ray, world_normal);
    rec.material = material;
    
    return true;
}

Point3 BezierSurface::evaluateAt(double u, double v) const {
    // Clamp parameters to valid range instead of throwing
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    Vector3 result = Vector3(0, 0, 0);
    int n = u_count - 1; // degree in u direction
    int m = v_count - 1; // degree in v direction
    
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            double basis_u = bernstein(i, n, u);
            double basis_v = bernstein(j, m, v);
            double weight = basis_u * basis_v;
            
            // Convert Point3 to Vector3, multiply by weight, then back to Point3
            Vector3 point_vec(control_points[j][i]);
            Vector3 weighted_vec = point_vec * weight;
            result = result + weighted_vec;
        }
    }
    
    return Point3(result);
}

Vector3 BezierSurface::normalAt(double u, double v) const {
    // Clamp parameters to valid range
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));
    
    const double h = 1e-5; // Small step for numerical differentiation
    
    // Compute partial derivatives using finite differences
    Point3 p_center = evaluateAt(u, v);
    
    // Partial derivative with respect to u
    Vector3 du;
    if (u + h <= 1.0) {
        Point3 p_u_plus = evaluateAt(u + h, v);
        du = p_u_plus - p_center;
    } else {
        Point3 p_u_minus = evaluateAt(u - h, v);
        du = p_center - p_u_minus;
    }
    du = du * (1.0 / h);
    
    // Partial derivative with respect to v
    Vector3 dv;
    if (v + h <= 1.0) {
        Point3 p_v_plus = evaluateAt(u, v + h);
        dv = p_v_plus - p_center;
    } else {
        Point3 p_v_minus = evaluateAt(u, v - h);
        dv = p_center - p_v_minus;
    }
    dv = dv * (1.0 / h);
    
    // Normal is cross product of partial derivatives
    Vector3 normal = du.cross(dv);
    double magnitude = normal.magnitude();
    if (magnitude < 1e-12) {
        // Return a default normal if surface is degenerate at this point
        return Vector3(0, 1, 0);
    }
    return normal.normalize();
}

double BezierSurface::bernstein(int i, int n, double t) const {
    if (i < 0 || i > n || t < 0.0 || t > 1.0) {
        return 0.0;
    }
    
    double coeff = static_cast<double>(binomial(n, i));
    double term1 = std::pow(1.0 - t, n - i);
    double term2 = std::pow(t, i);
    
    return coeff * term1 * term2;
}

int BezierSurface::binomial(int n, int k) const {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use symmetry property: C(n,k) = C(n,n-k)
    if (k > n - k) k = n - k;
    
    int result = 1;
    for (int i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    
    return result;
}

void BezierSurface::tessellateSurface() {
    triangulated_mesh.clear();
    
    // Generate points on the surface
    std::vector<std::vector<Point3>> surface_points(v_segments + 1, std::vector<Point3>(u_segments + 1));
    std::vector<std::vector<Vector3>> surface_normals(v_segments + 1, std::vector<Vector3>(u_segments + 1));
    
    for (int i = 0; i <= v_segments; ++i) {
        for (int j = 0; j <= u_segments; ++j) {
            double u = static_cast<double>(j) / u_segments;
            double v = static_cast<double>(i) / v_segments;
            
            surface_points[i][j] = evaluateAt(u, v);
            surface_normals[i][j] = normalAt(u, v);
        }
    }
    
    // Create triangles from the grid
    for (int i = 0; i < v_segments; ++i) {
        for (int j = 0; j < u_segments; ++j) {
            // Get the four corners of the current quad
            Point3 p00 = surface_points[i][j];
            Point3 p10 = surface_points[i][j + 1];
            Point3 p01 = surface_points[i + 1][j];
            Point3 p11 = surface_points[i + 1][j + 1];
            
            Vector3 n00 = surface_normals[i][j];
            Vector3 n10 = surface_normals[i][j + 1];
            Vector3 n01 = surface_normals[i + 1][j];
            Vector3 n11 = surface_normals[i + 1][j + 1];
            
            // Create two triangles per quad
            // Triangle 1: p00, p10, p01
            triangulated_mesh.emplace_back(
                std::make_shared<Point3>(p00),
                std::make_shared<Point3>(p10),
                std::make_shared<Point3>(p01),
                std::make_shared<Vector3>(n00),
                std::make_shared<Vector3>(n10),
                std::make_shared<Vector3>(n01)
            );
            
            // Triangle 2: p10, p11, p01
            triangulated_mesh.emplace_back(
                std::make_shared<Point3>(p10),
                std::make_shared<Point3>(p11),
                std::make_shared<Point3>(p01),
                std::make_shared<Vector3>(n10),
                std::make_shared<Vector3>(n11),
                std::make_shared<Vector3>(n01)
            );
        }
    }
}

void BezierSurface::computeBoundingBox() {
    if (triangulated_mesh.empty()) {
        bbox_min = Point3(0, 0, 0);
        bbox_max = Point3(0, 0, 0);
        return;
    }
    
    // Initialize with first point
    Point3 first_point = triangulated_mesh[0].getPoint1();
    bbox_min = first_point;
    bbox_max = first_point;
    
    // Find min/max coordinates from all triangle vertices
    for (const auto& triangle : triangulated_mesh) {
        Point3 points[3] = {
            triangle.getPoint1(),
            triangle.getPoint2(),
            triangle.getPoint3()
        };
        
        for (const Point3& p : points) {
            bbox_min.x = std::min(bbox_min.x, p.x);
            bbox_min.y = std::min(bbox_min.y, p.y);
            bbox_min.z = std::min(bbox_min.z, p.z);
            
            bbox_max.x = std::max(bbox_max.x, p.x);
            bbox_max.y = std::max(bbox_max.y, p.y);
            bbox_max.z = std::max(bbox_max.z, p.z);
        }
    }
}

bool BezierSurface::rayIntersectsBoundingBox(const Ray& ray, double t_min, double t_max) const {
    // Ray-AABB intersection using the slab method
    Vector3 inv_dir(1.0 / ray.direction().x, 1.0 / ray.direction().y, 1.0 / ray.direction().z);
    Vector3 origin = ray.origin();
    
    double t1 = (bbox_min.x - origin.x) * inv_dir.x;
    double t2 = (bbox_max.x - origin.x) * inv_dir.x;
    if (t1 > t2) std::swap(t1, t2);
    
    double tymin = (bbox_min.y - origin.y) * inv_dir.y;
    double tymax = (bbox_max.y - origin.y) * inv_dir.y;
    if (tymin > tymax) std::swap(tymin, tymax);
    
    if (t1 > tymax || tymin > t2) return false;
    
    t1 = std::max(t1, tymin);
    t2 = std::min(t2, tymax);
    
    double tzmin = (bbox_min.z - origin.z) * inv_dir.z;
    double tzmax = (bbox_max.z - origin.z) * inv_dir.z;
    if (tzmin > tzmax) std::swap(tzmin, tzmax);
    
    if (t1 > tzmax || tzmin > t2) return false;
    
    t1 = std::max(t1, tzmin);
    t2 = std::min(t2, tzmax);
    
    return t2 >= t_min && t1 <= t_max;
}

} // namespace Prism